// la-dml — local-ancestry-aware DSS beta-binomial association per CpG.
//
// For each CpG and ancestral population j, fits a DSS beta-binomial model
// (arcsin transform + iterative WLS, as in LEMUR dml) on ancestry-specific
// counts/coverage:
//
//   arcsin((k+0.1)/(n+0.2)) ~ phenotype + hapcount_j + covariates
//
// Inputs are the per-ancestry matrices from tracts:
//   <prefix>.anc<j>.counts.txt
//   <prefix>.anc<j>.cov.txt
//   <prefix>.anc<j>.hapcount.txt
//
// Usage example:
//   la-dml --tract-prefix chr1 --metadata samples.csv
//       --num-ancs 3 --output chr1.la_dml.tsv

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr const char* kVersion = "1.0.1";
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kArcsinC0 = 0.1;
constexpr double kPhiMin = 0.001;
constexpr double kPhiMax = 0.999;

constexpr int kPhenotypeIdx = 1;
constexpr int kHapcountIdx = 2;

struct Args {
    std::string tract_prefix;
    std::string metadata;
    int num_ancs = 0;
    std::string output;
    std::string case_label = "Case";
    std::string control_label = "Control";
    int num_threads = 1;
    std::size_t batch_size = 8192;
};

struct SampleMeta {
    std::string sample_id;
    int phenotype_bin = 0;
    int col_idx = -1;
    std::vector<double> covariates;
};

struct CohortMeta {
    std::vector<SampleMeta> samples;
    int n_params = 0;
};

struct SiteRow {
    std::string chr;
    int pos = 0;
    std::vector<std::vector<double>> counts;
    std::vector<std::vector<double>> cov;
    std::vector<std::vector<double>> hapcount;
};

struct CoefResult {
    double beta = kNaN;
    double se = kNaN;
    double pvalue = kNaN;
    double tval = kNaN;
};

struct AncSiteResult {
    int n = 0;
    double mean_frac = kNaN;
    double la_prop = kNaN;
    double phi = kNaN;
    CoefResult phenotype;
    CoefResult hapcount;
};

struct SiteResult {
    std::string chr;
    int pos = 0;
    std::vector<AncSiteResult> anc;
};

struct Obs {
    double k;
    double n;
    double hapcount;
    std::vector<double> x;
};

bool is_nan(double x) { return std::isnan(x); }

void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog
        << " --tract-prefix <prefix> --metadata <samples.csv>\n"
        << "    --num-ancs <N> --output <results.tsv>\n"
        << "    [--case-label Case] [--control-label Control]\n"
        << "    [-j N] [-b BATCH]\n\n"
        << "Fits per-ancestry DSS beta-binomial models on tracts output.\n";
}

std::vector<std::string> split_line(const std::string& line, char delim) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream iss(line);
    while (std::getline(iss, field, delim)) fields.push_back(field);
    return fields;
}

std::string trim(const std::string& s) {
    std::size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) ++start;
    std::size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) --end;
    return s.substr(start, end - start);
}

int column_index(const std::vector<std::string>& header, const std::string& name) {
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name) return static_cast<int>(i);
    }
    return -1;
}

double parse_field(const std::string& s) {
    if (s.empty() || s == ".") return kNaN;
    try {
        return std::stod(s);
    } catch (...) {
        return kNaN;
    }
}

bool file_exists(const std::string& path) {
    std::ifstream in(path);
    return static_cast<bool>(in);
}

double normal_sf(double z) {
    const double abs_z = std::abs(z);
    const double t = 1.0 / (1.0 + 0.2316419 * abs_z);
    const double d = 0.3989423 * std::exp(-0.5 * abs_z * abs_z);
    const double p = d * t *
                     (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
    return z >= 0 ? p : 1.0 - p;
}

int phenotype_bin_from_field(const std::string& phen, const std::string& case_label,
                             const std::string& control_label) {
    if (phen == case_label || phen == "2") return 1;
    if (phen == control_label || phen == "1") return 0;
    throw std::runtime_error("Unexpected phenotype value: " + phen);
}

Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need_value = [&](const char* flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(std::string("Missing value for ") + flag);
            return argv[++i];
        };

        if (arg == "--tract-prefix") {
            args.tract_prefix = need_value("--tract-prefix");
        } else if (arg == "--metadata") {
            args.metadata = need_value("--metadata");
        } else if (arg == "--num-ancs") {
            args.num_ancs = std::stoi(need_value("--num-ancs"));
        } else if (arg == "--output") {
            args.output = need_value("--output");
        } else if (arg == "--case-label") {
            args.case_label = need_value("--case-label");
        } else if (arg == "--control-label") {
            args.control_label = need_value("--control-label");
        } else if (arg == "-j") {
            args.num_threads = std::stoi(need_value("-j"));
        } else if (arg == "-b") {
            args.batch_size = static_cast<std::size_t>(std::stoul(need_value("-b")));
        } else if (arg == "-h" || arg == "--help") {
            usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (args.tract_prefix.empty() || args.metadata.empty() || args.output.empty() ||
        args.num_ancs <= 0) {
        usage(argv[0]);
        throw std::runtime_error("Required: --tract-prefix, --metadata, --num-ancs, --output");
    }
    if (!file_exists(args.metadata)) {
        throw std::runtime_error("Metadata file not found: " + args.metadata);
    }
    for (int anc = 0; anc < args.num_ancs; ++anc) {
        const std::string suffix = ".anc" + std::to_string(anc);
        if (!file_exists(args.tract_prefix + suffix + ".counts.txt")) {
            throw std::runtime_error("Missing file: " + args.tract_prefix + suffix + ".counts.txt");
        }
        if (!file_exists(args.tract_prefix + suffix + ".cov.txt")) {
            throw std::runtime_error("Missing file: " + args.tract_prefix + suffix + ".cov.txt");
        }
        if (!file_exists(args.tract_prefix + suffix + ".hapcount.txt")) {
            throw std::runtime_error("Missing file: " + args.tract_prefix + suffix + ".hapcount.txt");
        }
    }
    return args;
}

CohortMeta load_metadata(const std::string& path, const std::vector<std::string>& tract_header,
                         const std::string& case_label, const std::string& control_label) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open metadata: " + path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty metadata: " + path);

    const auto header = split_line(line, ',');
    const int id_idx = column_index(header, "sample_id");
    const int phen_idx = column_index(header, "phenotype");
    const int age_idx = column_index(header, "AGE");
    const int sex_idx = column_index(header, "SEX");
    const int bmi_idx = column_index(header, "BMI");
    const int amr_idx = column_index(header, "AMR");
    const int cov_mean_idx = column_index(header, "COVERAGE_MEAN");

    if (id_idx < 0 || phen_idx < 0 || age_idx < 0 || bmi_idx < 0) {
        throw std::runtime_error("Metadata must contain sample_id,phenotype,AGE,BMI");
    }

    const bool ancestry_design = amr_idx >= 0 && cov_mean_idx >= 0 && sex_idx >= 0;

    CohortMeta meta;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto fields = split_line(line, ',');
        if (static_cast<int>(fields.size()) <= std::max({id_idx, phen_idx, age_idx, bmi_idx})) continue;

        SampleMeta s;
        s.sample_id = trim(fields[static_cast<std::size_t>(id_idx)]);
        s.phenotype_bin =
            phenotype_bin_from_field(trim(fields[static_cast<std::size_t>(phen_idx)]), case_label, control_label);
        s.col_idx = column_index(tract_header, s.sample_id);
        if (s.col_idx < 0) {
            throw std::runtime_error("Sample " + s.sample_id + " missing from tract matrices");
        }

        const double age = std::stod(trim(fields[static_cast<std::size_t>(age_idx)]));
        const double bmi = std::stod(trim(fields[static_cast<std::size_t>(bmi_idx)]));
        if (ancestry_design) {
            const double sex = std::stod(trim(fields[static_cast<std::size_t>(sex_idx)]));
            const double coverage = std::stod(trim(fields[static_cast<std::size_t>(cov_mean_idx)]));
            const double amr = std::stod(trim(fields[static_cast<std::size_t>(amr_idx)]));
            s.covariates = {1.0, static_cast<double>(s.phenotype_bin), 0.0, age, sex, bmi, amr, coverage};
        } else {
            s.covariates = {1.0, static_cast<double>(s.phenotype_bin), 0.0, age, bmi};
        }
        meta.samples.push_back(std::move(s));
    }

    if (meta.samples.empty()) throw std::runtime_error("No samples loaded from metadata");
    meta.n_params = static_cast<int>(meta.samples.front().covariates.size());
    return meta;
}

double arcsin_methylation(double y, double n) {
    const double frac = (y + kArcsinC0) / (n + 2.0 * kArcsinC0);
    return std::asin(2.0 * frac - 1.0);
}

double predict(const Obs& o, const std::vector<double>& beta) {
    double eta = 0.0;
    for (std::size_t j = 0; j < beta.size(); ++j) eta += beta[j] * o.x[j];
    return eta;
}

bool cholesky_solve(std::vector<double>& a, int p, std::vector<double>& b) {
    std::vector<double> L(static_cast<std::size_t>(p * p), 0.0);
    auto at = [&](int r, int c) -> double& { return a[static_cast<std::size_t>(r * p + c)]; };
    auto lt = [&](int r, int c) -> double& { return L[static_cast<std::size_t>(r * p + c)]; };

    for (int i = 0; i < p; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = at(i, j);
            for (int k = 0; k < j; ++k) sum -= lt(i, k) * lt(j, k);
            if (i == j) {
                if (sum <= 1e-12) return false;
                lt(i, j) = std::sqrt(sum);
            } else {
                lt(i, j) = sum / lt(j, j);
            }
        }
    }

    for (int i = 0; i < p; ++i) {
        double sum = b[static_cast<std::size_t>(i)];
        for (int k = 0; k < i; ++k) sum -= lt(i, k) * b[static_cast<std::size_t>(k)];
        b[static_cast<std::size_t>(i)] = sum / lt(i, i);
    }
    for (int i = p - 1; i >= 0; --i) {
        double sum = b[static_cast<std::size_t>(i)];
        for (int k = i + 1; k < p; ++k) sum -= lt(k, i) * b[static_cast<std::size_t>(k)];
        b[static_cast<std::size_t>(i)] = sum / lt(i, i);
    }
    return true;
}

bool invert_spd(std::vector<double> a, int p, std::vector<double>& inv) {
    inv.assign(static_cast<std::size_t>(p * p), 0.0);
    std::vector<double> rhs(static_cast<std::size_t>(p), 0.0);
    std::vector<double> col(static_cast<std::size_t>(p), 0.0);
    for (int j = 0; j < p; ++j) {
        std::fill(rhs.begin(), rhs.end(), 0.0);
        rhs[static_cast<std::size_t>(j)] = 1.0;
        col = rhs;
        std::vector<double> work = a;
        if (!cholesky_solve(work, p, col)) return false;
        for (int i = 0; i < p; ++i) inv[static_cast<std::size_t>(i * p + j)] = col[static_cast<std::size_t>(i)];
    }
    return true;
}

bool wls(const std::vector<Obs>& obs, int p, const std::vector<double>& z,
         const std::vector<double>& weights, std::vector<double>& beta) {
    std::vector<double> xtx(static_cast<std::size_t>(p * p), 0.0);
    beta.assign(static_cast<std::size_t>(p), 0.0);
    std::vector<double> xtz(static_cast<std::size_t>(p), 0.0);

    for (std::size_t i = 0; i < obs.size(); ++i) {
        const double w = weights[i];
        const double zi = z[i];
        for (int j = 0; j < p; ++j) {
            const double xj = obs[i].x[static_cast<std::size_t>(j)];
            xtz[static_cast<std::size_t>(j)] += w * xj * zi;
            for (int k = 0; k < p; ++k) {
                xtx[static_cast<std::size_t>(j * p + k)] +=
                    w * xj * obs[i].x[static_cast<std::size_t>(k)];
            }
        }
    }
    beta = xtz;
    return cholesky_solve(xtx, p, beta);
}

std::optional<std::pair<std::vector<double>, double>> dss_fit(
    const std::vector<Obs>& obs, int p, std::vector<double>& var_out) {
    const int n = static_cast<int>(obs.size());
    if (n <= p + 1) return std::nullopt;

    std::vector<double> z(static_cast<std::size_t>(n));
    std::vector<double> N(static_cast<std::size_t>(n));
    std::vector<double> w1(static_cast<std::size_t>(n));

    for (int i = 0; i < n; ++i) {
        N[static_cast<std::size_t>(i)] = obs[static_cast<std::size_t>(i)].n;
        z[static_cast<std::size_t>(i)] = arcsin_methylation(obs[static_cast<std::size_t>(i)].k, N[static_cast<std::size_t>(i)]);
        w1[static_cast<std::size_t>(i)] = N[static_cast<std::size_t>(i)];
    }

    std::vector<double> beta;
    if (!wls(obs, p, z, w1, beta)) return std::nullopt;

    double sum_sq = 0.0;
    double sum_nm1 = 0.0;
    for (int i = 0; i < n; ++i) {
        const double resid = z[static_cast<std::size_t>(i)] - predict(obs[static_cast<std::size_t>(i)], beta);
        sum_sq += resid * resid * N[static_cast<std::size_t>(i)];
        sum_nm1 += N[static_cast<std::size_t>(i)] - 1.0;
    }
    if (sum_nm1 <= 0.0) return std::nullopt;

    const double df = static_cast<double>(n - p);
    double phi = (sum_sq - df) * static_cast<double>(n) / df / sum_nm1;
    phi = std::clamp(phi, kPhiMin, kPhiMax);

    std::vector<double> w2(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double ni = N[static_cast<std::size_t>(i)];
        w2[static_cast<std::size_t>(i)] = ni / (1.0 + (ni - 1.0) * phi);
    }

    std::vector<double> xtx(static_cast<std::size_t>(p * p), 0.0);
    std::vector<double> xtz(static_cast<std::size_t>(p), 0.0);
    for (std::size_t i = 0; i < obs.size(); ++i) {
        const double w = w2[i];
        const double zi = z[i];
        for (int j = 0; j < p; ++j) {
            const double xj = obs[i].x[static_cast<std::size_t>(j)];
            xtz[static_cast<std::size_t>(j)] += w * xj * zi;
            for (int k = 0; k < p; ++k) {
                xtx[static_cast<std::size_t>(j * p + k)] +=
                    w * xj * obs[i].x[static_cast<std::size_t>(k)];
            }
        }
    }

    beta = xtz;
    if (!cholesky_solve(xtx, p, beta)) return std::nullopt;
    if (!invert_spd(xtx, p, var_out)) return std::nullopt;

    return std::make_pair(beta, phi);
}

CoefResult coef_test(const std::vector<double>& beta, const std::vector<double>& var, int p,
                     int coef_idx) {
    CoefResult out;
    const double var_i = var[static_cast<std::size_t>(coef_idx * p + coef_idx)];
    if (var_i <= 1e-12 || !std::isfinite(var_i)) return out;
    out.beta = beta[static_cast<std::size_t>(coef_idx)];
    out.se = std::sqrt(var_i);
    if (out.se <= 0.0 || !std::isfinite(out.beta)) return out;
    out.tval = out.beta / out.se;
    out.pvalue = 2.0 * normal_sf(std::abs(out.tval));
    return out;
}

std::optional<AncSiteResult> fit_ancestry_site(const SiteRow& site, int anc_idx,
                                               const CohortMeta& meta) {
    std::vector<Obs> obs;
    obs.reserve(meta.samples.size());

    double sum_k = 0.0;
    double sum_n = 0.0;
    double sum_la = 0.0;
    int case_n = 0;
    int control_n = 0;

    for (std::size_t i = 0; i < meta.samples.size(); ++i) {
        const auto& sample = meta.samples[i];
        const double k = site.counts[anc_idx][i];
        const double n = site.cov[anc_idx][i];
        const double la = site.hapcount[anc_idx][i];
        if (is_nan(k) || is_nan(n) || n <= 0.0 || is_nan(la)) continue;

        Obs o;
        o.k = k;
        o.n = n;
        o.hapcount = la;
        o.x = sample.covariates;
        o.x[kHapcountIdx] = la;
        obs.push_back(std::move(o));

        sum_k += k;
        sum_n += n;
        sum_la += la;
        if (sample.phenotype_bin == 1) ++case_n;
        else ++control_n;
    }

    if (static_cast<int>(obs.size()) <= meta.n_params + 1 || case_n == 0 || control_n == 0) {
        return std::nullopt;
    }

    std::vector<double> var;
    const auto fit = dss_fit(obs, meta.n_params, var);
    if (!fit) return std::nullopt;

    AncSiteResult out;
    out.n = static_cast<int>(obs.size());
    out.mean_frac = sum_k / sum_n;
    out.la_prop = sum_la / (2.0 * static_cast<double>(out.n));
    out.phi = fit->second;
    out.phenotype = coef_test(fit->first, var, meta.n_params, kPhenotypeIdx);
    out.hapcount = coef_test(fit->first, var, meta.n_params, kHapcountIdx);
    return out;
}

std::optional<SiteResult> fit_site(const SiteRow& site, const CohortMeta& meta, int num_ancs) {
    SiteResult result;
    result.chr = site.chr;
    result.pos = site.pos;
    result.anc.resize(static_cast<std::size_t>(num_ancs));

    bool any = false;
    for (int anc = 0; anc < num_ancs; ++anc) {
        const auto fit = fit_ancestry_site(site, anc, meta);
        if (fit) {
            result.anc[static_cast<std::size_t>(anc)] = *fit;
            any = true;
        }
    }
    return any ? std::optional<SiteResult>(result) : std::nullopt;
}

struct TractReaders {
    std::vector<std::ifstream> counts;
    std::vector<std::ifstream> cov;
    std::vector<std::ifstream> hapcount;
};

void open_readers(TractReaders& readers, const Args& args) {
    readers.counts.resize(static_cast<std::size_t>(args.num_ancs));
    readers.cov.resize(static_cast<std::size_t>(args.num_ancs));
    readers.hapcount.resize(static_cast<std::size_t>(args.num_ancs));

    for (int anc = 0; anc < args.num_ancs; ++anc) {
        const std::string suffix = ".anc" + std::to_string(anc);
        readers.counts[static_cast<std::size_t>(anc)].open(args.tract_prefix + suffix + ".counts.txt");
        readers.cov[static_cast<std::size_t>(anc)].open(args.tract_prefix + suffix + ".cov.txt");
        readers.hapcount[static_cast<std::size_t>(anc)].open(args.tract_prefix + suffix + ".hapcount.txt");
        if (!readers.counts[static_cast<std::size_t>(anc)] ||
            !readers.cov[static_cast<std::size_t>(anc)] ||
            !readers.hapcount[static_cast<std::size_t>(anc)]) {
            throw std::runtime_error("Failed to open tract files for ancestry " + std::to_string(anc));
        }
    }
}

bool read_site_row(TractReaders& readers, const Args& args, const CohortMeta& meta, SiteRow& site) {
    std::string counts_line;
    if (!std::getline(readers.counts[0], counts_line)) return false;
    if (counts_line.empty()) return false;

    const auto counts_fields = split_line(counts_line, '\t');
    if (counts_fields.size() < 2) return false;

    site.chr = counts_fields[0];
    site.pos = std::stoi(counts_fields[1]);
    site.counts.assign(static_cast<std::size_t>(args.num_ancs), std::vector<double>(meta.samples.size(), kNaN));
    site.cov.assign(static_cast<std::size_t>(args.num_ancs), std::vector<double>(meta.samples.size(), kNaN));
    site.hapcount.assign(static_cast<std::size_t>(args.num_ancs), std::vector<double>(meta.samples.size(), kNaN));

    auto fill_from_fields = [&](const std::vector<std::string>& fields, std::vector<double>& dest) {
        for (std::size_t i = 0; i < meta.samples.size(); ++i) {
            const int idx = meta.samples[i].col_idx;
            dest[i] = idx < static_cast<int>(fields.size()) ? parse_field(fields[static_cast<std::size_t>(idx)]) : kNaN;
        }
    };

    fill_from_fields(counts_fields, site.counts[0]);

    for (int anc = 0; anc < args.num_ancs; ++anc) {
        if (anc == 0) continue;
        std::string line;
        if (!std::getline(readers.counts[static_cast<std::size_t>(anc)], line)) {
            throw std::runtime_error("Unexpected EOF in counts file for ancestry " + std::to_string(anc));
        }
        fill_from_fields(split_line(line, '\t'), site.counts[static_cast<std::size_t>(anc)]);
    }

    for (int anc = 0; anc < args.num_ancs; ++anc) {
        std::string line;
        if (!std::getline(readers.cov[static_cast<std::size_t>(anc)], line)) {
            throw std::runtime_error("Unexpected EOF in cov file for ancestry " + std::to_string(anc));
        }
        fill_from_fields(split_line(line, '\t'), site.cov[static_cast<std::size_t>(anc)]);
    }

    for (int anc = 0; anc < args.num_ancs; ++anc) {
        std::string line;
        if (!std::getline(readers.hapcount[static_cast<std::size_t>(anc)], line)) {
            throw std::runtime_error("Unexpected EOF in hapcount file for ancestry " + std::to_string(anc));
        }
        fill_from_fields(split_line(line, '\t'), site.hapcount[static_cast<std::size_t>(anc)]);
    }

    return true;
}

void write_header(std::ostream& out, int num_ancs) {
    out << "chr\tpos";
    for (int anc = 0; anc < num_ancs; ++anc) {
        const std::string a = "anc" + std::to_string(anc);
        out << "\tN_" << a
            << "\tmf_" << a
            << "\tLAprop_" << a
            << "\tphi_" << a
            << "\tbeta_pheno_" << a
            << "\tse_pheno_" << a
            << "\tpval_pheno_" << a
            << "\ttval_pheno_" << a
            << "\tbeta_LA_" << a
            << "\tse_LA_" << a
            << "\tpval_LA_" << a
            << "\ttval_LA_" << a;
    }
    out << '\n';
}

std::string format_field(double v, int precision = 6) {
    if (is_nan(v)) return ".";
    std::ostringstream oss;
    oss.setf(std::ios::fmtflags(0), std::ios::floatfield);
    oss.precision(precision);
    oss << v;
    return oss.str();
}

void write_field(std::ostream& out, const std::string& value) {
    out << '\t' << value;
}

void write_field(std::ostream& out, double v, int precision = 6) {
    write_field(out, format_field(v, precision));
}

void write_coef(std::ostream& out, const CoefResult& c) {
    write_field(out, c.beta, 6);
    write_field(out, c.se, 4);
    write_field(out, c.pvalue, 6);
    write_field(out, c.tval, 4);
}

void write_missing_ancestry(std::ostream& out) {
    constexpr int kFieldsPerAnc = 12;
    for (int k = 0; k < kFieldsPerAnc; ++k) write_field(out, ".");
}

void write_result_row(std::ostream& out, const SiteResult& r, int num_ancs) {
    out << r.chr << '\t' << r.pos;
    for (int anc = 0; anc < num_ancs; ++anc) {
        const auto& a = r.anc[static_cast<std::size_t>(anc)];
        if (a.n == 0) {
            write_missing_ancestry(out);
            continue;
        }
        write_field(out, std::to_string(a.n));
        write_field(out, a.mean_frac, 6);
        write_field(out, a.la_prop, 6);
        write_field(out, a.phi, 6);
        write_coef(out, a.phenotype);
        write_coef(out, a.hapcount);
    }
    out << '\n';
}

void run_analysis(const Args& args) {
#ifdef _OPENMP
    if (args.num_threads > 0) omp_set_num_threads(args.num_threads);
#endif

    TractReaders readers;
    open_readers(readers, args);

    std::string header_line;
    if (!std::getline(readers.counts[0], header_line)) {
        throw std::runtime_error("Empty counts file for ancestry 0");
    }
    const auto tract_header = split_line(header_line, '\t');
    const CohortMeta meta = load_metadata(args.metadata, tract_header, args.case_label, args.control_label);

    for (int anc = 1; anc < args.num_ancs; ++anc) {
        std::string line;
        if (!std::getline(readers.counts[static_cast<std::size_t>(anc)], line)) {
            throw std::runtime_error("Empty counts file for ancestry " + std::to_string(anc));
        }
    }
    for (int anc = 0; anc < args.num_ancs; ++anc) {
        std::string line;
        if (!std::getline(readers.cov[static_cast<std::size_t>(anc)], line)) {
            throw std::runtime_error("Empty cov file for ancestry " + std::to_string(anc));
        }
        if (!std::getline(readers.hapcount[static_cast<std::size_t>(anc)], line)) {
            throw std::runtime_error("Empty hapcount file for ancestry " + std::to_string(anc));
        }
    }

    std::ofstream out(args.output);
    if (!out) throw std::runtime_error("Cannot open output: " + args.output);
    write_header(out, args.num_ancs);

    std::vector<SiteRow> batch;
    batch.reserve(args.batch_size);
    std::size_t rows_read = 0;
    std::size_t rows_written = 0;

    while (true) {
        batch.clear();
        SiteRow site;
        while (batch.size() < args.batch_size && read_site_row(readers, args, meta, site)) {
            batch.push_back(site);
        }
        if (batch.empty()) break;
        rows_read += batch.size();

        std::vector<std::optional<SiteResult>> batch_out(batch.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i < batch.size(); ++i) {
            batch_out[i] = fit_site(batch[i], meta, args.num_ancs);
        }

        for (const auto& fit : batch_out) {
            if (fit) {
                write_result_row(out, *fit, args.num_ancs);
                ++rows_written;
            }
        }

        if (rows_read % (args.batch_size * 5) == 0) {
            std::cerr << "Processed " << rows_read << " loci, wrote " << rows_written << " results...\n";
        }
    }

    std::cerr << "la-dml v" << kVersion << '\n'
              << "Read " << rows_read << " loci, wrote " << rows_written << " tested CpGs\n"
              << "Output: " << args.output << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Args args = parse_args(argc, argv);
        run_analysis(args);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << '\n';
        return 1;
    }
}
