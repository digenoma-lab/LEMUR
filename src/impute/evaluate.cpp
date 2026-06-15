#include "impute_methylation/evaluate.hpp"

#include "impute_methylation/haplotype_target.hpp"
#include "impute_methylation/tsv_io.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace impute_methylation {

namespace {

uint64_t splitmix64(uint64_t& state) {
    uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

bool should_mask(const std::string& chr, int pos, const std::string& y_col, uint64_t seed,
                 double mask_fraction) {
    uint64_t state = seed;
    state ^= std::hash<std::string>{}(chr) + 0x9e3779b97f4a7c15ULL + (state << 6) + (state >> 2);
    state ^= static_cast<uint64_t>(pos) * 0x632be59bd9b4e019ULL;
    state ^= std::hash<std::string>{}(y_col) + 0x9e3779b97f4a7c15ULL;
    const uint64_t r = splitmix64(state);
    const uint64_t threshold = static_cast<uint64_t>(mask_fraction * 10000.0);
    return (r % 10000) < threshold;
}

struct MaskedHoldout {
    std::string chr;
    int pos = 0;
    std::string y_col;
    std::string out_col;
    double true_frac = kNaN;
};

std::string default_sidecar_path(const std::string& input, const std::string& suffix) {
    return input + suffix;
}

void write_masked_row(std::ostream& out, const std::vector<std::string>& input_header,
                      const std::vector<std::string>& fields, const HaplotypeTarget& target,
                      bool mask_this_row) {
    for (std::size_t i = 0; i < input_header.size(); ++i) {
        if (i) out << '\t';
        const std::string& col = input_header[i];
        const bool masked =
            mask_this_row &&
            (col == target.y_col || col == target.n_col || col == target.pct_col);
        if (masked) {
            out << '.';
        } else {
            out << (i < fields.size() ? fields[i] : ".");
        }
    }
    out << '\n';
}

struct MetricAccum {
    double sum_sq_err = 0.0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_xx = 0.0;
    double sum_yy = 0.0;
    double sum_xy = 0.0;
    std::size_t n = 0;
};

void add_pair(MetricAccum& acc, double pred, double truth) {
    const double err = pred - truth;
    acc.sum_sq_err += err * err;
    acc.sum_x += pred;
    acc.sum_y += truth;
    acc.sum_xx += pred * pred;
    acc.sum_yy += truth * truth;
    acc.sum_xy += pred * truth;
    ++acc.n;
}

EvaluateMetrics finalize_metrics(MetricAccum& acc, std::size_t n_masked) {
    EvaluateMetrics m;
    m.n_masked = n_masked;
    m.n_scored = acc.n;
    if (acc.n == 0) {
        m.mse = kNaN;
        m.pearson = kNaN;
        return m;
    }
    m.mse = acc.sum_sq_err / static_cast<double>(acc.n);
    if (acc.n < 2) {
        m.pearson = kNaN;
        return m;
    }
    const double n = static_cast<double>(acc.n);
    const double cov = acc.sum_xy - acc.sum_x * acc.sum_y / n;
    const double var_x = acc.sum_xx - acc.sum_x * acc.sum_x / n;
    const double var_y = acc.sum_yy - acc.sum_y * acc.sum_y / n;
    if (var_x <= 0.0 || var_y <= 0.0) {
        m.pearson = kNaN;
    } else {
        m.pearson = cov / std::sqrt(var_x * var_y);
    }
    return m;
}

std::string site_key(const std::string& chr, int pos) {
    return chr + '\t' + std::to_string(pos);
}

std::string format_metric(double value) {
    if (is_nan(value)) return "nan";
    std::ostringstream os;
    os << std::fixed << std::setprecision(6) << value;
    return os.str();
}

std::string cohort_work_dir(const std::string& output_path) {
    const std::filesystem::path out(output_path);
    return (out.parent_path() / (out.stem().string() + ".evaluate.work")).string();
}

std::string target_label(const HaplotypeTarget& target) {
    return target.hap.empty() ? target.sample_id : target.sample_id + "." + target.hap;
}

std::string target_sidecar_path(const std::string& work_dir, const HaplotypeTarget& target,
                                const std::string& suffix) {
    return work_dir + "/" + target_label(target) + suffix;
}

EvaluateCohortRow metrics_to_row(const HaplotypeTarget& target, const EvaluateOptions& opts,
                                 const EvaluateMetrics& metrics) {
    EvaluateCohortRow row;
    row.sample = target_label(target);
    row.window_size = opts.impute.window_bp;
    row.n_neighbors = opts.impute.min_neighbors;
    row.a = opts.impute.alpha;
    row.b = opts.impute.beta;
    row.chr = opts.chromosome;
    row.pearson = metrics.pearson;
    row.mse = metrics.mse;
    row.count_masked = metrics.n_masked;
    row.count_imputed = metrics.n_scored;
    return row;
}

void remove_if_exists(const std::string& path) {
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

}  // namespace

EvaluateMetrics run_evaluate(const std::string& input_path, const EvaluateOptions& opts) {
    if (opts.y_col.empty()) {
        throw std::runtime_error(opts.impute.sample_mode
                                     ? "EvaluateOptions.y_col is required (e.g. CHI08A.counts)"
                                     : "EvaluateOptions.y_col is required (e.g. CHI08A.hap1_counts)");
    }

    const std::string masked_path =
        opts.masked_path.empty() ? default_sidecar_path(input_path, ".masked.tsv") : opts.masked_path;
    const std::string imputed_path = opts.imputed_path.empty()
                                         ? default_sidecar_path(input_path, ".imputed.eval.tsv")
                                         : opts.imputed_path;

    std::ifstream in(input_path);
    if (!in) throw std::runtime_error("Cannot open input: " + input_path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty input: " + input_path);

    const std::vector<std::string> input_header = split_tab(line);
    const HaplotypeTarget target =
        opts.impute.sample_mode ? find_sample_target(input_header, opts.y_col)
                                : find_haplotype_target(input_header, opts.y_col);
    const std::vector<HaplotypeTarget> targets = {target};

    std::unordered_map<std::string, MaskedHoldout> holdouts;
    holdouts.reserve(10000);
    std::size_t n_rows_chr = 0;

    std::ofstream masked_out(masked_path);
    if (!masked_out) throw std::runtime_error("Cannot open masked output: " + masked_path);
    write_header_line(masked_out, input_header);

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto fields = split_tab(line);
        if (fields.size() < 2) continue;
        if (fields.size() < input_header.size()) fields.resize(input_header.size(), ".");

        const std::string& chr = fields[0];
        if (chr != opts.chromosome) continue;

        ++n_rows_chr;
        const int pos = std::stoi(fields[1]);

        const double y =
            parse_double_field(fields[static_cast<std::size_t>(target.y_idx)]);
        const double n =
            parse_double_field(fields[static_cast<std::size_t>(target.n_idx)]);
        double pct = kNaN;
        if (target.pct_idx >= 0) {
            pct = parse_double_field(fields[static_cast<std::size_t>(target.pct_idx)]);
        }
        const double truth = ground_truth_fraction(y, n, pct);

        bool mask_row = false;
        if (!is_nan(truth) &&
            should_mask(chr, pos, target.y_col, opts.seed, opts.mask_fraction)) {
            mask_row = true;
            MaskedHoldout h;
            h.chr = chr;
            h.pos = pos;
            h.y_col = target.y_col;
            h.out_col = target.out_col;
            h.true_frac = truth;
            holdouts[site_key(chr, pos)] = h;
        }

        write_masked_row(masked_out, input_header, fields, target, mask_row);
    }

    masked_out.close();
    in.close();

    if (n_rows_chr == 0) {
        throw std::runtime_error("No rows on chromosome: " + opts.chromosome);
    }

    const std::size_t n_masked = holdouts.size();

    if (!opts.quiet) {
        std::cerr << "Chromosome " << opts.chromosome << " (" << n_rows_chr << " sites), column "
                  << target.y_col << ": masked " << n_masked << " (~"
                  << opts.mask_fraction * 100 << "% of valid on this chr) → " << masked_path
                  << '\n';
    }

    ImputeOptions impute_opts = opts.impute;
    impute_opts.num_threads = 1;
    stream_beta_binomial_impute_targets(masked_path, imputed_path, targets, impute_opts);

    std::ifstream imp_in(imputed_path);
    if (!imp_in) throw std::runtime_error("Cannot open imputed: " + imputed_path);
    if (!std::getline(imp_in, line)) throw std::runtime_error("Empty imputed file");

    const std::vector<std::string> out_header = split_tab(line);
    const int out_idx = column_index(out_header, target.out_col);
    if (out_idx < 0) {
        throw std::runtime_error("Imputed column missing: " + target.out_col);
    }

    MetricAccum acc;

    while (std::getline(imp_in, line)) {
        if (line.empty()) continue;
        auto fields = split_tab(line);
        if (fields.size() < 2) continue;

        const std::string& chr = fields[0];
        const int pos = std::stoi(fields[1]);

        const auto it = holdouts.find(site_key(chr, pos));
        if (it == holdouts.end()) continue;

        if (static_cast<std::size_t>(out_idx) >= fields.size()) continue;

        const double pred = parse_fraction_field(fields[static_cast<std::size_t>(out_idx)]);
        if (is_nan(pred)) continue;

        add_pair(acc, pred, it->second.true_frac);
        holdouts.erase(it);
    }

    const EvaluateMetrics metrics = finalize_metrics(acc, n_masked);

    if (!opts.quiet) {
        std::cerr << "\n=== Evaluation (masked hold-out, MSE & Pearson) ===\n";
        std::cerr << std::fixed << std::setprecision(6);
        std::cerr << "chromosome=" << opts.chromosome << " column=" << target.y_col
                  << " mask_fraction=" << opts.mask_fraction << " seed=" << opts.seed
                  << " window=" << opts.impute.window_bp
                  << " min_neighbors=" << opts.impute.min_neighbors << "\n\n";
        std::cerr << target.out_col << "\tmasked=" << metrics.n_masked
                  << "\tscored=" << metrics.n_scored << "\tMSE=" << metrics.mse
                  << "\tPearson=" << metrics.pearson << '\n';
        std::cerr << "Imputed file: " << imputed_path << '\n';
    }

    return metrics;
}

std::vector<EvaluateCohortRow> run_evaluate_cohort(const std::string& input_path,
                                                   const EvaluateCohortOptions& opts) {
    if (opts.output_path.empty()) {
        throw std::runtime_error("EvaluateCohortOptions.output_path is required");
    }

    std::ifstream in(input_path);
    if (!in) throw std::runtime_error("Cannot open input: " + input_path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty input: " + input_path);
    in.close();

    const std::vector<std::string> input_header = split_tab(line);
    const std::vector<HaplotypeTarget> targets =
        opts.impute.sample_mode ? discover_sample_targets(input_header)
                               : discover_haplotype_targets(input_header);
    if (targets.empty()) {
        throw std::runtime_error(opts.impute.sample_mode
                                     ? "No {sample}.counts columns found in header"
                                     : "No {sample}.hap{1,2}_counts columns found in header");
    }

    const std::string work_dir = cohort_work_dir(opts.output_path);
    std::filesystem::create_directories(work_dir);

    std::vector<EvaluateCohortRow> rows(targets.size());

#ifdef _OPENMP
    const int prev_threads = omp_get_max_threads();
    if (opts.impute.num_threads > 0) {
        omp_set_num_threads(opts.impute.num_threads);
    }
#endif

#pragma omp parallel for schedule(dynamic) default(none) shared( \
    input_path, opts, targets, rows, work_dir) if (targets.size() > 1)
    for (std::size_t t = 0; t < targets.size(); ++t) {
        EvaluateOptions single;
        single.y_col = targets[t].y_col;
        single.chromosome = opts.chromosome;
        single.mask_fraction = opts.mask_fraction;
        single.seed = opts.seed;
        single.impute = opts.impute;
        single.impute.num_threads = 1;
        single.quiet = true;
        single.masked_path = target_sidecar_path(work_dir, targets[t], ".masked.tsv");
        single.imputed_path = target_sidecar_path(work_dir, targets[t], ".imputed.eval.tsv");

        const EvaluateMetrics metrics = run_evaluate(input_path, single);
        rows[t] = metrics_to_row(targets[t], single, metrics);

        if (!opts.keep_sidecars) {
            remove_if_exists(single.masked_path);
            remove_if_exists(single.imputed_path);
        }
    }

#ifdef _OPENMP
    if (opts.impute.num_threads > 0) {
        omp_set_num_threads(prev_threads);
    }
#endif

    if (!opts.keep_sidecars) {
        std::error_code ec;
        std::filesystem::remove_all(work_dir, ec);
    }

    write_evaluate_cohort_tsv(opts.output_path, rows);

    std::cerr << "Evaluated " << rows.size() << " "
              << (opts.impute.sample_mode ? "sample" : "haplotype") << " columns on "
              << opts.chromosome << " → " << opts.output_path << '\n';

    return rows;
}

void write_evaluate_cohort_tsv(const std::string& path,
                               const std::vector<EvaluateCohortRow>& rows) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open cohort output: " + path);

    out << "sample\twindow_size\tn_neighbors\ta\tb\tchr\tpearson\tmse\tcount_masked\tcount_imputed\n";
    for (const auto& r : rows) {
        out << r.sample << '\t' << r.window_size << '\t' << r.n_neighbors << '\t' << r.a << '\t'
            << r.b << '\t' << r.chr << '\t' << format_metric(r.pearson) << '\t'
            << format_metric(r.mse) << '\t' << r.count_masked << '\t' << r.count_imputed << '\n';
    }
}

}  // namespace impute_methylation
