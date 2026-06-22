// tracts — extract ancestry-specific methylation tracts from phased haplotype matrices.
//
// Methylation analogue of Tractor's extract_tracts.py: combines haplotype-level
// methylation (counts/coverage or imputed fractions) with MSP local-ancestry
// calls to produce per-ancestry matrices for downstream la-dml analysis.
//
// For each ancestral population i, writes:
//   <prefix>.anc<i>.counts.txt   summed modified counts on ancestry-i haplotypes
//   <prefix>.anc<i>.cov.txt      summed coverage on ancestry-i haplotypes
//   <prefix>.anc<i>.hapcount.txt number of haplotypes from ancestry i (0-2)
//
// Usage example:
//   tracts --methylation chr1.tsv --msp lai_chr1.msp.tsv
//       --num-ancs 3 --output-dir out/ [--output-prefix chr1]

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr const char* kVersion = "1.1.0";

struct Args {
    std::string methylation;
    std::string msp;
    int num_ancs = 0;
    std::string output_dir;
    std::string output_prefix;
};

struct MspWindow {
    std::string chrom;
    int start = 0;
    int end = 0;
    std::vector<int> calls;
    bool valid = false;
};

struct HapColumns {
    int counts = -1;
    int cov = -1;
    int frac = -1;
};

struct SampleColumns {
    std::string sample_id;
    HapColumns hap1;
    HapColumns hap2;
    bool use_fractions = false;
};

void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog
        << " --methylation <haplotype methylation TSV>\n"
        << "    --msp <local ancestry MSP file>\n"
        << "    --num-ancs <number of ancestral populations>\n"
        << "    [--output-dir <directory>]\n"
        << "    [--output-prefix <prefix>]\n\n"
        << "Options:\n"
        << "  --methylation       Filtered haplotype methylation matrix "
           "(counts/cov or hap*_frac columns)\n"
        << "  --msp               RFMix/G-Nomix MSP file (*.msp or *.msp.tsv)\n"
        << "  --num-ancs          Number of ancestral populations\n"
        << "  --output-dir        Output directory (default: methylation file directory)\n"
        << "  --output-prefix     Output filename prefix (default: methylation basename)\n";
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream iss(line);
    while (std::getline(iss, field, '\t')) fields.push_back(field);
    return fields;
}

bool is_missing(const std::string& s) {
    return s.empty() || s == ".";
}

bool file_exists(const std::string& path) {
    std::ifstream in(path);
    return static_cast<bool>(in);
}

std::string basename_no_ext(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    const std::string name = (slash == std::string::npos) ? path : path.substr(slash + 1);
    const std::size_t dot = name.find('.');
    return (dot == std::string::npos) ? name : name.substr(0, dot);
}

std::string parent_dir(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    if (slash == std::string::npos) return ".";
    return path.substr(0, slash);
}

std::string join_path(const std::string& dir, const std::string& name) {
    if (dir.empty() || dir == ".") return name;
    if (dir.back() == '/') return dir + name;
    return dir + "/" + name;
}

std::string msp_sample_id(const std::string& col) {
  // CHI01.0 -> CHI01
    const std::size_t dot = col.rfind('.');
    if (dot == std::string::npos || dot == 0) return col;
    const std::string suffix = col.substr(dot + 1);
    if (suffix != "0" && suffix != "1") return col;
    return col.substr(0, dot);
}

bool parse_hap_counts_column(const std::string& col, std::string& sample_id, int& hap_num) {
    const std::size_t h1 = col.rfind(".hap1_counts");
    const std::size_t h2 = col.rfind(".hap2_counts");
    if (h1 != std::string::npos && h1 + 12 == col.size()) {
        sample_id = col.substr(0, h1);
        hap_num = 1;
        return true;
    }
    if (h2 != std::string::npos && h2 + 12 == col.size()) {
        sample_id = col.substr(0, h2);
        hap_num = 2;
        return true;
    }
    return false;
}

bool parse_hap_cov_column(const std::string& col, std::string& sample_id, int& hap_num) {
    const std::size_t h1 = col.rfind(".hap1_cov");
    const std::size_t h2 = col.rfind(".hap2_cov");
    if (h1 != std::string::npos && h1 + 9 == col.size()) {
        sample_id = col.substr(0, h1);
        hap_num = 1;
        return true;
    }
    if (h2 != std::string::npos && h2 + 9 == col.size()) {
        sample_id = col.substr(0, h2);
        hap_num = 2;
        return true;
    }
    return false;
}

bool parse_hap_frac_column(const std::string& col, std::string& sample_id, int& hap_num) {
    const std::size_t h1 = col.rfind(".hap1_frac");
    const std::size_t h2 = col.rfind(".hap2_frac");
    if (h1 != std::string::npos && h1 + 10 == col.size()) {
        sample_id = col.substr(0, h1);
        hap_num = 1;
        return true;
    }
    if (h2 != std::string::npos && h2 + 10 == col.size()) {
        sample_id = col.substr(0, h2);
        hap_num = 2;
        return true;
    }
    return false;
}

std::string normalize_sample_id(const std::string& raw) {
    const std::size_t chr_pos = raw.find(".chr");
    if (chr_pos != std::string::npos) return raw.substr(0, chr_pos);
    return raw;
}

double parse_double(const std::string& s) {
    if (is_missing(s)) return std::numeric_limits<double>::quiet_NaN();
    try {
        return std::stod(s);
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

const std::string& field_as_string(const std::vector<std::string>& fields, int index) {
    static const std::string kMissing = ".";
    if (index < 0 || index >= static_cast<int>(fields.size())) return kMissing;
    return fields[index];
}

bool is_nan(double x) { return std::isnan(x); }

Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need_value = [&](const char* flag) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error(std::string("Missing value for ") + flag);
            }
            return argv[++i];
        };

        if (arg == "--methylation") {
            args.methylation = need_value("--methylation");
        } else if (arg == "--msp") {
            args.msp = need_value("--msp");
        } else if (arg == "--num-ancs") {
            args.num_ancs = std::stoi(need_value("--num-ancs"));
        } else if (arg == "--output-dir") {
            args.output_dir = need_value("--output-dir");
        } else if (arg == "--output-prefix") {
            args.output_prefix = need_value("--output-prefix");
        } else if (arg == "-h" || arg == "--help") {
            usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (args.methylation.empty() || args.msp.empty() || args.num_ancs <= 0) {
        usage(argv[0]);
        throw std::runtime_error("Required arguments: --methylation, --msp, --num-ancs");
    }
    if (!file_exists(args.methylation)) {
        throw std::runtime_error("Methylation file not found: " + args.methylation);
    }
    if (!file_exists(args.msp)) {
        throw std::runtime_error("MSP file not found: " + args.msp);
    }
    if (args.output_dir.empty()) args.output_dir = parent_dir(args.methylation);
    if (args.output_prefix.empty()) args.output_prefix = basename_no_ext(args.methylation);
    return args;
}

std::vector<std::string> parse_msp_samples(std::istream& msp_in) {
    std::string line;
    while (std::getline(msp_in, line)) {
        if (line.empty()) continue;
        if (line[0] != '#') continue;
        if (line.rfind("#chm", 0) != 0) continue;

        const auto fields = split_tab(line.substr(1));
        if (fields.size() < 7) {
            throw std::runtime_error("Malformed MSP header line");
        }

        std::vector<std::string> samples;
        samples.reserve((fields.size() - 6) / 2);
        for (std::size_t i = 6; i < fields.size(); i += 2) {
            const std::string sample = msp_sample_id(fields[i]);
            if (samples.empty() || samples.back() != sample) {
                samples.push_back(sample);
            }
        }
        if (samples.empty()) {
            throw std::runtime_error("No samples found in MSP header");
        }
        return samples;
    }
    throw std::runtime_error("MSP file is missing #chm header line");
}

std::unordered_map<std::string, SampleColumns> parse_methylation_header(
    const std::vector<std::string>& header) {
    if (header.size() < 3 || header[0] != "chr" || header[1] != "pos") {
        throw std::runtime_error("Methylation header must start with chr\\tpos");
    }

    std::unordered_map<std::string, SampleColumns> by_sample;
    for (std::size_t i = 2; i < header.size(); ++i) {
        const std::string& col = header[i];
        std::string raw_sample;
        int hap_num = 0;

        if (parse_hap_counts_column(col, raw_sample, hap_num)) {
            const std::string sample_id = normalize_sample_id(raw_sample);
            auto& sample = by_sample[sample_id];
            sample.sample_id = sample_id;
            if (hap_num == 1) sample.hap1.counts = static_cast<int>(i);
            else sample.hap2.counts = static_cast<int>(i);
        } else if (parse_hap_cov_column(col, raw_sample, hap_num)) {
            const std::string sample_id = normalize_sample_id(raw_sample);
            auto& sample = by_sample[sample_id];
            sample.sample_id = sample_id;
            if (hap_num == 1) sample.hap1.cov = static_cast<int>(i);
            else sample.hap2.cov = static_cast<int>(i);
        } else if (parse_hap_frac_column(col, raw_sample, hap_num)) {
            const std::string sample_id = normalize_sample_id(raw_sample);
            auto& sample = by_sample[sample_id];
            sample.sample_id = sample_id;
            sample.use_fractions = true;
            if (hap_num == 1) sample.hap1.frac = static_cast<int>(i);
            else sample.hap2.frac = static_cast<int>(i);
        }
    }

    if (by_sample.empty()) {
        throw std::runtime_error(
            "No haplotype methylation columns found "
            "(expected *.hap1_counts / *.hap1_frac style columns)");
    }
    return by_sample;
}

std::vector<SampleColumns> order_samples(
    const std::vector<std::string>& msp_samples,
    const std::unordered_map<std::string, SampleColumns>& by_sample) {
    std::vector<SampleColumns> ordered;
    ordered.reserve(msp_samples.size());
    for (const auto& sample_id : msp_samples) {
        const auto it = by_sample.find(sample_id);
        if (it == by_sample.end()) {
            throw std::runtime_error("Sample in MSP not found in methylation header: " + sample_id);
        }
        ordered.push_back(it->second);
    }
    return ordered;
}

bool read_msp_window(std::istream& msp_in, MspWindow& window) {
    window = MspWindow{};
    std::string line;
    while (std::getline(msp_in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tab(line);
        if (fields.size() < 7) continue;

        window.chrom = fields[0];
        window.start = std::stoi(fields[1]);
        window.end = std::stoi(fields[2]);
        window.calls.clear();
        window.calls.reserve(fields.size() - 6);
        for (std::size_t i = 6; i < fields.size(); ++i) {
            window.calls.push_back(std::stoi(fields[i]));
        }
        window.valid = true;
        return true;
    }
    return false;
}

bool position_in_window(const std::string& chrom, int pos, const MspWindow& window) {
    return window.valid && chrom == window.chrom && window.start <= pos && pos < window.end;
}

struct OutputStreams {
    std::vector<std::ofstream> counts;
    std::vector<std::ofstream> cov;
    std::vector<std::ofstream> hapcount;
};

OutputStreams open_outputs(const Args& args) {
    OutputStreams out;
    out.counts.resize(args.num_ancs);
    out.cov.resize(args.num_ancs);
    out.hapcount.resize(args.num_ancs);

    for (int anc = 0; anc < args.num_ancs; ++anc) {
        const std::string base = join_path(args.output_dir, args.output_prefix + ".anc" + std::to_string(anc));
        out.counts[anc].open(base + ".counts.txt");
        out.cov[anc].open(base + ".cov.txt");
        out.hapcount[anc].open(base + ".hapcount.txt");
        if (!out.counts[anc] || !out.cov[anc] || !out.hapcount[anc]) {
            throw std::runtime_error("Failed to open output files for ancestry " + std::to_string(anc));
        }
    }
    return out;
}

void write_headers(OutputStreams& out,
                   const Args& args,
                   const std::vector<std::string>& sample_ids) {
    const std::string locus_header = "chr\tpos";
    std::string sample_header = locus_header;
    for (const auto& sample_id : sample_ids) {
        sample_header += '\t';
        sample_header += sample_id;
    }
    sample_header += '\n';

    for (int anc = 0; anc < args.num_ancs; ++anc) {
        out.counts[anc] << sample_header;
        out.cov[anc] << sample_header;
        out.hapcount[anc] << sample_header;
    }
}

std::string format_count(double value) {
    if (is_nan(value)) return ".";
    const auto rounded = static_cast<long long>(std::llround(value));
    if (rounded < 0) return ".";
    return std::to_string(rounded);
}

void append_hap_fraction_as_counts(const std::string& frac_field,
                                   const std::string& cov_field,
                                   bool ancestry_match,
                                   double& count_sum,
                                   double& cov_sum) {
    if (!ancestry_match) return;
    const double frac = parse_double(frac_field);
    const double cov = parse_double(cov_field);
    if (is_nan(frac) || is_nan(cov) || cov <= 0.0) return;
    count_sum += frac * cov;
    cov_sum += cov;
}

void append_hap_counts(const std::string& counts_field,
                       const std::string& cov_field,
                       bool ancestry_match,
                       double& count_sum,
                       double& cov_sum) {
    if (!ancestry_match) return;
    const double counts = parse_double(counts_field);
    const double cov = parse_double(cov_field);
    if (is_nan(counts) || is_nan(cov) || cov <= 0.0) return;
    count_sum += counts;
    cov_sum += cov;
}

void extract_tracts_methylation(const Args& args) {
    std::ifstream msp_in(args.msp);
    if (!msp_in) throw std::runtime_error("Failed to open MSP file: " + args.msp);

    const std::vector<std::string> msp_samples = parse_msp_samples(msp_in);

    std::ifstream methyl_in(args.methylation);
    if (!methyl_in) throw std::runtime_error("Failed to open methylation file: " + args.methylation);

    std::string header_line;
    if (!std::getline(methyl_in, header_line)) {
        throw std::runtime_error("Methylation file is empty");
    }
    const auto header = split_tab(header_line);
    const auto by_sample = parse_methylation_header(header);
    const std::vector<SampleColumns> samples = order_samples(msp_samples, by_sample);

    OutputStreams out = open_outputs(args);
    write_headers(out, args, msp_samples);

    MspWindow window;
    bool msp_eof = !read_msp_window(msp_in, window);

    std::size_t rows_written = 0;
    std::size_t rows_skipped = 0;
    std::string line;

    while (std::getline(methyl_in, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() < 2) continue;

        const std::string chrom = fields[0];
        const int pos = std::stoi(fields[1]);
        bool skip_line = false;

        while (!position_in_window(chrom, pos, window)) {
            if (window.valid && chrom == window.chrom && window.start > pos) {
                skip_line = true;
                break;
            }
            if (msp_eof) break;
            msp_eof = !read_msp_window(msp_in, window);
            if (!msp_eof && chrom == window.chrom && window.start > pos) {
                skip_line = true;
                break;
            }
        }

        if (skip_line || !position_in_window(chrom, pos, window)) {
            rows_skipped++;
            continue;
        }

        if (static_cast<int>(window.calls.size()) != static_cast<int>(samples.size()) * 2) {
            throw std::runtime_error(
                "MSP call count does not match number of samples at position " + std::to_string(pos));
        }

        std::vector<std::string> count_values(args.num_ancs);
        std::vector<std::string> cov_values(args.num_ancs);
        std::vector<std::string> hap_values(args.num_ancs);

        for (std::size_t sample_idx = 0; sample_idx < samples.size(); ++sample_idx) {
            const SampleColumns& sample = samples[sample_idx];
            const int call_hap1 = window.calls[2 * sample_idx];
            const int call_hap2 = window.calls[2 * sample_idx + 1];

            std::vector<int> hap_count(args.num_ancs, 0);
            std::vector<double> count_sum(args.num_ancs, 0.0);
            std::vector<double> cov_sum(args.num_ancs, 0.0);

            if (sample.use_fractions) {
                for (int anc = 0; anc < args.num_ancs; ++anc) {
                    if (sample.hap1.frac >= 0 && sample.hap1.cov >= 0) {
                        append_hap_fraction_as_counts(field_as_string(fields, sample.hap1.frac),
                                                      field_as_string(fields, sample.hap1.cov),
                                                      call_hap1 == anc,
                                                      count_sum[anc],
                                                      cov_sum[anc]);
                    }
                    if (sample.hap2.frac >= 0 && sample.hap2.cov >= 0) {
                        append_hap_fraction_as_counts(field_as_string(fields, sample.hap2.frac),
                                                      field_as_string(fields, sample.hap2.cov),
                                                      call_hap2 == anc,
                                                      count_sum[anc],
                                                      cov_sum[anc]);
                    }

                    if (call_hap1 == anc) hap_count[anc] += 1;
                    if (call_hap2 == anc) hap_count[anc] += 1;
                }
            } else {
                for (int anc = 0; anc < args.num_ancs; ++anc) {
                    if (sample.hap1.counts >= 0 && sample.hap1.cov >= 0) {
                        append_hap_counts(field_as_string(fields, sample.hap1.counts),
                                          field_as_string(fields, sample.hap1.cov),
                                          call_hap1 == anc,
                                          count_sum[anc],
                                          cov_sum[anc]);
                    }
                    if (sample.hap2.counts >= 0 && sample.hap2.cov >= 0) {
                        append_hap_counts(field_as_string(fields, sample.hap2.counts),
                                          field_as_string(fields, sample.hap2.cov),
                                          call_hap2 == anc,
                                          count_sum[anc],
                                          cov_sum[anc]);
                    }

                    if (call_hap1 == anc) hap_count[anc] += 1;
                    if (call_hap2 == anc) hap_count[anc] += 1;
                }
            }

            for (int anc = 0; anc < args.num_ancs; ++anc) {
                if (sample_idx > 0) {
                    count_values[anc] += '\t';
                    cov_values[anc] += '\t';
                    hap_values[anc] += '\t';
                }

                if (cov_sum[anc] > 0.0) {
                    count_values[anc] += format_count(count_sum[anc]);
                    cov_values[anc] += format_count(cov_sum[anc]);
                } else {
                    count_values[anc] += ".";
                    cov_values[anc] += ".";
                }
                hap_values[anc] += std::to_string(hap_count[anc]);
            }
        }

        const std::string locus_prefix = chrom + '\t' + std::to_string(pos);
        for (int anc = 0; anc < args.num_ancs; ++anc) {
            out.counts[anc] << locus_prefix << '\t' << count_values[anc] << '\n';
            out.cov[anc] << locus_prefix << '\t' << cov_values[anc] << '\n';
            out.hapcount[anc] << locus_prefix << '\t' << hap_values[anc] << '\n';
        }
        rows_written++;
    }

    std::cerr << "tracts v" << kVersion << '\n'
              << "Wrote " << rows_written << " loci across " << args.num_ancs << " ancestries\n"
              << "Skipped " << rows_skipped << " loci outside MSP windows\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Args args = parse_args(argc, argv);
        extract_tracts_methylation(args);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << '\n';
        return 1;
    }
}
