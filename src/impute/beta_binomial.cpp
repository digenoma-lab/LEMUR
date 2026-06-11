#include "impute_methylation/beta_binomial.hpp"

#include "impute_methylation/haplotype_target.hpp"
#include "impute_methylation/tsv_io.hpp"

#include <deque>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace impute_methylation {

namespace {

struct SampleImputeGroup {
    std::vector<std::size_t> target_indices;
    std::vector<std::deque<WindowSite>> windows;
};

std::vector<SampleImputeGroup> group_targets_by_sample(
    const std::vector<HaplotypeTarget>& targets) {
    std::vector<SampleImputeGroup> groups;
    std::unordered_map<std::string, std::size_t> sample_to_group;
    sample_to_group.reserve(targets.size() / 2);
    groups.reserve(targets.size() / 2);

    for (std::size_t t = 0; t < targets.size(); ++t) {
        const auto [it, inserted] = sample_to_group.emplace(targets[t].sample_id, groups.size());
        if (inserted) {
            groups.emplace_back();
        }
        SampleImputeGroup& group = groups[it->second];
        group.target_indices.push_back(t);
        group.windows.emplace_back();
    }
    return groups;
}

void process_sample_group(SampleImputeGroup& group, const std::vector<HaplotypeTarget>& targets,
                          const ImputeOptions& opts, const std::vector<std::string>& fields,
                          int pos, bool new_chromosome,
                          std::vector<ImputeResult>& imputed_values) {
    if (new_chromosome) {
        for (auto& window : group.windows) {
            window.clear();
        }
    }

    for (std::size_t i = 0; i < group.target_indices.size(); ++i) {
        const std::size_t t = group.target_indices[i];
        const auto& target = targets[t];
        auto& window = group.windows[i];

        while (!window.empty() && pos - window.front().pos > opts.window_bp) {
            window.pop_front();
        }

        WindowSite site;
        site.pos = pos;
        site.y = parse_double_field(fields[static_cast<std::size_t>(target.y_idx)]);
        site.n = parse_double_field(fields[static_cast<std::size_t>(target.n_idx)]);
        window.push_back(site);

        double pct = kNaN;
        if (target.pct_idx >= 0) {
            pct = parse_double_field(fields[static_cast<std::size_t>(target.pct_idx)]);
        }
        imputed_values[t] = impute_from_window(window, opts, site.y, site.n, pct);
    }
}

}  // namespace

std::string format_fraction(double value) {
    if (is_nan(value)) return ".";
    constexpr double kEps = 1e-9;
    if (value <= kEps) return "0";
    if (value >= 1.0 - kEps) return "1";
    std::ostringstream os;
    os << std::fixed << std::setprecision(4) << value;
    std::string s = os.str();
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    return s;
}

void stream_beta_binomial_impute_targets(const std::string& input_path,
                                         const std::string& output_path,
                                         const std::vector<HaplotypeTarget>& targets,
                                         const ImputeOptions& opts) {
    if (targets.empty()) {
        throw std::runtime_error("No haplotype targets to impute");
    }

    std::ifstream in(input_path);
    if (!in) throw std::runtime_error("Cannot open input: " + input_path);

    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("Cannot open output: " + output_path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty input: " + input_path);

    const std::vector<std::string> input_header = split_tab(line);

    const OutputColumnPlan plan = build_output_plan(input_header, targets, opts.mode);
    write_header_line(out, plan.header);

    std::vector<SampleImputeGroup> sample_groups = group_targets_by_sample(targets);
    std::vector<ImputeResult> imputed_values(targets.size());
    std::string current_chr;

#ifdef _OPENMP
    const int prev_threads = omp_get_max_threads();
    if (opts.num_threads > 0) {
        omp_set_num_threads(opts.num_threads);
    }
#endif

    const bool parallel_samples =
        opts.num_threads != 1 && sample_groups.size() > 1;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        auto fields = split_tab(line);
        if (fields.size() < 2) continue;
        if (fields.size() < input_header.size()) {
            fields.resize(input_header.size(), ".");
        }

        const std::string& chr = fields[0];
        const int pos = std::stoi(fields[1]);
        const bool new_chromosome = chr != current_chr;
        if (new_chromosome) {
            current_chr = chr;
        }

        if (parallel_samples) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (std::size_t g = 0; g < sample_groups.size(); ++g) {
                process_sample_group(sample_groups[g], targets, opts, fields, pos, new_chromosome,
                                     imputed_values);
            }
        } else {
            for (auto& group : sample_groups) {
                process_sample_group(group, targets, opts, fields, pos, new_chromosome,
                                     imputed_values);
            }
        }

        std::unordered_map<std::string, ImputeResult> imputed;
        imputed.reserve(targets.size());
        for (std::size_t t = 0; t < targets.size(); ++t) {
            imputed[targets[t].y_col] = imputed_values[t];
        }

        write_projected_row_all(out, input_header, fields, plan, imputed);
    }

#ifdef _OPENMP
    if (opts.num_threads > 0) {
        omp_set_num_threads(prev_threads);
    }
#endif
}

void stream_beta_binomial_impute_all(const std::string& input_path,
                                     const std::string& output_path,
                                     const ImputeOptions& opts) {
    std::ifstream in(input_path);
    if (!in) throw std::runtime_error("Cannot open input: " + input_path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty input: " + input_path);

    const std::vector<std::string> input_header = split_tab(line);
    const std::vector<HaplotypeTarget> targets =
        opts.hap_mode ? discover_haplotype_targets(input_header, opts.mode)
                      : discover_sample_targets(input_header, opts.mode);
    if (targets.empty()) {
        throw std::runtime_error(opts.hap_mode
                                     ? "No {sample}.hap{1,2}_counts columns found in header"
                                     : "No {sample}.counts columns found in header");
    }

    stream_beta_binomial_impute_targets(input_path, output_path, targets, opts);
}

}  // namespace impute_methylation
