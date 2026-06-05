#include "impute_methylation/beta_binomial.hpp"

#include "impute_methylation/haplotype_target.hpp"
#include "impute_methylation/tsv_io.hpp"

#include <deque>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace impute_methylation {

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

    const OutputColumnPlan plan = build_output_plan(input_header, targets);
    write_header_line(out, plan.header);

    std::vector<std::deque<WindowSite>> windows(targets.size());
    std::string current_chr;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        auto fields = split_tab(line);
        if (fields.size() < 2) continue;
        if (fields.size() < input_header.size()) {
            fields.resize(input_header.size(), ".");
        }

        const std::string& chr = fields[0];
        const int pos = std::stoi(fields[1]);

        if (chr != current_chr) {
            for (auto& w : windows) w.clear();
            current_chr = chr;
        }

        std::unordered_map<std::string, std::string> imputed;
        imputed.reserve(targets.size());

        for (std::size_t t = 0; t < targets.size(); ++t) {
            const auto& target = targets[t];
            auto& window = windows[t];

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
                pct = parse_double_field(
                    fields[static_cast<std::size_t>(target.pct_idx)]);
            }
            imputed[target.y_col] =
                impute_from_window(window, opts, site.y, site.n, pct);
        }

        write_projected_row_all(out, input_header, fields, plan, imputed);
    }
}

void stream_beta_binomial_impute_all(const std::string& input_path,
                                     const std::string& output_path,
                                     const ImputeOptions& opts) {
    std::ifstream in(input_path);
    if (!in) throw std::runtime_error("Cannot open input: " + input_path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty input: " + input_path);

    const std::vector<std::string> input_header = split_tab(line);
    const std::vector<HaplotypeTarget> targets = discover_haplotype_targets(input_header);
    if (targets.empty()) {
        throw std::runtime_error("No {sample}.hap{1,2}_counts columns found in header");
    }

    stream_beta_binomial_impute_targets(input_path, output_path, targets, opts);
}

}  // namespace impute_methylation
