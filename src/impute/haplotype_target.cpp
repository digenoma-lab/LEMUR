#include "impute_methylation/haplotype_target.hpp"

#include "impute_methylation/beta_binomial.hpp"
#include "impute_methylation/tsv_io.hpp"

#include <ostream>

#include <cmath>
#include <stdexcept>

namespace impute_methylation {

namespace {

bool parse_sample_counts_column(const std::string& col, std::string& sample_id) {
    const std::string suffix = ".counts";
    if (col.size() <= suffix.size()) return false;
    if (col.compare(col.size() - suffix.size(), suffix.size(), suffix) != 0) return false;
    sample_id = col.substr(0, col.size() - suffix.size());
    return !sample_id.empty();
}

bool parse_counts_column(const std::string& col, std::string& sample_id, std::string& hap) {
    const std::string suffix = "_counts";
    if (col.size() <= suffix.size()) return false;
    if (col.compare(col.size() - suffix.size(), suffix.size(), suffix) != 0) return false;

    const std::size_t h1 = col.rfind(".hap1_counts");
    const std::size_t h2 = col.rfind(".hap2_counts");
    if (h1 != std::string::npos && h1 + 12 == col.size()) {
        sample_id = col.substr(0, h1);
        hap = "hap1";
        return !sample_id.empty();
    }
    if (h2 != std::string::npos && h2 + 12 == col.size()) {
        sample_id = col.substr(0, h2);
        hap = "hap2";
        return !sample_id.empty();
    }
    return false;
}

bool is_dropped(const std::unordered_set<std::string>& drop, const std::string& col) {
    return drop.count(col) > 0;
}

void set_fraction_outputs(HaplotypeTarget& t) {
    if (!t.hap.empty()) {
        t.out_col = t.sample_id + "." + t.hap + "_frac_imputed";
    } else {
        t.out_col = t.sample_id + ".frac_imputed";
    }
}

void set_counts_cov_outputs(HaplotypeTarget& t) {
    if (!t.hap.empty()) {
        t.out_y_col = t.sample_id + "." + t.hap + "_counts_imputed";
        t.out_n_col = t.sample_id + "." + t.hap + "_cov_imputed";
    } else {
        t.out_y_col = t.sample_id + ".counts_imputed";
        t.out_n_col = t.sample_id + ".cov_imputed";
    }
}

std::string format_count(double value) {
    if (is_nan(value)) return ".";
    const auto rounded = static_cast<long long>(std::llround(value));
    if (rounded < 0) return "0";
    return std::to_string(rounded);
}

std::string format_cov(double value) {
    if (is_nan(value)) return ".";
    const auto rounded = static_cast<long long>(std::llround(value));
    if (rounded <= 0) return ".";
    return std::to_string(rounded);
}

std::string observed_fraction(double y, double n, double pct) {
    const double frac = ground_truth_fraction(y, n, pct);
    if (is_nan(frac)) return ".";
    return format_fraction(frac);
}

ImputeResult observed_values(double y, double n, double pct, ImputeMode mode) {
    ImputeResult result;
    if (mode == ImputeMode::Fraction) {
        result.fraction = observed_fraction(y, n, pct);
        return result;
    }

    if (is_nan(n) || n <= 0.0) {
        result.counts = ".";
        result.cov = ".";
        return result;
    }

    result.cov = format_cov(n);
    if (!is_nan(y)) {
        result.counts = format_count(y);
    } else if (!is_nan(pct)) {
        result.counts = format_count((pct / 100.0) * n);
    } else {
        result.counts = ".";
    }
    return result;
}

}  // namespace

std::vector<HaplotypeTarget> discover_haplotype_targets(const std::vector<std::string>& header,
                                                        ImputeMode mode) {
    std::vector<HaplotypeTarget> targets;
    targets.reserve(header.size() / 6);

    for (std::size_t i = 0; i < header.size(); ++i) {
        std::string sample_id;
        std::string hap;
        if (!parse_counts_column(header[i], sample_id, hap)) continue;

        const std::string n_col = sample_id + "." + hap + "_cov";
        const int n_idx = column_index(header, n_col);
        if (n_idx < 0) continue;

        HaplotypeTarget t;
        t.sample_id = sample_id;
        t.hap = hap;
        t.y_idx = static_cast<int>(i);
        t.n_idx = n_idx;
        t.y_col = header[i];
        t.n_col = n_col;
        t.pct_col = sample_id + "." + hap + "_percentage";
        t.pct_idx = column_index(header, t.pct_col);
        if (mode == ImputeMode::CountsCov) {
            set_counts_cov_outputs(t);
        } else {
            set_fraction_outputs(t);
        }
        targets.push_back(std::move(t));
    }
    return targets;
}

std::vector<HaplotypeTarget> discover_sample_targets(const std::vector<std::string>& header,
                                                     ImputeMode mode) {
    std::vector<HaplotypeTarget> targets;
    targets.reserve(header.size() / 3);

    for (std::size_t i = 0; i < header.size(); ++i) {
        std::string sample_id;
        if (!parse_sample_counts_column(header[i], sample_id)) continue;

        const std::string n_col = sample_id + ".cov";
        const int n_idx = column_index(header, n_col);
        if (n_idx < 0) continue;

        HaplotypeTarget t;
        t.sample_id = sample_id;
        t.y_idx = static_cast<int>(i);
        t.n_idx = n_idx;
        t.y_col = header[i];
        t.n_col = n_col;
        t.pct_col = sample_id + ".percentage";
        t.pct_idx = column_index(header, t.pct_col);
        if (mode == ImputeMode::CountsCov) {
            set_counts_cov_outputs(t);
        } else {
            set_fraction_outputs(t);
        }
        targets.push_back(std::move(t));
    }
    return targets;
}

HaplotypeTarget find_haplotype_target(const std::vector<std::string>& header,
                                      const std::string& y_col) {
    for (const auto& t : discover_haplotype_targets(header)) {
        if (t.y_col == y_col) return t;
    }
    throw std::runtime_error("Counts column not found in header: " + y_col);
}

HaplotypeTarget find_sample_target(const std::vector<std::string>& header,
                                   const std::string& y_col) {
    for (const auto& t : discover_sample_targets(header)) {
        if (t.y_col == y_col) return t;
    }
    throw std::runtime_error("Counts column not found in header: " + y_col);
}

OutputColumnPlan build_output_plan(const std::vector<std::string>& input_header,
                                   const std::vector<HaplotypeTarget>& targets,
                                   ImputeMode mode) {
    OutputColumnPlan plan;
    plan.mode = mode;
    for (const auto& t : targets) {
        plan.drop_cols.insert(t.y_col);
        plan.drop_cols.insert(t.n_col);
        plan.drop_cols.insert(t.pct_col);
        plan.drop_cols.insert(t.out_col);
        plan.drop_cols.insert(t.out_y_col);
        plan.drop_cols.insert(t.out_n_col);
        plan.impute_at_counts_cols.insert(t.y_col);
    }

    for (const auto& col : input_header) {
        if (is_dropped(plan.drop_cols, col)) {
            if (plan.impute_at_counts_cols.count(col) > 0) {
                for (const auto& t : targets) {
                    if (t.y_col != col) continue;
                    if (mode == ImputeMode::CountsCov) {
                        plan.header.push_back(t.out_y_col);
                        plan.header.push_back(t.out_n_col);
                    } else {
                        plan.header.push_back(t.out_col);
                    }
                    break;
                }
            }
            continue;
        }
        plan.header.push_back(col);
    }
    return plan;
}

void write_projected_row_all(std::ostream& out, const std::vector<std::string>& input_header,
                             const std::vector<std::string>& fields,
                             const OutputColumnPlan& plan,
                             const std::unordered_map<std::string, ImputeResult>& imputed) {
    bool first = true;
    for (std::size_t i = 0; i < input_header.size(); ++i) {
        const std::string& col = input_header[i];
        if (is_dropped(plan.drop_cols, col)) {
            if (plan.impute_at_counts_cols.count(col) > 0) {
                auto it = imputed.find(col);
                if (it != imputed.end()) {
                    const ImputeResult& values = it->second;
                    if (plan.mode == ImputeMode::CountsCov) {
                        if (!first) out << '\t';
                        out << (values.counts.empty() ? "." : values.counts);
                        first = false;
                        out << '\t';
                        out << (values.cov.empty() ? "." : values.cov);
                    } else {
                        if (!first) out << '\t';
                        out << (values.fraction.empty() ? "." : values.fraction);
                        first = false;
                    }
                } else {
                    if (!first) out << '\t';
                    out << '.';
                    first = false;
                }
            }
            continue;
        }
        if (!first) out << '\t';
        out << (i < fields.size() ? fields[i] : ".");
        first = false;
    }
    out << '\n';
}

ImputeResult impute_from_window(const std::deque<WindowSite>& window, const ImputeOptions& opts,
                                double current_y, double current_n, double current_pct) {
    double y_sum = 0.0;
    double n_sum = 0.0;
    std::size_t valid = 0;

    for (const auto& site : window) {
        if (is_nan(site.n) || site.n <= 0.0) continue;
        ++valid;
        n_sum += site.n;
        if (!is_nan(site.y)) y_sum += site.y;
    }

    if (valid < static_cast<std::size_t>(opts.min_neighbors)) {
        return observed_values(current_y, current_n, current_pct, opts.mode);
    }

    const double p_hat = (y_sum + opts.alpha) / (n_sum + opts.alpha + opts.beta);
    ImputeResult result;

    if (opts.mode == ImputeMode::CountsCov) {
        const double cov_hat = n_sum / static_cast<double>(valid);
        double counts_hat = std::llround(p_hat * cov_hat);
        if (counts_hat < 0.0) counts_hat = 0.0;
        if (counts_hat > cov_hat) counts_hat = std::floor(cov_hat);
        result.cov = format_cov(cov_hat);
        result.counts = format_count(counts_hat);
    } else {
        result.fraction = format_fraction(p_hat);
    }

    return result;
}

}  // namespace impute_methylation
