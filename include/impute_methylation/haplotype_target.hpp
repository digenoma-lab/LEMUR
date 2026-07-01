#pragma once

#include <deque>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace impute_methylation {

enum class ImputeMode {
    Fraction,   // output {id}.frac_imputed or {id}.frac_smoothed (default)
    CountsCov,  // output {id}.counts and {id}.cov (sample) or {id}.hap{1,2}_counts/_cov
};

enum class ProcessOperation {
    Impute,
    Smooth,
};

struct ImputeOptions;

struct WindowSite {
    int pos = 0;
    double y = 0;
    double n = 0;
};

struct HaplotypeTarget {
    std::string sample_id;
    std::string hap;
    int y_idx = -1;
    int n_idx = -1;
    int pct_idx = -1;
    std::string y_col;
    std::string n_col;
    std::string pct_col;
    std::string out_col;      // fraction mode
    std::string out_y_col;    // counts/cov mode
    std::string out_n_col;    // counts/cov mode
};

struct ImputeResult {
    std::string fraction;  // fraction mode
    std::string counts;    // counts/cov mode
    std::string cov;
};

std::vector<HaplotypeTarget> discover_haplotype_targets(
    const std::vector<std::string>& header, ImputeMode mode = ImputeMode::Fraction,
    ProcessOperation operation = ProcessOperation::Impute);

std::vector<HaplotypeTarget> discover_sample_targets(
    const std::vector<std::string>& header, ImputeMode mode = ImputeMode::Fraction,
    ProcessOperation operation = ProcessOperation::Impute);

// y_col must be a {sample}.hap{1,2}_counts column name from the header.
HaplotypeTarget find_haplotype_target(const std::vector<std::string>& header,
                                      const std::string& y_col,
                                      ImputeMode mode = ImputeMode::Fraction);

// y_col must be a {sample}.counts column name from the header.
HaplotypeTarget find_sample_target(const std::vector<std::string>& header,
                                   const std::string& y_col,
                                   ImputeMode mode = ImputeMode::Fraction);

struct OutputColumnPlan {
    ImputeMode mode = ImputeMode::Fraction;
    std::vector<std::string> header;
    std::unordered_set<std::string> drop_cols;
    // At each *_counts column, emit imputed value (key = counts column name).
    std::unordered_set<std::string> impute_at_counts_cols;
};

OutputColumnPlan build_output_plan(const std::vector<std::string>& input_header,
                                   const std::vector<HaplotypeTarget>& targets,
                                   ImputeMode mode);

void write_projected_row_all(
    std::ostream& out, const std::vector<std::string>& input_header,
    const std::vector<std::string>& fields, const OutputColumnPlan& plan,
    const std::unordered_map<std::string, ImputeResult>& imputed_by_counts_col);

ImputeResult impute_from_window(const std::deque<WindowSite>& window, const ImputeOptions& opts,
                                double current_y, double current_n, double current_pct);

ImputeResult smooth_from_window(const std::deque<WindowSite>& window, const ImputeOptions& opts,
                                double current_y, double current_n, double current_pct,
                                int current_pos);

}  // namespace impute_methylation
