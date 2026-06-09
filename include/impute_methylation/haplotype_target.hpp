#pragma once

#include <deque>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace impute_methylation {

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
    std::string out_col;
};

std::vector<HaplotypeTarget> discover_haplotype_targets(
    const std::vector<std::string>& header);

std::vector<HaplotypeTarget> discover_sample_targets(
    const std::vector<std::string>& header);

// y_col must be a {sample}.hap{1,2}_counts column name from the header.
HaplotypeTarget find_haplotype_target(const std::vector<std::string>& header,
                                      const std::string& y_col);

struct OutputColumnPlan {
    std::vector<std::string> header;
    std::unordered_set<std::string> drop_cols;
    // At each *_counts column, emit imputed value (key = counts column name).
    std::unordered_set<std::string> impute_at_counts_cols;
};

OutputColumnPlan build_output_plan(const std::vector<std::string>& input_header,
                                   const std::vector<HaplotypeTarget>& targets);

void write_projected_row_all(
    std::ostream& out, const std::vector<std::string>& input_header,
    const std::vector<std::string>& fields, const OutputColumnPlan& plan,
    const std::unordered_map<std::string, std::string>& imputed_by_counts_col);

std::string impute_from_window(const std::deque<WindowSite>& window, const ImputeOptions& opts,
                               double current_y, double current_n, double current_pct);

}  // namespace impute_methylation
