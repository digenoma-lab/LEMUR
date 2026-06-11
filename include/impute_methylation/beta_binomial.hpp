#pragma once

#include "impute_methylation/haplotype_target.hpp"

#include <string>
#include <vector>

namespace impute_methylation {

struct ImputeOptions {
    int window_bp = 200;
    double alpha = 1.0;
    double beta = 1.0;
    int min_neighbors = 5;
    int num_threads = 1;  // 0 = all available cores (OpenMP)
    bool hap_mode = false;
    ImputeMode mode = ImputeMode::Fraction;
};

std::string format_fraction(double value);

// Impute every {sample}.hap{1,2} found in the header (streaming, one pass).
void stream_beta_binomial_impute_all(const std::string& input_path,
                                     const std::string& output_path,
                                     const ImputeOptions& opts);

// Impute only the listed targets; other haplotype columns are copied unchanged.
void stream_beta_binomial_impute_targets(
    const std::string& input_path, const std::string& output_path,
    const std::vector<HaplotypeTarget>& targets, const ImputeOptions& opts);

}  // namespace impute_methylation
