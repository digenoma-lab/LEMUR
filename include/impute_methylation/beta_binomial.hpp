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
    bool sample_mode = false;  // false = haplotype columns (default); true = {id}.counts/.cov
    ImputeMode mode = ImputeMode::Fraction;
    ProcessOperation operation = ProcessOperation::Impute;
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

// Smooth every {sample}.hap{1,2} found in the header (streaming, one pass).
// Observed sites are beta-binomial smoothed; missing sites are left unchanged.
void stream_beta_binomial_smooth_all(const std::string& input_path,
                                     const std::string& output_path,
                                     const ImputeOptions& opts);

// Smooth only the listed targets; other haplotype columns are copied unchanged.
void stream_beta_binomial_smooth_targets(
    const std::string& input_path, const std::string& output_path,
    const std::vector<HaplotypeTarget>& targets, const ImputeOptions& opts);

}  // namespace impute_methylation
