#pragma once

#include "impute_methylation/beta_binomial.hpp"

#include <cstdint>
#include <string>

namespace impute_methylation {

struct EvaluateOptions {
    std::string y_col;       // e.g. CHI08A.hap1_counts (required)
    std::string chromosome = "chr1";  // only rows on this chr are masked/imputed/scored
    double mask_fraction = 0.2;
    uint64_t seed = 42;
    ImputeOptions impute;
    std::string masked_path;   // empty → input + ".masked.tsv"
    std::string imputed_path;  // empty → input + ".imputed.eval.tsv"
};

struct EvaluateMetrics {
    std::size_t n_masked = 0;
    std::size_t n_scored = 0;
    double mse = 0.0;
    double pearson = 0.0;
};

// Mask fraction of valid sites in y_col, impute that haplotype, score vs ground truth.
EvaluateMetrics run_evaluate(const std::string& input_path, const EvaluateOptions& opts);

}  // namespace impute_methylation
