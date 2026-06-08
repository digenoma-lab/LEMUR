#pragma once

#include "impute_methylation/beta_binomial.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace impute_methylation {

struct EvaluateOptions {
    std::string y_col;       // e.g. CHI08A.hap1_counts (required for single-target mode)
    std::string chromosome = "chr1";  // only rows on this chr are masked/imputed/scored
    double mask_fraction = 0.2;
    uint64_t seed = 42;
    ImputeOptions impute;
    std::string masked_path;   // empty → input + ".masked.tsv"
    std::string imputed_path;  // empty → input + ".imputed.eval.tsv"
    bool quiet = false;
};

struct EvaluateMetrics {
    std::size_t n_masked = 0;
    std::size_t n_scored = 0;
    double mse = 0.0;
    double pearson = 0.0;
};

struct EvaluateCohortOptions {
    std::string chromosome = "chr1";
    double mask_fraction = 0.2;
    uint64_t seed = 42;
    ImputeOptions impute;
    std::string output_path;
    bool keep_sidecars = false;
};

struct EvaluateCohortRow {
    std::string sample;
    int window_size = 0;
    int n_neighbors = 0;
    double a = 1.0;
    double b = 1.0;
    std::string chr;
    double pearson = 0.0;
    double mse = 0.0;
    std::size_t count_masked = 0;
    std::size_t count_imputed = 0;
};

// Mask fraction of valid sites in y_col, impute that haplotype, score vs ground truth.
EvaluateMetrics run_evaluate(const std::string& input_path, const EvaluateOptions& opts);

// Evaluate every {sample}.hap{1,2} in the header (parallel over targets when OpenMP enabled).
std::vector<EvaluateCohortRow> run_evaluate_cohort(const std::string& input_path,
                                                   const EvaluateCohortOptions& opts);

void write_evaluate_cohort_tsv(const std::string& path,
                               const std::vector<EvaluateCohortRow>& rows);

}  // namespace impute_methylation
