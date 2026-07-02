#pragma once

#include "dml/model.hpp"

#include <string>
#include <vector>

namespace dml {

struct SmoothingOptions {
    bool enabled = false;
    int span_bp = 500;
};

// Load all CpG sites from a methylation TSV (required before smoothing).
std::vector<SiteInput> load_all_sites(const std::string& path, const CohortMeta& meta);

// DSS-style moving average on counts and cov per sample, per chromosome.
// Non-missing sites are smoothed; missing values stay missing. Z is filled from smoothed values.
void apply_dss_smoothing(std::vector<SiteInput>& sites, int span_bp);

}  // namespace dml
