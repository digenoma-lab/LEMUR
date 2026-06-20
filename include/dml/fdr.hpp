#pragma once

#include "dml/model.hpp"

#include <vector>

namespace dml {

// DSS/BSseq report 1-based CpG positions (BED start + 1); see DifferentialMethylationRegions dss.R.
constexpr int kDssOutputCoordOffset = 1;

void apply_fdr_bh(std::vector<SiteResult>& results);

void write_results_csv(const std::string& path, std::vector<SiteResult>& results);

}  // namespace dml
