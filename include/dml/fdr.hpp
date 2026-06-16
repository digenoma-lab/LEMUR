#pragma once

#include "dml/model.hpp"

#include <vector>

namespace dml {

void apply_fdr_bh(std::vector<SiteResult>& results);

void write_results_csv(const std::string& path, std::vector<SiteResult>& results);

}  // namespace dml
