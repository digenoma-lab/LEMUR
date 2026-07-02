#pragma once

#include "dml/metadata.hpp"
#include "dml/util.hpp"

#include <optional>
#include <string>
#include <vector>

namespace dml {

struct SiteInput {
    std::string chr;
    int pos = 0;
    std::vector<double> counts;
    std::vector<double> cov;
    // Precomputed arcsin(Z) from smoothed counts/cov (DSS smoothing); empty = use raw counts/cov.
    std::vector<double> z;
};

struct SiteResult {
    std::string chr;
    int pos = 0;
    double beta_phenotype = kNaN;
    double se_phenotype = kNaN;
    double pvalue = kNaN;
    double phi = kNaN;
    double mean_case = kNaN;
    double mean_control = kNaN;
    double delta_beta = kNaN;
    int n_samples = 0;
    double fdr = kNaN;
    bool significant = false;
};

// Fit DSS-style DML (arcsin WLS) for one CpG; returns nullopt if not enough data or fit fails.
std::optional<SiteResult> fit_site(const SiteInput& site, const CohortMeta& meta);

}  // namespace dml
