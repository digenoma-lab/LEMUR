#pragma once

#include <string>
#include <vector>

namespace dml {

struct SampleMeta {
    std::string sample_id;
    int phenotype_bin = 0;
    int count_idx = -1;
    int cov_idx = -1;
    std::vector<double> xrow;
};

struct CohortMeta {
    std::vector<SampleMeta> samples;
    int n_params = 0;
    bool sample_mode = false;
};

// Reads metadata CSV and maps count/cov column indices in the methylation TSV.
// sample_mode: {id}.counts_imputed / {id}.cov_imputed (requires --sample).
CohortMeta load_metadata(const std::string& path, const std::vector<std::string>& meth_header,
                         bool sample_mode, const std::string& case_label,
                         const std::string& control_label);

}  // namespace dml
