#pragma once

#include <string>
#include <vector>

namespace dml {

struct CallDmrOptions {
    double p_threshold = 1e-5;
    int dis_merge = 100;
    int min_cg = 3;
    int min_len = 50;
    double pct_sig = 0.5;
    int num_threads = 1;  // 0 = all cores (OpenMP)
    bool sample_mode = false;
};

struct DmrRecord {
    std::string chr;
    int start = 0;
    int end = 0;
    int length = 0;
    int n_cg = 0;
    int n_cg_sig = 0;
    double pct_sig = 0.0;
    double area_stat = 0.0;
    double mean_stat = 0.0;
    double mean_diff = 0.0;
    std::string direction;
    double min_p = 0.0;
    double mean_p = 0.0;
};

// DSS callDMR on sorted DML CSV (chr, pos). Streams by chromosome; parallel per chr.
std::vector<DmrRecord> call_dmr_file(const std::string& dml_csv, const CallDmrOptions& opts);

void write_dmr_csv(const std::string& path, const std::vector<DmrRecord>& dmrs);

}  // namespace dml
