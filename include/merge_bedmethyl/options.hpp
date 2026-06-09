#pragma once

#include "impute_methylation/beta_binomial.hpp"

#include <string>
#include <vector>

namespace merge_bedmethyl {

struct SamplePair {
    std::string label;
    std::string hp1_path;
    std::string hp2_path;
};

struct Options {
    int min_coverage = 3;
    int min_samples = -1;  // unset: use N-1
    bool impute = false;
    impute_methylation::ImputeOptions impute_options;
};

struct ParsedArgs {
    Options options;
    std::string output_path;
    std::vector<SamplePair> samples;
    bool show_help = false;
};

// Returns exit code: 0 ok, 1 error, 2 help shown
int parse_arguments(int argc, char* argv[], ParsedArgs& out);
void print_usage(const char* prog);

}  // namespace merge_bedmethyl
