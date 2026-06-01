#pragma once

#include "merge_bedmethyl/options.hpp"

namespace merge_bedmethyl {

// Returns process exit code (0 success, 1 error).
int run_merge(const ParsedArgs& args);

}  // namespace merge_bedmethyl
