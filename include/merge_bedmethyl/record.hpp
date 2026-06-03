#pragma once

#include "merge_bedmethyl/locus.hpp"

namespace merge_bedmethyl {

struct Record {
    Locus locus;
    int coverage = 0;
    int n_modified = 0;
    double meth_percent = 0.0;
    bool valid = false;
};

}  // namespace merge_bedmethyl
