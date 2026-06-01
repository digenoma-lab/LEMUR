#pragma once

#include "merge_bedmethyl/locus.hpp"

namespace merge_bedmethyl {

struct Record {
    Locus locus;
    int coverage = 0;
    double meth_frac = 0.0;
    bool valid = false;
};

}  // namespace merge_bedmethyl
