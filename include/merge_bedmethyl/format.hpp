#pragma once

#include "merge_bedmethyl/locus.hpp"
#include "merge_bedmethyl/reader.hpp"

#include <string>

namespace merge_bedmethyl {

std::string format_fraction(double value);
std::string format_pair_cell(const Reader& a, const Reader& b, const Locus& target,
                             int min_coverage);

}  // namespace merge_bedmethyl
