#pragma once

#include "merge_bedmethyl/locus.hpp"
#include "merge_bedmethyl/reader.hpp"

#include <ostream>

namespace merge_bedmethyl {

void append_haplotype_columns(std::ostream& out, const Reader& hp1, const Reader& hp2,
                              const Locus& target, int min_coverage);

void append_sample_columns(std::ostream& out, const Reader& reader, const Locus& target,
                           int min_coverage);

}  // namespace merge_bedmethyl
