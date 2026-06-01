#include "merge_bedmethyl/format.hpp"

#include <iomanip>
#include <sstream>

namespace merge_bedmethyl {

std::string format_fraction(double value) {
    constexpr double kEps = 1e-9;
    if (value <= kEps) return "0";
    if (value >= 1.0 - kEps) return "1";
    std::ostringstream os;
    os << std::fixed << std::setprecision(1) << value;
    return os.str();
}

std::string format_pair_cell(const Reader& a, const Reader& b, const Locus& target,
                             int min_coverage) {
    const bool has1 = passes_coverage(a, target, min_coverage);
    const bool has2 = passes_coverage(b, target, min_coverage);

    std::ostringstream os;
    os << (has1 ? format_fraction(a.current.meth_frac) : ".");
    os << '|';
    os << (has2 ? format_fraction(b.current.meth_frac) : ".");
    return os.str();
}

}  // namespace merge_bedmethyl
