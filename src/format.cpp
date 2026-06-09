#include "merge_bedmethyl/format.hpp"

#include <iomanip>
#include <sstream>

namespace merge_bedmethyl {

namespace {

template <typename T>
void append_field(std::ostream& out, bool has, T value) {
    out << '\t';
    if (!has) {
        out << '.';
        return;
    }
    out << value;
}

void append_percent(std::ostream& out, bool has, double percent) {
    out << '\t';
    if (!has) {
        out << '.';
        return;
    }
    constexpr double kEps = 1e-9;
    if (percent <= kEps) {
        out << "0";
        return;
    }
    if (percent >= 100.0 - kEps) {
        out << "100";
        return;
    }
    std::ostringstream os;
    os << std::fixed << std::setprecision(2) << percent;
    std::string s = os.str();
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    out << s;
}

}  // namespace

void append_haplotype_columns(std::ostream& out, const Reader& hp1, const Reader& hp2,
                              const Locus& target, int min_coverage) {
    const bool has1 = passes_coverage(hp1, target, min_coverage);
    const bool has2 = passes_coverage(hp2, target, min_coverage);

    append_field(out, has1, hp1.current.n_modified);
    append_field(out, has2, hp2.current.n_modified);
    append_field(out, has1, hp1.current.coverage);
    append_field(out, has2, hp2.current.coverage);
    append_percent(out, has1, hp1.current.meth_percent);
    append_percent(out, has2, hp2.current.meth_percent);
}

void append_sample_columns(std::ostream& out, const Reader& reader, const Locus& target,
                           int min_coverage) {
    const bool has = passes_coverage(reader, target, min_coverage);

    append_field(out, has, reader.current.n_modified);
    append_field(out, has, reader.current.coverage);
    append_percent(out, has, reader.current.meth_percent);
}

}  // namespace merge_bedmethyl
