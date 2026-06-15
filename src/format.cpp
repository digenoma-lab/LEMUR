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

void append_sample_columns(std::ostream& out, const Reader& hp1, const Reader& hp2,
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

void append_sample_columns_aggregated(std::ostream& out, const Reader& hp1, const Reader& hp2,
                                      const Locus& target, int min_coverage) {
    const bool has1 = passes_coverage(hp1, target, min_coverage);
    const bool has2 = passes_coverage(hp2, target, min_coverage);

    if (!has1 && !has2) {
        out << "\t.\t.\t.";
        return;
    }

    long long counts = 0;
    long long cov = 0;
    if (has1) {
        counts += hp1.current.n_modified;
        cov += hp1.current.coverage;
    }
    if (has2) {
        counts += hp2.current.n_modified;
        cov += hp2.current.coverage;
    }

    out << '\t' << counts << '\t' << cov;
    const double pct = (cov > 0) ? (100.0 * static_cast<double>(counts) / static_cast<double>(cov))
                                 : 0.0;
    append_percent(out, true, pct);
}

}  // namespace merge_bedmethyl
