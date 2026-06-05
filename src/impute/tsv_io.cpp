#include "impute_methylation/tsv_io.hpp"

#include <sstream>

namespace impute_methylation {

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream iss(line);
    while (std::getline(iss, field, '\t')) fields.push_back(field);
    return fields;
}

int column_index(const std::vector<std::string>& header, const std::string& name) {
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name) return static_cast<int>(i);
    }
    return -1;
}

bool has_column(const std::vector<std::string>& header, const std::string& name) {
    return column_index(header, name) >= 0;
}

double parse_double_field(const std::string& s) {
    if (is_missing(s)) return kNaN;
    try {
        return std::stod(s);
    } catch (...) {
        return kNaN;
    }
}

double parse_fraction_field(const std::string& s) {
    return parse_double_field(s);
}

double ground_truth_fraction(double y, double n, double pct) {
    if (is_nan(n) || n <= 0.0) return kNaN;
    if (!is_nan(y)) return y / n;
    if (!is_nan(pct)) return pct / 100.0;
    return kNaN;
}

void write_header_line(std::ostream& out, const std::vector<std::string>& header) {
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (i) out << '\t';
        out << header[i];
    }
    out << '\n';
}

}  // namespace impute_methylation
