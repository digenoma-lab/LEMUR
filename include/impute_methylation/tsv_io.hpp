#pragma once

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace impute_methylation {

inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
inline bool is_nan(double x) { return std::isnan(x); }
inline bool is_missing(const std::string& s) { return s.empty() || s == "."; }

std::vector<std::string> split_tab(const std::string& line);

int column_index(const std::vector<std::string>& header, const std::string& name);

bool has_column(const std::vector<std::string>& header, const std::string& name);

double parse_double_field(const std::string& s);

double parse_fraction_field(const std::string& s);

double ground_truth_fraction(double y, double n, double pct);

void write_header_line(std::ostream& out, const std::vector<std::string>& header);

}  // namespace impute_methylation
