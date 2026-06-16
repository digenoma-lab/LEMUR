#pragma once

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace dml {

inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
inline bool is_nan(double x) { return std::isnan(x); }

std::vector<std::string> split_csv_line(const std::string& line, char delim);

int column_index(const std::vector<std::string>& header, const std::string& name);

double parse_double_field(const std::string& s);

double sigmoid(double x);

// Standard normal survival function P(Z > z).
double normal_sf(double z);

// Digamma function psi(x) = d/dx log Gamma(x).
double digamma(double x);

// Trigamma function psi'(x).
double trigamma(double x);

}  // namespace dml
