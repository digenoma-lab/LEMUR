#include "dml/util.hpp"

#include <sstream>
#include <stdexcept>

namespace dml {

std::vector<std::string> split_csv_line(const std::string& line, char delim) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream iss(line);
    while (std::getline(iss, field, delim)) {
        fields.push_back(field);
    }
    return fields;
}

int column_index(const std::vector<std::string>& header, const std::string& name) {
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name) return static_cast<int>(i);
    }
    return -1;
}

double parse_double_field(const std::string& s) {
    if (s.empty() || s == ".") return kNaN;
    try {
        return std::stod(s);
    } catch (...) {
        return kNaN;
    }
}

double sigmoid(double x) {
    if (x >= 0) {
        const double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    }
    const double z = std::exp(x);
    return z / (1.0 + z);
}

double normal_sf(double z) {
    // Abramowitz & Stegun 26.2.17 approximation for Phi(z).
    const double abs_z = std::abs(z);
    const double t = 1.0 / (1.0 + 0.2316419 * abs_z);
    const double d = 0.3989423 * std::exp(-0.5 * abs_z * abs_z);
    const double p = d * t *
                     (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
    return z >= 0 ? p : 1.0 - p;
}

double digamma(double x) {
    if (x <= 0.0) return kNaN;
    double result = 0.0;
    while (x < 8.0) {
        result -= 1.0 / x;
        x += 1.0;
    }
    const double inv = 1.0 / x;
    const double inv2 = inv * inv;
    result += std::log(x) - 0.5 * inv - inv2 / 12.0 + inv2 * inv2 / 120.0;
    return result;
}

double trigamma(double x) {
    if (x <= 0.0) return kNaN;
    double result = 0.0;
    while (x < 8.0) {
        result += 1.0 / (x * x);
        x += 1.0;
    }
    const double inv = 1.0 / x;
    const double inv2 = inv * inv;
    result += 0.5 * inv2 + inv2 * inv / 6.0 - inv2 * inv2 / 30.0;
    return result;
}

}  // namespace dml
