#include "dml/smoothing.hpp"

#include "dml/util.hpp"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace dml {

namespace {

double arcsin_methylation(double y, double n) {
    constexpr double c0 = 0.1;
    const double frac = (y + c0) / (n + 2.0 * c0);
    return std::asin(2.0 * frac - 1.0);
}

void smooth_chromosome_moving_avg(std::vector<double>& values, const std::vector<int>& positions,
                                  int span_bp) {
    const int n = static_cast<int>(values.size());
    if (n == 0) return;

    std::vector<int> valid;
    valid.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        if (!is_nan(values[static_cast<std::size_t>(i)])) {
            valid.push_back(i);
        }
    }

    const int nv = static_cast<int>(valid.size());
    std::vector<double> out(static_cast<std::size_t>(n), kNaN);
    const int half_span = span_bp / 2;

    for (int vi = 0; vi < nv; ++vi) {
        const int i = valid[static_cast<std::size_t>(vi)];
        const int pos_i = positions[static_cast<std::size_t>(i)];

        int left = vi;
        while (left > 0 &&
               positions[static_cast<std::size_t>(valid[static_cast<std::size_t>(left - 1)])] >=
                   pos_i - half_span) {
            --left;
        }
        int right = vi;
        while (right + 1 < nv &&
               positions[static_cast<std::size_t>(valid[static_cast<std::size_t>(right + 1)])] <=
                   pos_i + half_span) {
            ++right;
        }

        double sum = 0.0;
        int count = 0;
        for (int k = left; k <= right; ++k) {
            sum += values[static_cast<std::size_t>(valid[static_cast<std::size_t>(k)])];
            ++count;
        }
        out[static_cast<std::size_t>(i)] = std::round(sum / static_cast<double>(count));
    }

    values = std::move(out);
}

}  // namespace

std::vector<SiteInput> load_all_sites(const std::string& path, const CohortMeta& meta) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open methylation TSV: " + path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty methylation TSV: " + path);

    const std::vector<std::string> header = split_csv_line(line, '\t');
    if (column_index(header, "chr") < 0 || column_index(header, "pos") < 0) {
        throw std::runtime_error("Methylation TSV must contain chr and pos columns");
    }

    std::vector<SiteInput> sites;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto fields = split_csv_line(line, '\t');
        if (fields.size() < 2) continue;

        SiteInput site;
        site.chr = fields[0];
        site.pos = std::stoi(fields[1]);
        site.counts.reserve(meta.samples.size());
        site.cov.reserve(meta.samples.size());

        for (const auto& sample : meta.samples) {
            const std::size_t ci = static_cast<std::size_t>(sample.count_idx);
            const std::size_t vi = static_cast<std::size_t>(sample.cov_idx);
            site.counts.push_back(ci < fields.size() ? parse_double_field(fields[ci]) : kNaN);
            site.cov.push_back(vi < fields.size() ? parse_double_field(fields[vi]) : kNaN);
        }
        sites.push_back(std::move(site));
    }

    return sites;
}

void apply_dss_smoothing(std::vector<SiteInput>& sites, int span_bp) {
    if (sites.empty() || span_bp <= 0) return;

    std::unordered_map<std::string, std::vector<std::size_t>> chr_indices;
    chr_indices.reserve(sites.size() / 1000 + 1);
    for (std::size_t i = 0; i < sites.size(); ++i) {
        chr_indices[sites[i].chr].push_back(i);
    }

    const std::size_t num_samples = sites[0].counts.size();
    for (auto& [chr, indices] : chr_indices) {
        (void)chr;
        const int n = static_cast<int>(indices.size());
        std::vector<int> positions(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            positions[static_cast<std::size_t>(i)] = sites[indices[static_cast<std::size_t>(i)]].pos;
        }

        for (std::size_t s = 0; s < num_samples; ++s) {
            std::vector<double> counts(static_cast<std::size_t>(n));
            std::vector<double> cov(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i) {
                const std::size_t idx = indices[static_cast<std::size_t>(i)];
                counts[static_cast<std::size_t>(i)] = sites[idx].counts[s];
                cov[static_cast<std::size_t>(i)] = sites[idx].cov[s];
            }

            smooth_chromosome_moving_avg(counts, positions, span_bp);
            smooth_chromosome_moving_avg(cov, positions, span_bp);

            for (int i = 0; i < n; ++i) {
                const std::size_t idx = indices[static_cast<std::size_t>(i)];
                const double y_raw = sites[idx].counts[s];
                const double n_raw = sites[idx].cov[s];
                const double y_sm = counts[static_cast<std::size_t>(i)];
                const double n_sm = cov[static_cast<std::size_t>(i)];

                if (is_nan(y_raw) || is_nan(n_raw) || n_raw <= 0.0) continue;

                if (sites[idx].z.empty()) {
                    sites[idx].z.assign(num_samples, kNaN);
                }
                if (!is_nan(y_sm) && !is_nan(n_sm) && n_sm > 0.0) {
                    sites[idx].z[s] = arcsin_methylation(y_sm, n_sm);
                }
            }
        }
    }
}

}  // namespace dml
