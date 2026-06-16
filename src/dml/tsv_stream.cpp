#include "dml/tsv_stream.hpp"

#include "dml/util.hpp"

#include <stdexcept>

namespace dml {

std::vector<std::string> read_tsv_header(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open methylation TSV: " + path);
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty methylation TSV: " + path);
    return split_csv_line(line, '\t');
}

MethylationStream::MethylationStream(const std::string& path, const CohortMeta& meta)
    : in_(path), meta_(meta) {
    if (!in_) throw std::runtime_error("Cannot open methylation TSV: " + path);
}

bool MethylationStream::read_batch(std::vector<SiteInput>& batch, std::size_t max_rows) {
    batch.clear();
    batch.reserve(max_rows);

    if (!header_read_) {
        std::string line;
        if (!std::getline(in_, line)) return false;
        header_ = split_csv_line(line, '\t');
        if (column_index(header_, "chr") < 0 || column_index(header_, "pos") < 0) {
            throw std::runtime_error("Methylation TSV must contain chr and pos columns");
        }
        header_read_ = true;
    }

    std::string line;
    while (batch.size() < max_rows && std::getline(in_, line)) {
        if (line.empty()) continue;
        auto fields = split_csv_line(line, '\t');
        if (fields.size() < 2) continue;

        SiteInput site;
        site.chr = fields[0];
        site.pos = std::stoi(fields[1]);
        site.counts.reserve(meta_.samples.size());
        site.cov.reserve(meta_.samples.size());

        for (const auto& sample : meta_.samples) {
            const std::size_t ci = static_cast<std::size_t>(sample.count_idx);
            const std::size_t vi = static_cast<std::size_t>(sample.cov_idx);
            site.counts.push_back(ci < fields.size() ? parse_double_field(fields[ci]) : kNaN);
            site.cov.push_back(vi < fields.size() ? parse_double_field(fields[vi]) : kNaN);
        }
        batch.push_back(std::move(site));
    }

    return !batch.empty();
}

}  // namespace dml
