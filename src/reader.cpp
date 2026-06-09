#include "merge_bedmethyl/bedmethyl_columns.hpp"
#include "merge_bedmethyl/reader.hpp"

#include <sstream>
#include <utility>

namespace merge_bedmethyl {

Reader::Reader(std::string path, std::string sample_label)
    : in(std::move(path)), label(std::move(sample_label)) {
    advance();
}

bool Reader::eof() const { return !current.valid; }

void Reader::advance() {
    current.valid = false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (parse_line(line, current)) return;
    }
}

bool Reader::parse_line(const std::string& line, Record& rec) {
    std::vector<std::string> fields;
    fields.reserve(18);
    std::istringstream iss(line);
    std::string field;
    while (std::getline(iss, field, '\t')) fields.push_back(field);

    if (static_cast<int>(fields.size()) <= kNModCol) return false;

    try {
        rec.locus.chr = fields[kChrCol];
        rec.locus.pos = std::stoi(fields[kStartCol]);
        rec.coverage = std::stoi(fields[kCovCol]);
        rec.meth_percent = std::stod(fields[kMethCol]);
        rec.n_modified = std::stoi(fields[kNModCol]);
        rec.valid = true;
        return true;
    } catch (...) {
        return false;
    }
}

PairReaders::PairReaders(std::string sample_label, std::string path1, std::string path2)
    : label(std::move(sample_label)),
      hp1(std::move(path1), label + "_hp1"),
      hp2(std::move(path2), label + "_hp2") {}

SingleReader::SingleReader(std::string sample_label, std::string path)
    : label(std::move(sample_label)), reader(std::move(path), label) {}

void skip_until(Reader& r, const Locus& target) {
    while (!r.eof() && r.current.locus < target) r.advance();
}

Locus min_locus(const std::vector<PairReaders>& pairs) {
    Locus min_l;
    bool first = true;
    for (const auto& pair : pairs) {
        for (const Reader* r : {&pair.hp1, &pair.hp2}) {
            if (r->eof()) continue;
            if (first || r->current.locus < min_l) {
                min_l = r->current.locus;
                first = false;
            }
        }
    }
    return min_l;
}

Locus min_locus(const std::vector<SingleReader>& samples) {
    Locus min_l;
    bool first = true;
    for (const auto& sample : samples) {
        if (sample.reader.eof()) continue;
        if (first || sample.reader.current.locus < min_l) {
            min_l = sample.reader.current.locus;
            first = false;
        }
    }
    return min_l;
}

bool any_valid(const std::vector<PairReaders>& pairs) {
    for (const auto& pair : pairs) {
        if (!pair.hp1.eof() || !pair.hp2.eof()) return true;
    }
    return false;
}

bool any_valid(const std::vector<SingleReader>& samples) {
    for (const auto& sample : samples) {
        if (!sample.reader.eof()) return true;
    }
    return false;
}

void advance_readers_at_locus(std::vector<PairReaders>& pairs, const Locus& target) {
    for (auto& pair : pairs) {
        if (!pair.hp1.eof() && pair.hp1.current.locus == target) pair.hp1.advance();
        if (!pair.hp2.eof() && pair.hp2.current.locus == target) pair.hp2.advance();
    }
}

void advance_readers_at_locus(std::vector<SingleReader>& samples, const Locus& target) {
    for (auto& sample : samples) {
        if (!sample.reader.eof() && sample.reader.current.locus == target) {
            sample.reader.advance();
        }
    }
}

bool passes_coverage(const Reader& r, const Locus& target, int min_coverage) {
    return !r.eof() && r.current.locus == target && r.current.coverage > min_coverage;
}

bool sample_has_information(const PairReaders& pair, const Locus& target, int min_coverage) {
    return passes_coverage(pair.hp1, target, min_coverage) ||
           passes_coverage(pair.hp2, target, min_coverage);
}

bool sample_has_information(const SingleReader& sample, const Locus& target, int min_coverage) {
    return passes_coverage(sample.reader, target, min_coverage);
}

int count_samples_with_information(const std::vector<PairReaders>& pairs, const Locus& target,
                                   int min_coverage) {
    int count = 0;
    for (const auto& pair : pairs) {
        if (sample_has_information(pair, target, min_coverage)) ++count;
    }
    return count;
}

int count_samples_with_information(const std::vector<SingleReader>& samples,
                                   const Locus& target, int min_coverage) {
    int count = 0;
    for (const auto& sample : samples) {
        if (sample_has_information(sample, target, min_coverage)) ++count;
    }
    return count;
}

}  // namespace merge_bedmethyl
