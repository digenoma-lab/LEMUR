#pragma once

#include "merge_bedmethyl/record.hpp"

#include <fstream>
#include <string>
#include <vector>

namespace merge_bedmethyl {

struct Reader {
    std::ifstream in;
    Record current;
    std::string label;

    Reader(std::string path, std::string sample_label);

    bool eof() const;
    void advance();
    static bool parse_line(const std::string& line, Record& rec);
};

struct PairReaders {
    std::string label;
    Reader hp1;
    Reader hp2;

    PairReaders(std::string sample_label, std::string path1, std::string path2);
};

void skip_until(Reader& r, const Locus& target);
Locus min_locus(const std::vector<PairReaders>& pairs);
bool any_valid(const std::vector<PairReaders>& pairs);
void advance_readers_at_locus(std::vector<PairReaders>& pairs, const Locus& target);

bool passes_coverage(const Reader& r, const Locus& target, int min_coverage);
bool sample_has_information(const PairReaders& pair, const Locus& target, int min_coverage);
int count_samples_with_information(const std::vector<PairReaders>& pairs, const Locus& target,
                                   int min_coverage);

}  // namespace merge_bedmethyl
