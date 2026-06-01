#pragma once

#include <string>

namespace merge_bedmethyl {

struct Locus {
    std::string chr;
    int pos = 0;

    bool operator<(const Locus& other) const;
    bool operator==(const Locus& other) const;
};

}  // namespace merge_bedmethyl
