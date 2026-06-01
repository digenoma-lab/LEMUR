#include "merge_bedmethyl/locus.hpp"

namespace merge_bedmethyl {

bool Locus::operator<(const Locus& other) const {
    if (chr != other.chr) return chr < other.chr;
    return pos < other.pos;
}

bool Locus::operator==(const Locus& other) const {
    return chr == other.chr && pos == other.pos;
}

}  // namespace merge_bedmethyl
