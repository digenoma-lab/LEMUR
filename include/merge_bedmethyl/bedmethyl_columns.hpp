#pragma once

namespace merge_bedmethyl {

// modkit bedMethyl column indices (0-based)
inline constexpr int kChrCol = 0;
inline constexpr int kStartCol = 1;
inline constexpr int kCovCol = 9;    // valid coverage
inline constexpr int kMethCol = 10;  // percent_modified (0–100)

}  // namespace merge_bedmethyl
