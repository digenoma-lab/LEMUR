#include "merge_bedmethyl/format.hpp"
#include "merge_bedmethyl/merger.hpp"
#include "merge_bedmethyl/reader.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>

namespace merge_bedmethyl {

int run_merge(const ParsedArgs& args) {
    std::vector<PairReaders> pairs;
    pairs.reserve(args.samples.size());
    for (const auto& sample : args.samples) {
        pairs.emplace_back(sample.label, sample.hp1_path, sample.hp2_path);
    }

    std::ofstream out(args.output_path);
    if (!out) {
        std::cerr << "Cannot open output: " << args.output_path << '\n';
        return 1;
    }

    out << "chr\tpos";
    for (const auto& pair : pairs) out << '\t' << pair.label;
    out << '\n';

    const int num_samples = static_cast<int>(pairs.size());
    int min_samples_with_info =
        args.options.min_samples < 0 ? std::max(0, num_samples - 1) : args.options.min_samples;

    if (min_samples_with_info < 0 || min_samples_with_info > num_samples) {
        std::cerr << "Error: --min-samples must be between 0 and " << num_samples << " (got "
                  << min_samples_with_info << ")\n";
        return 1;
    }

    const int min_coverage = args.options.min_coverage;
    std::size_t rows = 0;

    while (any_valid(pairs)) {
        const Locus target = min_locus(pairs);

        for (auto& pair : pairs) {
            skip_until(pair.hp1, target);
            skip_until(pair.hp2, target);
        }

        const int samples_with_info =
            count_samples_with_information(pairs, target, min_coverage);
        if (samples_with_info < min_samples_with_info) {
            advance_readers_at_locus(pairs, target);
            continue;
        }

        out << target.chr << '\t' << target.pos;
        for (auto& pair : pairs) {
            out << '\t' << format_pair_cell(pair.hp1, pair.hp2, target, min_coverage);
        }
        out << '\n';
        ++rows;

        advance_readers_at_locus(pairs, target);
    }

    std::cerr << "Wrote " << rows << " rows to " << args.output_path
              << " (min coverage > " << min_coverage << ", >= " << min_samples_with_info << "/"
              << num_samples << " samples per row)\n";
    return 0;
}

}  // namespace merge_bedmethyl
