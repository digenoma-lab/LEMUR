#include "impute_methylation/beta_binomial.hpp"
#include "merge_bedmethyl/format.hpp"
#include "merge_bedmethyl/merger.hpp"
#include "merge_bedmethyl/reader.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>

namespace merge_bedmethyl {

namespace {

int run_hap_merge(const ParsedArgs& args, const std::string& merge_output_path) {
    std::vector<PairReaders> pairs;
    pairs.reserve(args.samples.size());
    for (const auto& sample : args.samples) {
        pairs.emplace_back(sample.label, sample.path1, sample.path2);
    }

    std::ofstream out(merge_output_path);
    if (!out) {
        std::cerr << "Cannot open output: " << merge_output_path << '\n';
        return 1;
    }

    out << "chr\tpos";
    for (const auto& pair : pairs) {
        const std::string& id = pair.label;
        out << '\t' << id << ".hap1_counts" << '\t' << id << ".hap2_counts" << '\t' << id
            << ".hap1_cov" << '\t' << id << ".hap2_cov" << '\t' << id << ".hap1_percentage"
            << '\t' << id << ".hap2_percentage";
    }
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
            append_haplotype_columns(out, pair.hp1, pair.hp2, target, min_coverage);
        }
        out << '\n';
        ++rows;

        advance_readers_at_locus(pairs, target);
    }

    out.close();

    std::cerr << "Wrote " << rows << " rows to " << merge_output_path
              << " (min coverage > " << min_coverage << ", >= " << min_samples_with_info << "/"
              << num_samples << " samples per row)\n";
    return 0;
}

int run_sample_merge(const ParsedArgs& args, const std::string& merge_output_path) {
    std::vector<SingleReader> samples;
    samples.reserve(args.samples.size());
    for (const auto& sample : args.samples) {
        samples.emplace_back(sample.label, sample.path1);
    }

    std::ofstream out(merge_output_path);
    if (!out) {
        std::cerr << "Cannot open output: " << merge_output_path << '\n';
        return 1;
    }

    out << "chr\tpos";
    for (const auto& sample : samples) {
        const std::string& id = sample.label;
        out << '\t' << id << ".counts" << '\t' << id << ".cov" << '\t' << id << ".percentage";
    }
    out << '\n';

    const int num_samples = static_cast<int>(samples.size());
    int min_samples_with_info =
        args.options.min_samples < 0 ? std::max(0, num_samples - 1) : args.options.min_samples;

    if (min_samples_with_info < 0 || min_samples_with_info > num_samples) {
        std::cerr << "Error: --min-samples must be between 0 and " << num_samples << " (got "
                  << min_samples_with_info << ")\n";
        return 1;
    }

    const int min_coverage = args.options.min_coverage;
    std::size_t rows = 0;

    while (any_valid(samples)) {
        const Locus target = min_locus(samples);

        for (auto& sample : samples) {
            skip_until(sample.reader, target);
        }

        const int samples_with_info =
            count_samples_with_information(samples, target, min_coverage);
        if (samples_with_info < min_samples_with_info) {
            advance_readers_at_locus(samples, target);
            continue;
        }

        out << target.chr << '\t' << target.pos;
        for (auto& sample : samples) {
            append_sample_columns(out, sample.reader, target, min_coverage);
        }
        out << '\n';
        ++rows;

        advance_readers_at_locus(samples, target);
    }

    out.close();

    std::cerr << "Wrote " << rows << " rows to " << merge_output_path
              << " (min coverage > " << min_coverage << ", >= " << min_samples_with_info << "/"
              << num_samples << " samples per row)\n";
    return 0;
}

}  // namespace

int run_merge(const ParsedArgs& args) {
    const std::string merge_output_path =
        args.options.impute ? args.output_path + ".merge.tmp" : args.output_path;

    const int merge_rc =
        args.options.hap_mode ? run_hap_merge(args, merge_output_path)
                              : run_sample_merge(args, merge_output_path);
    if (merge_rc != 0) return merge_rc;

    if (args.options.impute) {
        if (!args.options.hap_mode) {
            std::cerr << "Error: --impute requires --hap (haplotype input)\n";
            std::remove(merge_output_path.c_str());
            return 1;
        }
        try {
            impute_methylation::stream_beta_binomial_impute_all(
                merge_output_path, args.output_path, args.options.impute_options);
        } catch (const std::exception& e) {
            std::cerr << e.what() << '\n';
            std::remove(merge_output_path.c_str());
            return 1;
        }
        std::remove(merge_output_path.c_str());
        std::cerr << "Imputed to " << args.output_path << " (window="
                  << args.options.impute_options.window_bp
                  << " bp, min_neighbors=" << args.options.impute_options.min_neighbors
                  << ")\n";
    }

    return 0;
}

}  // namespace merge_bedmethyl
