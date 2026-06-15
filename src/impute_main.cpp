#include "impute_methylation/beta_binomial.hpp"

#include <cstring>
#include <iostream>
#include <string>

namespace {

void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " [-w BP] [-a ALPHA] [-b BETA] [-n MIN_NEIGHBORS] [-j N] [--sample] [--counts-cov]\n"
        << " <input.tsv> <output.tsv>\n\n"
        << "Local beta-binomial imputation for all samples in the TSV (streaming).\n"
        << "  -w   Genomic window in bp (default 200)\n"
        << "  -a   Beta-binomial prior alpha (default 1)\n"
        << "  -b   Beta-binomial prior beta (default 1)\n"
        << "  -n   Minimum valid neighbors in window (default 5)\n"
        << "  -j   Parallel samples (default 1; 0 = all cores)\n"
        << "  --sample  Input has {id}.counts / {id}.cov columns (default: haplotype)\n"
        << "  --counts-cov  Impute counts and coverage instead of fraction (default: fraction)\n\n"
        << "Haplotype mode (default): writes {id}.hap{1,2}_frac_imputed.\n"
        << "Sample mode (--sample): writes {id}.frac_imputed "
        << "(or {id}.counts_imputed / {id}.cov_imputed with --counts-cov).\n"
        << "Drops source counts/cov/percentage columns (fallback: observed values).\n"
        << "Memory: O(num_targets * window sites), not file size.\n";
}

bool parse_args(int argc, char* argv[], impute_methylation::ImputeOptions& opts,
                std::string& input, std::string& output) {
    int argi = 1;
    while (argi < argc && std::strncmp(argv[argi], "-", 1) == 0) {
        if (std::strcmp(argv[argi], "-w") == 0 && argi + 1 < argc) {
            opts.window_bp = std::stoi(argv[++argi]);
        } else if (std::strcmp(argv[argi], "-a") == 0 && argi + 1 < argc) {
            opts.alpha = std::stod(argv[++argi]);
        } else if (std::strcmp(argv[argi], "-b") == 0 && argi + 1 < argc) {
            opts.beta = std::stod(argv[++argi]);
        } else if (std::strcmp(argv[argi], "-n") == 0 && argi + 1 < argc) {
            opts.min_neighbors = std::stoi(argv[++argi]);
        } else if (std::strcmp(argv[argi], "-j") == 0 && argi + 1 < argc) {
            opts.num_threads = std::stoi(argv[++argi]);
        } else if (std::strcmp(argv[argi], "--sample") == 0) {
            opts.sample_mode = true;
        } else if (std::strcmp(argv[argi], "--counts-cov") == 0) {
            opts.mode = impute_methylation::ImputeMode::CountsCov;
        } else if (std::strcmp(argv[argi], "-h") == 0 ||
                   std::strcmp(argv[argi], "--help") == 0) {
            print_usage(argv[0]);
            return false;
        } else {
            std::cerr << "Unknown option: " << argv[argi] << '\n';
            print_usage(argv[0]);
            return false;
        }
        ++argi;
    }

    if (argc - argi != 2) {
        print_usage(argv[0]);
        return false;
    }

    input = argv[argi];
    output = argv[argi + 1];
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    impute_methylation::ImputeOptions opts;
    std::string input, output;
    if (!parse_args(argc, argv, opts, input, output)) return 1;

    try {
        impute_methylation::stream_beta_binomial_impute_all(input, output, opts);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }

    const char* impute_mode =
        opts.mode == impute_methylation::ImputeMode::CountsCov ? "counts/cov" : "fraction";
    std::cerr << "Wrote " << output << " ("
              << (opts.sample_mode ? "sample" : "haplotype") << " mode, impute=" << impute_mode
              << ", window=" << opts.window_bp << " bp, min_neighbors=" << opts.min_neighbors
              << ", streaming)\n";
    return 0;
}
