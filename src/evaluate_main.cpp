#include "impute_methylation/evaluate.hpp"

#include <cstring>
#include <iostream>
#include <string>

namespace {

void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " -c COUNTS_COL [-chr CHR] [-m FRAC] [-s SEED] [-w BP] [-a A] [-b B] [-n N] <input.tsv>\n\n"
        << "Hold-out evaluation on one haplotype column (e.g. CHI08A.hap1_counts)\n"
        << "and one chromosome (default chr1): mask, impute, score vs ground truth.\n\n"
        << "  -c   Counts column to evaluate (required)\n"
        << "  -chr Chromosome to use (default chr1)\n"
        << "  -m   Mask fraction (default 0.2)\n"
        << "  -s   RNG seed for reproducible mask (default 42)\n"
        << "  -w,-a,-b,-n   Same as impute_methylation (default 200bp, 5 neighbors)\n\n"
        << "Writes <input>.masked.tsv and <input>.imputed.eval.tsv\n"
        << "Metrics: MSE and Pearson on masked sites for that haplotype only.\n";
}

bool parse_args(int argc, char* argv[], impute_methylation::EvaluateOptions& opts,
                std::string& input) {
    int argi = 1;
    while (argi < argc && std::strncmp(argv[argi], "-", 1) == 0) {
        if (std::strcmp(argv[argi], "-c") == 0 && argi + 1 < argc) {
            opts.y_col = argv[++argi];
        } else if (std::strcmp(argv[argi], "-chr") == 0 && argi + 1 < argc) {
            opts.chromosome = argv[++argi];
        } else if (std::strcmp(argv[argi], "-m") == 0 && argi + 1 < argc) {
            opts.mask_fraction = std::stod(argv[++argi]);
        } else if (std::strcmp(argv[argi], "-s") == 0 && argi + 1 < argc) {
            opts.seed = static_cast<uint64_t>(std::stoull(argv[++argi]));
        } else if (std::strcmp(argv[argi], "-w") == 0 && argi + 1 < argc) {
            opts.impute.window_bp = std::stoi(argv[++argi]);
        } else if (std::strcmp(argv[argi], "-a") == 0 && argi + 1 < argc) {
            opts.impute.alpha = std::stod(argv[++argi]);
        } else if (std::strcmp(argv[argi], "-b") == 0 && argi + 1 < argc) {
            opts.impute.beta = std::stod(argv[++argi]);
        } else if (std::strcmp(argv[argi], "-n") == 0 && argi + 1 < argc) {
            opts.impute.min_neighbors = std::stoi(argv[++argi]);
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

    if (argc - argi != 1) {
        print_usage(argv[0]);
        return false;
    }

    input = argv[argi];
    if (opts.y_col.empty()) {
        std::cerr << "Missing required -c <counts_column>\n";
        return false;
    }
    if (opts.mask_fraction <= 0.0 || opts.mask_fraction >= 1.0) {
        std::cerr << "Mask fraction must be in (0, 1)\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    impute_methylation::EvaluateOptions opts;
    std::string input;
    if (!parse_args(argc, argv, opts, input)) return 1;

    try {
        impute_methylation::run_evaluate(input, opts);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
