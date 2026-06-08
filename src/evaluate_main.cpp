#include "impute_methylation/evaluate.hpp"

#include <cstring>
#include <iostream>
#include <string>

namespace {

void print_usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog
        << " -o OUT.tsv [-chr CHR] [-m FRAC] [-s SEED] [-w BP] [-a A] [-b B] [-n N] [-j N] "
           "<input.tsv>\n"
        << "      Evaluate all {sample}.hap{1,2} columns in parallel.\n\n"
        << "  " << prog
        << " -c COUNTS_COL [-chr CHR] [-m FRAC] [-s SEED] [-w BP] [-a A] [-b B] [-n N] "
           "<input.tsv>\n"
        << "      Hold-out evaluation on one haplotype column.\n\n"
        << "  -o   Cohort summary TSV (required without -c)\n"
        << "  -c   Counts column for single-target mode (e.g. CHI08A.hap1_counts)\n"
        << "  -chr Chromosome to use (default chr1)\n"
        << "  -m   Mask fraction (default 0.2)\n"
        << "  -s   RNG seed for reproducible mask (default 42)\n"
        << "  -w,-a,-b,-n   Imputation parameters (default 200bp, 1, 1, 5 neighbors)\n"
        << "  -j   Parallel evaluation threads for cohort mode (default 1; 0 = all cores)\n\n"
        << "Cohort output columns:\n"
        << "  sample,window_size,n_neighbors,a,b,chr,pearson,mse,count_masked,count_imputed\n\n"
        << "Single-target mode writes <input>.masked.tsv and <input>.imputed.eval.tsv\n";
}

struct ParsedArgs {
    impute_methylation::EvaluateOptions single;
    impute_methylation::EvaluateCohortOptions cohort;
    std::string input;
    bool cohort_mode = false;
};

bool parse_args(int argc, char* argv[], ParsedArgs& args) {
    int argi = 1;
    while (argi < argc && std::strncmp(argv[argi], "-", 1) == 0) {
        if (std::strcmp(argv[argi], "-c") == 0 && argi + 1 < argc) {
            args.single.y_col = argv[++argi];
        } else if (std::strcmp(argv[argi], "-o") == 0 && argi + 1 < argc) {
            args.cohort.output_path = argv[++argi];
        } else if (std::strcmp(argv[argi], "-chr") == 0 && argi + 1 < argc) {
            const std::string chr = argv[++argi];
            args.single.chromosome = chr;
            args.cohort.chromosome = chr;
        } else if (std::strcmp(argv[argi], "-m") == 0 && argi + 1 < argc) {
            const double frac = std::stod(argv[++argi]);
            args.single.mask_fraction = frac;
            args.cohort.mask_fraction = frac;
        } else if (std::strcmp(argv[argi], "-s") == 0 && argi + 1 < argc) {
            const uint64_t seed = static_cast<uint64_t>(std::stoull(argv[++argi]));
            args.single.seed = seed;
            args.cohort.seed = seed;
        } else if (std::strcmp(argv[argi], "-w") == 0 && argi + 1 < argc) {
            const int window = std::stoi(argv[++argi]);
            args.single.impute.window_bp = window;
            args.cohort.impute.window_bp = window;
        } else if (std::strcmp(argv[argi], "-a") == 0 && argi + 1 < argc) {
            const double alpha = std::stod(argv[++argi]);
            args.single.impute.alpha = alpha;
            args.cohort.impute.alpha = alpha;
        } else if (std::strcmp(argv[argi], "-b") == 0 && argi + 1 < argc) {
            const double beta = std::stod(argv[++argi]);
            args.single.impute.beta = beta;
            args.cohort.impute.beta = beta;
        } else if (std::strcmp(argv[argi], "-n") == 0 && argi + 1 < argc) {
            const int neighbors = std::stoi(argv[++argi]);
            args.single.impute.min_neighbors = neighbors;
            args.cohort.impute.min_neighbors = neighbors;
        } else if (std::strcmp(argv[argi], "-j") == 0 && argi + 1 < argc) {
            args.cohort.impute.num_threads = std::stoi(argv[++argi]);
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

    args.input = argv[argi];
    args.cohort_mode = args.single.y_col.empty();

    if (args.cohort_mode) {
        if (args.cohort.output_path.empty()) {
            std::cerr << "Cohort mode requires -o <output.tsv>\n";
            return false;
        }
    }

    const double mask_fraction =
        args.cohort_mode ? args.cohort.mask_fraction : args.single.mask_fraction;
    if (mask_fraction <= 0.0 || mask_fraction >= 1.0) {
        std::cerr << "Mask fraction must be in (0, 1)\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    ParsedArgs args;
    if (!parse_args(argc, argv, args)) return 1;

    try {
        if (args.cohort_mode) {
            impute_methylation::run_evaluate_cohort(args.input, args.cohort);
        } else {
            impute_methylation::run_evaluate(args.input, args.single);
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
