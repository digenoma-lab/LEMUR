#include "merge_bedmethyl/options.hpp"

#include <cstring>
#include <iostream>

namespace merge_bedmethyl {

void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " [-c N] [-s M] <output.tsv> <label1> <hp1> <hp2> ...\n\n"
        << "  -c, --min-cov N      Minimum coverage to report (column 10; default 3).\n"
        << "                       Values are included only when coverage > N.\n"
        << "  -s, --min-samples M  Minimum samples with data required per row.\n"
        << "                       Default: N-1 (N = number of sample pairs).\n\n"
        << "Writes chr, pos, then per sample: hap1_counts, hap2_counts, hap1_cov,\n"
        << "hap2_cov, hap1_percentage, hap2_percentage ('.' if missing).\n"
        << "Streams all inputs line-by-line; loci are merged on chr + start.\n\n"
        << "Imputation (--impute):\n"
        << "  --impute             After merge, run beta-binomial imputation on output.\n"
        << "                       Default: {id}.hap{1,2}_frac_imputed columns.\n"
        << "  --counts-cov         With --impute: output counts/cov imputed columns instead.\n"
        << "  -w BP                 Genomic window in bp (default 200)\n"
        << "  -a ALPHA              Beta-binomial prior alpha (default 1)\n"
        << "  -b BETA               Beta-binomial prior beta (default 1)\n"
        << "  -n MIN_NEIGHBORS      Minimum valid neighbors in window (default 5)\n"
        << "  -j N                  Parallel imputation samples (default 1; 0 = all cores)\n";
}

int parse_arguments(int argc, char* argv[], ParsedArgs& out) {
    out = ParsedArgs{};
    int argi = 1;

    while (argi < argc && std::strncmp(argv[argi], "-", 1) == 0) {
        if (std::strcmp(argv[argi], "-c") == 0 || std::strcmp(argv[argi], "--min-cov") == 0) {
            if (argi + 1 >= argc) {
                print_usage(argv[0]);
                return 1;
            }
            out.options.min_coverage = std::stoi(argv[++argi]);
            ++argi;
        } else if (std::strcmp(argv[argi], "-s") == 0 ||
                   std::strcmp(argv[argi], "--min-samples") == 0) {
            if (argi + 1 >= argc) {
                print_usage(argv[0]);
                return 1;
            }
            out.options.min_samples = std::stoi(argv[++argi]);
            ++argi;
        } else if (std::strcmp(argv[argi], "--impute") == 0) {
            out.options.impute = true;
            ++argi;
        } else if (std::strcmp(argv[argi], "--counts-cov") == 0) {
            out.options.impute_options.mode = impute_methylation::ImputeMode::CountsCov;
            ++argi;
        } else if (std::strcmp(argv[argi], "-w") == 0) {
            if (argi + 1 >= argc) {
                print_usage(argv[0]);
                return 1;
            }
            out.options.impute_options.window_bp = std::stoi(argv[++argi]);
            ++argi;
        } else if (std::strcmp(argv[argi], "-a") == 0) {
            if (argi + 1 >= argc) {
                print_usage(argv[0]);
                return 1;
            }
            out.options.impute_options.alpha = std::stod(argv[++argi]);
            ++argi;
        } else if (std::strcmp(argv[argi], "-b") == 0) {
            if (argi + 1 >= argc) {
                print_usage(argv[0]);
                return 1;
            }
            out.options.impute_options.beta = std::stod(argv[++argi]);
            ++argi;
        } else if (std::strcmp(argv[argi], "-n") == 0) {
            if (argi + 1 >= argc) {
                print_usage(argv[0]);
                return 1;
            }
            out.options.impute_options.min_neighbors = std::stoi(argv[++argi]);
            ++argi;
        } else if (std::strcmp(argv[argi], "-j") == 0) {
            if (argi + 1 >= argc) {
                print_usage(argv[0]);
                return 1;
            }
            out.options.impute_options.num_threads = std::stoi(argv[++argi]);
            ++argi;
        } else if (std::strcmp(argv[argi], "-h") == 0 ||
                   std::strcmp(argv[argi], "--help") == 0) {
            print_usage(argv[0]);
            out.show_help = true;
            return 0;
        } else {
            std::cerr << "Unknown option: " << argv[argi] << '\n';
            print_usage(argv[0]);
            return 1;
        }
    }

    const int positional = argc - argi;
    if (positional < 4 || (positional - 1) % 3 != 0) {
        print_usage(argv[0]);
        return 1;
    }

    out.output_path = argv[argi];
    out.samples.reserve((positional - 1) / 3);

    for (int i = argi + 1; i < argc; i += 3) {
        out.samples.push_back({argv[i], argv[i + 1], argv[i + 2]});
    }

    return 0;
}

}  // namespace merge_bedmethyl
