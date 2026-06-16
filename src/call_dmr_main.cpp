#include "dml/call_dmr.hpp"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Options {
    std::string dml_csv;
    std::string output_csv;
    double p_threshold = 1e-5;
    int dis_merge = 100;
    int min_cg = 3;
    int min_len = 50;
    double pct_sig = 0.5;
    int num_threads = 1;
    bool sample_mode = false;
};

void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " --sample [-j N] [--p-threshold P] [--dis-merge BP] [--minCG N] [--minlen BP] "
           "[--pct-sig F]\n"
        << "       <dml.csv> <output.csv>\n\n"
        << "DSS callDMR: merge significant DML sites into DMRs (streaming by chromosome).\n"
        << "Input: sorted DML CSV from dml --sample.\n"
        << "  --sample         Sample-mode pipeline (required; haplotype not supported)\n"
        << "  -j N              OpenMP threads per chromosome (default 1; 0 = all cores)\n"
        << "  --p-threshold P   Significant CpG p-value cutoff (default 1e-5)\n"
        << "  --dis-merge BP    Max gap between significant CpGs to merge (default 100)\n"
        << "  --minCG N         Minimum significant CpGs per DMR (default 3)\n"
        << "  --minlen BP       Minimum genomic span of DMR (default 50)\n"
        << "  --pct-sig F       Minimum fraction of significant CpGs in span (default 0.5)\n";
}

bool parse_args(int argc, char* argv[], Options& opts) {
    int argi = 1;
    while (argi < argc && std::strncmp(argv[argi], "-", 1) == 0) {
        if (std::strcmp(argv[argi], "-j") == 0 && argi + 1 < argc) {
            opts.num_threads = std::stoi(argv[++argi]);
        } else if (std::strcmp(argv[argi], "--p-threshold") == 0 && argi + 1 < argc) {
            opts.p_threshold = std::stod(argv[++argi]);
        } else if (std::strcmp(argv[argi], "--dis-merge") == 0 && argi + 1 < argc) {
            opts.dis_merge = std::stoi(argv[++argi]);
        } else if (std::strcmp(argv[argi], "--minCG") == 0 && argi + 1 < argc) {
            opts.min_cg = std::stoi(argv[++argi]);
        } else if (std::strcmp(argv[argi], "--minlen") == 0 && argi + 1 < argc) {
            opts.min_len = std::stoi(argv[++argi]);
        } else if (std::strcmp(argv[argi], "--pct-sig") == 0 && argi + 1 < argc) {
            opts.pct_sig = std::stod(argv[++argi]);
        } else if (std::strcmp(argv[argi], "--sample") == 0) {
            opts.sample_mode = true;
        } else if (std::strcmp(argv[argi], "-h") == 0 || std::strcmp(argv[argi], "--help") == 0) {
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

    opts.dml_csv = argv[argi];
    opts.output_csv = argv[argi + 1];
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    Options opts;
    if (!parse_args(argc, argv, opts)) return 1;

    try {
        if (!opts.sample_mode) {
            throw std::runtime_error(
                "call_dmr requires --sample (input must be DML CSV from dml --sample)");
        }

        dml::CallDmrOptions call_opts;
        call_opts.p_threshold = opts.p_threshold;
        call_opts.dis_merge = opts.dis_merge;
        call_opts.min_cg = opts.min_cg;
        call_opts.min_len = opts.min_len;
        call_opts.pct_sig = opts.pct_sig;
        call_opts.num_threads = opts.num_threads;
        call_opts.sample_mode = opts.sample_mode;

        const auto dmrs = dml::call_dmr_file(opts.dml_csv, call_opts);
        dml::write_dmr_csv(opts.output_csv, dmrs);

        std::cerr << "DMRs found: " << dmrs.size() << '\n';
        std::cerr << "Wrote " << opts.output_csv << '\n';
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }

    return 0;
}
