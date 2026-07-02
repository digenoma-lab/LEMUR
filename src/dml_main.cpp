#include "dml/fdr.hpp"
#include "dml/metadata.hpp"
#include "dml/model.hpp"
#include "dml/smoothing.hpp"
#include "dml/tsv_stream.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Options {
    std::string methylation_tsv;
    std::string metadata_csv;
    std::string output_csv;
    std::string case_label = "Case";
    std::string control_label = "Control";
    bool sample_mode = false;
    bool smoothing = false;
    int smoothing_span = 500;
    int num_threads = 1;
    std::size_t batch_size = 16384;
};

void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " [-j N] [-b BATCH] [--sample] [--smoothing] [--smoothing-span BP]\n"
        << "       [--case-label L] [--control-label L]\n"
        << "       <methylation.tsv> <metadata.csv> <output.csv>\n\n"
        << "DSS DML multiFactor: per-CpG differential methylation (arcsin + WLS).\n"
        << "  -j N                 OpenMP threads (default 1; 0 = all cores)\n"
        << "  -b BATCH             Sites per read/fit batch (default 16384)\n"
        << "  --sample             Input has {id}.counts / {id}.cov (required)\n"
        << "  --smoothing          Apply DSS moving-average smoothing before fitting\n"
        << "  --smoothing-span BP  Smoothing window size in bp (default 500)\n"
        << "  --case-label         Case phenotype label (default Case)\n"
        << "  --control-label      Control phenotype label (default Control)\n";
}

bool parse_args(int argc, char* argv[], Options& opts) {
    int argi = 1;
    while (argi < argc && std::strncmp(argv[argi], "-", 1) == 0) {
        if (std::strcmp(argv[argi], "-j") == 0 && argi + 1 < argc) {
            opts.num_threads = std::stoi(argv[++argi]);
        } else if (std::strcmp(argv[argi], "-b") == 0 && argi + 1 < argc) {
            opts.batch_size = static_cast<std::size_t>(std::stoul(argv[++argi]));
        } else if (std::strcmp(argv[argi], "--case-label") == 0 && argi + 1 < argc) {
            opts.case_label = argv[++argi];
        } else if (std::strcmp(argv[argi], "--control-label") == 0 && argi + 1 < argc) {
            opts.control_label = argv[++argi];
        } else if (std::strcmp(argv[argi], "--sample") == 0) {
            opts.sample_mode = true;
        } else if (std::strcmp(argv[argi], "--smoothing") == 0) {
            opts.smoothing = true;
        } else if (std::strcmp(argv[argi], "--smoothing-span") == 0 && argi + 1 < argc) {
            opts.smoothing_span = std::stoi(argv[++argi]);
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

    if (argc - argi != 3) {
        print_usage(argv[0]);
        return false;
    }

    opts.methylation_tsv = argv[argi];
    opts.metadata_csv = argv[argi + 1];
    opts.output_csv = argv[argi + 2];
    return true;
}

void fit_sites(const std::vector<dml::SiteInput>& sites, const dml::CohortMeta& meta,
               std::vector<dml::SiteResult>& results) {
    std::vector<std::optional<dml::SiteResult>> batch_out(sites.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < sites.size(); ++i) {
        batch_out[i] = dml::fit_site(sites[i], meta);
    }

    for (const auto& fit : batch_out) {
        if (fit) results.push_back(*fit);
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    Options opts;
    if (!parse_args(argc, argv, opts)) return 1;

#ifdef _OPENMP
    if (opts.num_threads == 0) {
        opts.num_threads = omp_get_max_threads();
    }
    omp_set_num_threads(opts.num_threads);
#endif

    try {
        if (!opts.sample_mode) {
            throw std::runtime_error(
                "dml requires --sample (input must have {id}.counts / {id}.cov "
                "columns from impute_methylation --sample --counts-cov)");
        }
        if (opts.smoothing && opts.smoothing_span <= 0) {
            throw std::runtime_error("--smoothing-span must be positive");
        }

        const std::vector<std::string> header = dml::read_tsv_header(opts.methylation_tsv);
        const dml::CohortMeta meta = dml::load_metadata(opts.metadata_csv, header, opts.sample_mode,
                                                      opts.case_label, opts.control_label);

        std::vector<dml::SiteResult> results;
        results.reserve(1'000'000);
        std::size_t rows_read = 0;

        if (opts.smoothing) {
            std::vector<dml::SiteInput> sites =
                dml::load_all_sites(opts.methylation_tsv, meta);
            rows_read = sites.size();
            if (sites.empty()) throw std::runtime_error("No sites in methylation TSV");

            dml::apply_dss_smoothing(sites, opts.smoothing_span);
            fit_sites(sites, meta, results);
        } else {
            dml::MethylationStream stream(opts.methylation_tsv, meta);
            std::vector<dml::SiteInput> batch;

            while (stream.read_batch(batch, opts.batch_size)) {
                rows_read += batch.size();
                fit_sites(batch, meta, results);

                if (rows_read % (opts.batch_size * 10) == 0) {
                    std::cerr << "Processed " << rows_read << " sites, kept " << results.size()
                              << " tests...\n";
                }
            }
        }

        if (results.empty()) {
            throw std::runtime_error("No sites passed model fitting");
        }

        dml::write_results_csv(opts.output_csv, results);

        std::size_t sig = 0;
        for (const auto& r : results) {
            if (r.significant) ++sig;
        }

        std::cerr << "Read " << rows_read << " sites, tested " << results.size() << " CpGs";
        if (opts.smoothing) {
            std::cerr << " (smoothing span=" << opts.smoothing_span << " bp)";
        }
        std::cerr << '\n';
        std::cerr << "Significant DMCs (FDR<0.05): " << sig << '\n';
        std::cerr << "Wrote " << opts.output_csv << '\n';
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }

    return 0;
}
