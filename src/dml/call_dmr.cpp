#include "dml/call_dmr.hpp"

#include "dml/util.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace dml {

namespace {

struct ChrSortKey {
    int group = 2;
    int num = 0;
    std::string tail;
};

ChrSortKey chr_sort_key(const std::string& chr) {
    ChrSortKey key;
    std::string s = chr;
    if (s.size() >= 3 && s.compare(0, 3, "chr") == 0) {
        s = s.substr(3);
    }
    if (!s.empty() && std::all_of(s.begin(), s.end(), ::isdigit)) {
        key.group = 0;
        key.num = std::stoi(s);
        return key;
    }
    std::string upper = s;
    std::transform(upper.begin(), upper.end(), upper.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    if (upper == "X") {
        key.group = 1;
        key.num = 23;
    } else if (upper == "Y") {
        key.group = 1;
        key.num = 24;
    } else if (upper == "M" || upper == "MT") {
        key.group = 1;
        key.num = 25;
    } else {
        key.tail = s;
    }
    return key;
}

bool dmr_genomic_order(const DmrRecord& a, const DmrRecord& b) {
    const ChrSortKey ka = chr_sort_key(a.chr);
    const ChrSortKey kb = chr_sort_key(b.chr);
    if (ka.group != kb.group) return ka.group < kb.group;
    if (ka.num != kb.num) return ka.num < kb.num;
    if (ka.tail != kb.tail) return ka.tail < kb.tail;
    if (a.chr != b.chr) return a.chr < b.chr;
    return a.start < b.start;
}

struct ChrSpan {
    std::string chr;
    std::streamoff start = 0;
};

struct DmlSite {
    int pos = 0;
    double pvalue = kNaN;
    double stat = kNaN;
    double delta_beta = kNaN;
};

struct SigSite {
    int pos = 0;
    double pvalue = 0.0;
    double stat = 0.0;
    double delta_beta = 0.0;
    int direction = 0;
};

struct ColumnMap {
    int chr = -1;
    int pos = -1;
    int beta = -1;
    int se = -1;
    int pvalue = -1;
    int delta_beta = -1;
};

ColumnMap parse_header(const std::string& header_line) {
    const auto fields = split_csv_line(header_line, ',');
    ColumnMap m;
    m.chr = column_index(fields, "chr");
    m.pos = column_index(fields, "pos");
    m.beta = column_index(fields, "beta_phenotype");
    m.se = column_index(fields, "se_phenotype");
    m.pvalue = column_index(fields, "pvalue");
    m.delta_beta = column_index(fields, "delta_beta");
    if (m.chr < 0 || m.pos < 0 || m.beta < 0 || m.se < 0 || m.pvalue < 0 || m.delta_beta < 0) {
        throw std::runtime_error(
            "DML CSV must contain chr,pos,beta_phenotype,se_phenotype,pvalue,delta_beta");
    }
    return m;
}

std::optional<DmlSite> parse_site_line(const std::string& line, const ColumnMap& cols) {
    const auto fields = split_csv_line(line, ',');
    const std::size_t need = static_cast<std::size_t>(
        std::max({cols.chr, cols.pos, cols.beta, cols.se, cols.pvalue, cols.delta_beta}));
    if (fields.size() <= need) return std::nullopt;

    const double beta = parse_double_field(fields[static_cast<std::size_t>(cols.beta)]);
    const double se = parse_double_field(fields[static_cast<std::size_t>(cols.se)]);
    const double pvalue = parse_double_field(fields[static_cast<std::size_t>(cols.pvalue)]);
    const double delta = parse_double_field(fields[static_cast<std::size_t>(cols.delta_beta)]);

    if (is_nan(pvalue)) return std::nullopt;
    if (is_nan(se) || se <= 0.0) return std::nullopt;
    if (is_nan(beta)) return std::nullopt;

    const double stat = beta / se;
    if (!std::isfinite(stat)) return std::nullopt;

    DmlSite site;
    try {
        site.pos = std::stoi(fields[static_cast<std::size_t>(cols.pos)]);
    } catch (...) {
        return std::nullopt;
    }
    site.pvalue = pvalue;
    site.stat = stat;
    site.delta_beta = is_nan(delta) ? 0.0 : delta;
    return site;
}

std::vector<ChrSpan> index_chromosome_spans(const std::string& path, const ColumnMap& cols) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open DML CSV: " + path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty DML CSV: " + path);

    std::vector<ChrSpan> spans;
    std::string cur_chr;
    std::streamoff block_start = in.tellg();
    std::streamoff pos = block_start;

    while (std::getline(in, line)) {
        if (line.empty()) {
            pos = in.tellg();
            continue;
        }
        const auto fields = split_csv_line(line, ',');
        if (static_cast<int>(fields.size()) <= cols.chr) {
            pos = in.tellg();
            continue;
        }
        const std::string& chr = fields[static_cast<std::size_t>(cols.chr)];

        if (chr != cur_chr) {
            if (!cur_chr.empty()) {
                spans.push_back(ChrSpan{cur_chr, block_start});
            }
            cur_chr = chr;
            block_start = pos;
        }
        pos = in.tellg();
    }

    if (!cur_chr.empty()) {
        spans.push_back(ChrSpan{cur_chr, block_start});
    }
    return spans;
}

std::vector<DmlSite> read_chromosome_sites(const std::string& path, const ChrSpan& span,
                                            const ColumnMap& cols) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open DML CSV: " + path);
    in.seekg(span.start);

    std::vector<DmlSite> sites;
    sites.reserve(65536);

    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto fields = split_csv_line(line, ',');
        if (static_cast<int>(fields.size()) <= cols.chr) continue;

        const std::string& chr = fields[static_cast<std::size_t>(cols.chr)];
        if (!first && chr != span.chr) break;
        first = false;

        if (auto site = parse_site_line(line, cols)) {
            sites.push_back(*site);
        }
    }
    return sites;
}

int count_sites_in_range(const std::vector<DmlSite>& sites, int start, int end) {
    const auto lo = std::lower_bound(
        sites.begin(), sites.end(), start,
        [](const DmlSite& s, int p) { return s.pos < p; });
    const auto hi = std::upper_bound(
        sites.begin(), sites.end(), end,
        [](int p, const DmlSite& s) { return p < s.pos; });
    return static_cast<int>(hi - lo);
}

void append_clusters(const std::vector<SigSite>& dir_sig, const std::vector<DmlSite>& sites,
                     const std::string& chr, const CallDmrOptions& opts,
                     std::vector<DmrRecord>& out) {
    if (dir_sig.empty()) return;

    std::size_t i = 0;
    while (i < dir_sig.size()) {
        std::vector<const SigSite*> cluster;
        cluster.push_back(&dir_sig[i]);
        ++i;
        while (i < dir_sig.size() && dir_sig[i].pos - cluster.back()->pos <= opts.dis_merge) {
            cluster.push_back(&dir_sig[i]);
            ++i;
        }

        int start = cluster.front()->pos;
        int end = start;
        for (const SigSite* r : cluster) {
            end = std::max(end, r->pos);
        }

        const int region_len = end - start + 1;
        const int n_sig = static_cast<int>(cluster.size());
        const int n_total = count_sites_in_range(sites, start, end);
        if (n_total <= 0) continue;

        const double prop_sig = static_cast<double>(n_sig) / static_cast<double>(n_total);
        if (n_sig < opts.min_cg) continue;
        if (region_len < opts.min_len) continue;
        if (prop_sig < opts.pct_sig) continue;

        double area_stat = 0.0;
        double sum_stat = 0.0;
        double sum_diff = 0.0;
        double sum_p = 0.0;
        double min_p = cluster.front()->pvalue;
        for (const SigSite* r : cluster) {
            area_stat += r->stat;
            sum_stat += r->stat;
            sum_diff += r->delta_beta;
            sum_p += r->pvalue;
            min_p = std::min(min_p, r->pvalue);
        }
        const double mean_stat = sum_stat / static_cast<double>(n_sig);

        DmrRecord dmr;
        dmr.chr = chr;
        dmr.start = start;
        dmr.end = end;
        dmr.length = region_len;
        dmr.n_cg = n_total;
        dmr.n_cg_sig = n_sig;
        dmr.pct_sig = prop_sig;
        dmr.area_stat = area_stat;
        dmr.mean_stat = mean_stat;
        dmr.mean_diff = sum_diff / static_cast<double>(n_sig);
        dmr.direction = mean_stat > 0.0 ? "hyper" : "hypo";
        dmr.min_p = min_p;
        dmr.mean_p = sum_p / static_cast<double>(n_sig);
        out.push_back(std::move(dmr));
    }
}

std::vector<DmrRecord> process_chromosome(const std::string& chr,
                                          const std::vector<DmlSite>& sites,
                                          const CallDmrOptions& opts) {
    std::vector<DmrRecord> out;
    if (sites.empty()) return out;

    std::vector<SigSite> hyper;
    std::vector<SigSite> hypo;
    hyper.reserve(sites.size() / 100 + 4);
    hypo.reserve(sites.size() / 100 + 4);

    for (const auto& s : sites) {
        if (s.pvalue >= opts.p_threshold) continue;
        SigSite row{s.pos, s.pvalue, s.stat, s.delta_beta, s.stat >= 0.0 ? 1 : -1};
        if (row.direction > 0) {
            hyper.push_back(row);
        } else {
            hypo.push_back(row);
        }
    }

    append_clusters(hyper, sites, chr, opts, out);
    append_clusters(hypo, sites, chr, opts, out);
    return out;
}

}  // namespace

std::vector<DmrRecord> call_dmr_file(const std::string& dml_csv, const CallDmrOptions& opts) {
#ifdef _OPENMP
    if (opts.num_threads == 0) {
        omp_set_num_threads(omp_get_max_threads());
    } else if (opts.num_threads > 0) {
        omp_set_num_threads(opts.num_threads);
    }
#endif

    std::ifstream header_in(dml_csv);
    if (!header_in) throw std::runtime_error("Cannot open DML CSV: " + dml_csv);
    std::string header_line;
    if (!std::getline(header_in, header_line)) {
        throw std::runtime_error("Empty DML CSV: " + dml_csv);
    }
    const ColumnMap cols = parse_header(header_line);
    const std::vector<ChrSpan> spans = index_chromosome_spans(dml_csv, cols);

    std::vector<std::vector<DmrRecord>> per_chr(spans.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (std::size_t i = 0; i < spans.size(); ++i) {
        const std::vector<DmlSite> sites = read_chromosome_sites(dml_csv, spans[i], cols);
        per_chr[i] = process_chromosome(spans[i].chr, sites, opts);
    }

    std::size_t total = 0;
    for (const auto& v : per_chr) total += v.size();

    std::vector<DmrRecord> dmrs;
    dmrs.reserve(total);
    for (auto& v : per_chr) {
        dmrs.insert(dmrs.end(), std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
    }

    std::sort(dmrs.begin(), dmrs.end(), dmr_genomic_order);
    return dmrs;
}

void write_dmr_csv(const std::string& path, const std::vector<DmrRecord>& dmrs) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open output: " + path);

    out << "chr,start,end,length,nCG,nCG.sig,pct.sig,areaStat,meanStat,meanDiff,direction,minP,"
           "meanP\n";
    out << std::setprecision(10);

    for (const auto& r : dmrs) {
        out << r.chr << ',' << r.start << ',' << r.end << ',' << r.length << ',' << r.n_cg << ','
            << r.n_cg_sig << ',' << r.pct_sig << ',' << r.area_stat << ',' << r.mean_stat << ','
            << r.mean_diff << ',' << r.direction << ',' << r.min_p << ',' << r.mean_p << '\n';
    }
}

}  // namespace dml
