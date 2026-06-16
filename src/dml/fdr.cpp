#include "dml/fdr.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

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

bool site_genomic_order(const SiteResult& a, const SiteResult& b) {
    const ChrSortKey ka = chr_sort_key(a.chr);
    const ChrSortKey kb = chr_sort_key(b.chr);
    if (ka.group != kb.group) return ka.group < kb.group;
    if (ka.num != kb.num) return ka.num < kb.num;
    if (ka.tail != kb.tail) return ka.tail < kb.tail;
    if (a.chr != b.chr) return a.chr < b.chr;
    return a.pos < b.pos;
}

}  // namespace

void apply_fdr_bh(std::vector<SiteResult>& results) {
    if (results.empty()) return;

    const std::size_t m = results.size();
    std::vector<std::size_t> order(m);
    for (std::size_t i = 0; i < m; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return results[a].pvalue < results[b].pvalue;
    });

    std::vector<double> adj(m, 1.0);
    double running_min = 1.0;
    for (std::size_t rank = m; rank >= 1; --rank) {
        const std::size_t idx = order[rank - 1];
        const double p = results[idx].pvalue;
        const double q = p * static_cast<double>(m) / static_cast<double>(rank);
        running_min = std::min(running_min, q);
        adj[rank - 1] = running_min;
    }

    for (std::size_t rank = 0; rank < m; ++rank) {
        const std::size_t idx = order[rank];
        results[idx].fdr = std::min(adj[rank], 1.0);
        results[idx].significant = results[idx].fdr < 0.05;
    }
}

void write_results_csv(const std::string& path, std::vector<SiteResult>& results) {
    apply_fdr_bh(results);

    std::sort(results.begin(), results.end(), site_genomic_order);

    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open output: " + path);

    out << "chr,pos,beta_phenotype,se_phenotype,pvalue,phi,mean_case,mean_control,delta_beta,"
           "n_samples,FDR,significant\n";
    out << std::setprecision(10);

    for (const auto& r : results) {
        auto w = [&](double v) {
            if (std::isfinite(v)) out << v;
        };
        out << r.chr << ',' << r.pos << ',';
        w(r.beta_phenotype);
        out << ',';
        w(r.se_phenotype);
        out << ',';
        w(r.pvalue);
        out << ',';
        w(r.phi);
        out << ',';
        w(r.mean_case);
        out << ',';
        w(r.mean_control);
        out << ',';
        w(r.delta_beta);
        out << ',' << r.n_samples << ',';
        w(r.fdr);
        out << ',' << (r.significant ? "True" : "False") << '\n';
    }
}

}  // namespace dml
