#include "dml/model.hpp"

#include "dml/util.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace dml {

namespace {

constexpr int kPhenotypeIdx = 1;
constexpr double kArcsinC0 = 0.1;
constexpr double kPhiMin = 0.001;
constexpr double kPhiMax = 0.999;

struct Obs {
    double k;
    double n;
    std::vector<double> x;
};

double arcsin_methylation(double y, double n) {
    const double frac = (y + kArcsinC0) / (n + 2.0 * kArcsinC0);
    return std::asin(2.0 * frac - 1.0);
}

double predict(const Obs& o, const std::vector<double>& beta) {
    double eta = 0.0;
    for (std::size_t j = 0; j < beta.size(); ++j) {
        eta += beta[j] * o.x[j];
    }
    return eta;
}

bool cholesky_solve(std::vector<double>& a, int p, std::vector<double>& b) {
    std::vector<double> L(static_cast<std::size_t>(p * p), 0.0);
    auto at = [&](int r, int c) -> double& { return a[static_cast<std::size_t>(r * p + c)]; };
    auto lt = [&](int r, int c) -> double& { return L[static_cast<std::size_t>(r * p + c)]; };

    for (int i = 0; i < p; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = at(i, j);
            for (int k = 0; k < j; ++k) sum -= lt(i, k) * lt(j, k);
            if (i == j) {
                if (sum <= 1e-12) return false;
                lt(i, j) = std::sqrt(sum);
            } else {
                lt(i, j) = sum / lt(j, j);
            }
        }
    }

    for (int i = 0; i < p; ++i) {
        double sum = b[static_cast<std::size_t>(i)];
        for (int k = 0; k < i; ++k) sum -= lt(i, k) * b[static_cast<std::size_t>(k)];
        b[static_cast<std::size_t>(i)] = sum / lt(i, i);
    }
    for (int i = p - 1; i >= 0; --i) {
        double sum = b[static_cast<std::size_t>(i)];
        for (int k = i + 1; k < p; ++k) sum -= lt(k, i) * b[static_cast<std::size_t>(k)];
        b[static_cast<std::size_t>(i)] = sum / lt(i, i);
    }
    return true;
}

bool invert_spd(std::vector<double> a, int p, std::vector<double>& inv) {
    inv.assign(static_cast<std::size_t>(p * p), 0.0);
    std::vector<double> rhs(static_cast<std::size_t>(p), 0.0);
    std::vector<double> col(static_cast<std::size_t>(p), 0.0);
    for (int j = 0; j < p; ++j) {
        std::fill(rhs.begin(), rhs.end(), 0.0);
        rhs[static_cast<std::size_t>(j)] = 1.0;
        col = rhs;
        std::vector<double> work = a;
        if (!cholesky_solve(work, p, col)) return false;
        for (int i = 0; i < p; ++i) {
            inv[static_cast<std::size_t>(i * p + j)] = col[static_cast<std::size_t>(i)];
        }
    }
    return true;
}

bool wls(const std::vector<Obs>& obs, int p, const std::vector<double>& z,
         const std::vector<double>& weights, std::vector<double>& beta) {
    std::vector<double> xtx(static_cast<std::size_t>(p * p), 0.0);
    beta.assign(static_cast<std::size_t>(p), 0.0);
    std::vector<double> xtz(static_cast<std::size_t>(p), 0.0);

    for (std::size_t i = 0; i < obs.size(); ++i) {
        const double w = weights[i];
        const double zi = z[i];
        for (int j = 0; j < p; ++j) {
            const double xj = obs[i].x[static_cast<std::size_t>(j)];
            xtz[static_cast<std::size_t>(j)] += w * xj * zi;
            for (int k = 0; k < p; ++k) {
                xtx[static_cast<std::size_t>(j * p + k)] +=
                    w * xj * obs[i].x[static_cast<std::size_t>(k)];
            }
        }
    }
    beta = xtz;
    return cholesky_solve(xtx, p, beta);
}

bool dss_fit_one_cpg(const std::vector<Obs>& obs, int p, std::vector<double>& beta, double& phi,
                     double& se_phenotype) {
    const int n = static_cast<int>(obs.size());
    if (n <= p + 1) return false;
    if (p <= kPhenotypeIdx) return false;

    std::vector<double> z(static_cast<std::size_t>(n));
    std::vector<double> N(static_cast<std::size_t>(n));
    std::vector<double> w1(static_cast<std::size_t>(n));

    for (int i = 0; i < n; ++i) {
        N[static_cast<std::size_t>(i)] = obs[static_cast<std::size_t>(i)].n;
        z[static_cast<std::size_t>(i)] =
            arcsin_methylation(obs[static_cast<std::size_t>(i)].k, N[static_cast<std::size_t>(i)]);
        w1[static_cast<std::size_t>(i)] = N[static_cast<std::size_t>(i)];
    }

    if (!wls(obs, p, z, w1, beta)) return false;

    double sum_sq = 0.0;
    double sum_nm1 = 0.0;
    for (int i = 0; i < n; ++i) {
        const double resid = z[static_cast<std::size_t>(i)] - predict(obs[static_cast<std::size_t>(i)], beta);
        sum_sq += resid * resid * N[static_cast<std::size_t>(i)];
        sum_nm1 += N[static_cast<std::size_t>(i)] - 1.0;
    }
    if (sum_nm1 <= 0.0) return false;

    const double df = static_cast<double>(n - p);
    phi = (sum_sq - df) * static_cast<double>(n) / df / sum_nm1;
    phi = std::clamp(phi, kPhiMin, kPhiMax);

    std::vector<double> w2(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double ni = N[static_cast<std::size_t>(i)];
        w2[static_cast<std::size_t>(i)] = ni / (1.0 + (ni - 1.0) * phi);
    }

    std::vector<double> xtx(static_cast<std::size_t>(p * p), 0.0);
    std::vector<double> xtz(static_cast<std::size_t>(p), 0.0);
    for (std::size_t i = 0; i < obs.size(); ++i) {
        const double w = w2[i];
        const double zi = z[i];
        for (int j = 0; j < p; ++j) {
            const double xj = obs[i].x[static_cast<std::size_t>(j)];
            xtz[static_cast<std::size_t>(j)] += w * xj * zi;
            for (int k = 0; k < p; ++k) {
                xtx[static_cast<std::size_t>(j * p + k)] +=
                    w * xj * obs[i].x[static_cast<std::size_t>(k)];
            }
        }
    }

    beta = xtz;
    if (!cholesky_solve(xtx, p, beta)) return false;

    std::vector<double> xtx_inv;
    if (!invert_spd(xtx, p, xtx_inv)) return false;

    const double var = xtx_inv[static_cast<std::size_t>(kPhenotypeIdx * p + kPhenotypeIdx)];
    if (var <= 1e-12 || !std::isfinite(var)) return false;

    se_phenotype = std::sqrt(var);
    return std::isfinite(beta[static_cast<std::size_t>(kPhenotypeIdx)]) && se_phenotype > 0.0;
}

}  // namespace

std::optional<SiteResult> fit_site(const SiteInput& site, const CohortMeta& meta) {
    const int p = meta.n_params;
    std::vector<Obs> obs;
    obs.reserve(meta.samples.size());

    int case_n = 0;
    int control_n = 0;
    double case_sum = 0.0;
    double control_sum = 0.0;

    for (std::size_t i = 0; i < meta.samples.size(); ++i) {
        const double m = site.counts[i];
        const double c = site.cov[i];
        if (is_nan(m) || is_nan(c) || c <= 0.0) continue;

        obs.push_back(Obs{m, c, meta.samples[i].xrow});
        const double frac = m / c;
        if (meta.samples[i].phenotype_bin == 1) {
            case_sum += frac;
            ++case_n;
        } else {
            control_sum += frac;
            ++control_n;
        }
    }

    if (static_cast<int>(obs.size()) <= p + 1 || case_n == 0 || control_n == 0) {
        return std::nullopt;
    }

    std::vector<double> beta;
    double phi = kNaN;
    double se1 = kNaN;
    if (!dss_fit_one_cpg(obs, p, beta, phi, se1)) return std::nullopt;

    const double stat = beta[static_cast<std::size_t>(kPhenotypeIdx)] / se1;

    SiteResult out;
    out.chr = site.chr;
    out.pos = site.pos;
    out.beta_phenotype = beta[static_cast<std::size_t>(kPhenotypeIdx)];
    out.se_phenotype = se1;
    out.pvalue = 2.0 * normal_sf(std::abs(stat));
    out.phi = phi;
    out.mean_case = case_sum / static_cast<double>(case_n);
    out.mean_control = control_sum / static_cast<double>(control_n);
    out.delta_beta = out.mean_case - out.mean_control;
    out.n_samples = static_cast<int>(obs.size());
    return out;
}

}  // namespace dml
