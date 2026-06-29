#pragma once
// ============================================================================
//  StatTest.hpp -- statistical goodness-of-fit utilities for the ensemble tests.
//
//  Header-only, no Boost: the chi-square critical value is obtained from the
//  Wilson-Hilferty cube-root normal approximation combined with Acklam's inverse
//  normal CDF, accurate to a few parts in 10^4 in the upper tail -- far tighter
//  than any sampling test needs.
//
//  TOLERANCE PHILOSOPHY (Rule 8 -- a test encodes WHY behavior matters):
//   * For MEANS, assert |observed - expected| <= k * stderr with k = 4
//     (two-sided false-fail ~6e-5), and report the z-score on failure.
//   * For DISTRIBUTIONS, use chi-square / G goodness-of-fit at alpha = 1e-4, so
//     the whole suite's expected false-fail rate stays well under 1% across the
//     ~15 statistical assertions. The sample size N for each test is chosen so the
//     test has >0.99 power against the SPECIFIC wrong hypothesis it must reject
//     (e.g. Haar vs Euler-flat) -- documented at each call site.
//
//  Both chiSquareStatistic and gStatistic take OBSERVED integer counts and
//  EXPECTED real counts (same total). The G-test (likelihood ratio) is preferred
//  when some bins are lightly populated; chi-square is the classic Pearson form.
// ============================================================================

#include <cmath>
#include <vector>

namespace rtest {
namespace stat {

// ---- inverse normal CDF (Acklam), upper-tail quantile z s.t. P(Z>z)=p --------
inline double invNormalCdf(double p) {
    // lower-tail quantile: returns z with P(Z<=z) = p, 0<p<1.
    static const double a[] = {-3.969683028665376e+01,
                              2.209460984245205e+02,
                              -2.759285104469687e+02,
                              1.383577518672690e+02,
                              -3.066479806614716e+01,
                              2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01,
                              1.615858368580409e+02,
                              -1.556989798598866e+02,
                              6.680131188771972e+01,
                              -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03,
                              -3.223964580411365e-01,
                              -2.400758277161838e+00,
                              -2.549732539343734e+00,
                              4.374664141464968e+00,
                              2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03,
                              3.224671290700398e-01,
                              2.445134137142996e+00,
                              3.754408661907416e+00};
    const double plow = 0.02425, phigh = 1.0 - 0.02425;
    if (p < plow) {
        const double q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
               / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    if (p <= phigh) {
        const double q = p - 0.5, r = q * q;
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
               / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    }
    const double q = std::sqrt(-2.0 * std::log(1.0 - p));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
           / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
}

// Upper-tail normal quantile z with P(Z>z) = alpha.
inline double normalUpperQuantile(double alpha) {
    return invNormalCdf(1.0 - alpha);
}

// ---- chi-square upper critical value (Wilson-Hilferty) -----------------------
// chi^2_{1-alpha}(dof): the value the statistic exceeds with probability alpha
// under H0. dof >= 1.
inline double chiSquareCritical(int dof, double alpha) {
    const double z = normalUpperQuantile(alpha);
    const double k = static_cast<double>(dof);
    const double t = 1.0 - 2.0 / (9.0 * k) + z * std::sqrt(2.0 / (9.0 * k));
    return k * t * t * t;
}

// ---- goodness-of-fit statistics ----------------------------------------------
// Pearson chi-square: sum (O-E)^2 / E over bins with E>0.
inline double chiSquareStatistic(const std::vector<long>& observed, const std::vector<double>& expected) {
    double chi2 = 0.0;
    const std::size_t n = observed.size();
    for (std::size_t i = 0; i < n; ++i) {
        const double e = expected[i];
        if (e <= 0.0) {
            continue;
        }
        const double d = static_cast<double>(observed[i]) - e;
        chi2 += d * d / e;
    }
    return chi2;
}

// G-test (likelihood ratio): 2 sum O ln(O/E). Better-behaved for small E.
inline double gStatistic(const std::vector<long>& observed, const std::vector<double>& expected) {
    double g = 0.0;
    const std::size_t n = observed.size();
    for (std::size_t i = 0; i < n; ++i) {
        const long o = observed[i];
        const double e = expected[i];
        if (o > 0 && e > 0.0) {
            g += static_cast<double>(o) * std::log(static_cast<double>(o) / e);
        }
    }
    return 2.0 * g;
}

// ---- binning helpers ---------------------------------------------------------
// Build expected COUNTS for a uniform distribution on [lo,hi] over nbins, total N.
inline std::vector<double> uniformExpected(int nbins, long total) {
    return std::vector<double>(static_cast<std::size_t>(nbins),
                               static_cast<double>(total) / static_cast<double>(nbins));
}

// Build expected COUNTS proportional to a per-bin density weight, normalized to N.
inline std::vector<double> expectedFromWeights(const std::vector<double>& weights, long total) {
    double sum = 0.0;
    for (double w : weights) {
        sum += w;
    }
    std::vector<double> e(weights.size(), 0.0);
    if (sum <= 0.0) {
        return e;
    }
    for (std::size_t i = 0; i < weights.size(); ++i) {
        e[i] = static_cast<double>(total) * weights[i] / sum;
    }
    return e;
}

// A fixed-width histogram on [lo,hi]. Out-of-range samples are clamped into the
// edge bins (callers ensure the support is covered, so this is a guard only).
class Histogram {
    public:
    Histogram(double lo, double hi, int nbins)
        : lo_(lo)
        , hi_(hi)
        , nbins_(nbins)
        , counts_(static_cast<std::size_t>(nbins), 0) {
    }
    void add(double x) {
        int b = static_cast<int>((x - lo_) / (hi_ - lo_) * nbins_);
        if (b < 0) {
            b = 0;
        }
        if (b >= nbins_) {
            b = nbins_ - 1;
        }
        ++counts_[static_cast<std::size_t>(b)];
        ++total_;
    }
    [[nodiscard]] const std::vector<long>& counts() const {
        return counts_;
    }
    [[nodiscard]] long total() const {
        return total_;
    }
    [[nodiscard]] int nbins() const {
        return nbins_;
    }
    // Bin centre value, for evaluating a reference density at the bin.
    [[nodiscard]] double center(int b) const {
        return lo_ + (static_cast<double>(b) + 0.5) * (hi_ - lo_) / nbins_;
    }

    private:
    double lo_, hi_;
    int nbins_;
    std::vector<long> counts_;
    long total_ = 0;
};

// ---- running mean / standard error -------------------------------------------
class MeanAccumulator {
    public:
    void add(double x) {
        ++n_;
        const double d = x - mean_;
        mean_ += d / static_cast<double>(n_);
        m2_ += d * (x - mean_);
    }
    [[nodiscard]] long n() const {
        return n_;
    }
    [[nodiscard]] double mean() const {
        return mean_;
    }
    [[nodiscard]] double variance() const {
        return n_ > 1 ? m2_ / static_cast<double>(n_ - 1) : 0.0;
    }
    [[nodiscard]] double stderrMean() const {
        return n_ > 0 ? std::sqrt(variance() / static_cast<double>(n_)) : 0.0;
    }

    private:
    long n_ = 0;
    double mean_ = 0.0;
    double m2_ = 0.0;
};

} // namespace stat
} // namespace rtest
