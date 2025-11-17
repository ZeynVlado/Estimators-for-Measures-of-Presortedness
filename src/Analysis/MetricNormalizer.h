

#define METRICNORMALIZER_H

/// \class MetricNormalizer
/// \brief Normalizes disorder metrics to [0,1] (singleton).
///
/// Provides scaling of raw metric values to [0,1] using natural upper bounds.
class MetricNormalizer {
public:
    /// \brief Access the singleton instance (thread-safe).
    static MetricNormalizer& instance() {
        static MetricNormalizer inst;
        return inst;
    }

    MetricNormalizer(const MetricNormalizer&) = delete;
    MetricNormalizer& operator=(const MetricNormalizer&) = delete;
    MetricNormalizer(MetricNormalizer&&) = delete;
    MetricNormalizer& operator=(MetricNormalizer&&) = delete;

    /**
     * \brief Inversions normalization.
     * \details Formula:  inv_norm = invCount / C(n,2) = invCount / (n(n-1)/2).
     */
    double normalizeInversions(long long invCount, long long n) const {
        if (n < 2) return 0.0;
        long double nd = static_cast<long double>(n);
        long double denom = nd * (nd - 1.0L) / 2.0L;
        return static_cast<double>(static_cast<long double>(invCount) / denom);
    }

    /**
     * \brief Runs normalization.
     * \details Formula: runs_norm = (runsCount - 1) / (n - 1), for n >= 2.
     */
    double normalizeRuns(long long runsCount, long long n) const {
        if (n < 2) return 0.0;
        return static_cast<double>(
            (static_cast<long double>(runsCount) - 1.0L) /
            (static_cast<long double>(n) - 1.0L)
        );
    }

    /**
     * \brief Run-length (Shannon) entropy normalization.
     * \details Formula: H_norm = H / ln(R), where R = runsCount and H is natural-log entropy.
     */
    double normalizeRunsEntropy(double H, long long runsCount) const {
        if (runsCount <= 1) return 0.0;
        return H / std::log(static_cast<long double>(runsCount)); // [0,1]
    }

    /**
     * \brief Distinct-count normalization.
     * \details Formula: distinct_norm = (U - 1) / (n - 1), where U is the number of distinct values.
     */
    double normalizeDistinct(long long distinctCount, long long n) const {
        if (n < 2) return 0.0;
        long double num = static_cast<long double>(distinctCount) - 1.0L;
        long double den = static_cast<long double>(n) - 1.0L;
        if (den <= 0.0L) return 0.0;
        long double v = num / den;
        if (v < 0.0L) v = 0.0L;
        if (v > 1.0L) v = 1.0L;
        return static_cast<double>(v);
    }

    /**
     * \brief Value-frequency (Shannon) entropy normalization.
     * \details Formula: H_norm = H / ln(U), where U is the number of distinct values.
     */
    double normalizeValueFrequencyEntropy(double H, long long distinctCount) const {
        if (distinctCount <= 1) return 0.0;
        const long double denom = std::log(static_cast<long double>(distinctCount));
        if (!(denom > 0.0L)) return 0.0;
        long double v = static_cast<long double>(H) / denom;
        if (v < 0.0L) v = 0.0L;
        if (v > 1.0L) v = 1.0L;
        return static_cast<double>(v);
    }

    /**
     * \brief Hartley entropy normalization.
     * \details Formula: H0_norm = H0 / ln(n), where H0 = ln(U) and U is distinct count.
     */
    double normalizeDistinctElementEntropy(double H0, long long n) const {
        if (n <= 1) return 0.0;
        const long double denom = std::log(static_cast<long double>(n));
        return (denom > 0.0L) ? static_cast<double>(H0 / denom) : 0.0;
    }

    /**
     * \brief Removals-to-LIS normalization.
     * \details Formula: rem_norm = removals / (n - 1), where removals = n - LIS.
     */
    double normalizeRem(long long removals, long long n) const {
        if (n < 2) return 0.0;
        return static_cast<double>(
            static_cast<long double>(removals) / (static_cast<long double>(n) - 1.0L)
        );
    }

    /**
     * \brief Oscillations normalization.
     * \details Formula: osc_norm = oscCount / (n - 2), for n >= 3.
     */
    double normalizeOsc(long long oscCount, long long n) const {
        if (n < 3) return 0.0;
        return static_cast<double>(
            static_cast<long double>(oscCount) / (static_cast<long double>(n) - 2.0L)
        );
    }

    /**
     * \brief Spearman footrule (sum of absolute rank displacements) normalization.
     * \details Formula:
     *  - If n even:   dis_norm = disSum / (n^2 / 2)
     *  - If n odd:    dis_norm = disSum / ((n^2 - 1) / 2)
     */
    double normalizeDis(long long disSum, long long n) const {
        if (n <= 1) return 0.0;
        const long double nd = static_cast<long double>(n);
        const long double denom = (n % 2 == 0)
            ? (nd * nd / 2.0L)
            : ((nd * nd - 1.0L) / 2.0L);
        return static_cast<double>(static_cast<long double>(disSum) / denom);
    }

    /**
     * \brief Hamming distance-to-sorted normalization.
     * \details Formula: ham_norm = hamCount / n.
     */
    double normalizeHam(long long hamCount, long long n) const {
        if (n <= 0) return 0.0;
        return static_cast<double>(static_cast<long double>(hamCount) / static_cast<long double>(n));
    }

    /**
     * \brief Maximum displacement normalization.
     * \details Formula: max_norm = maxDisplacement / (n - 1), for n >= 2.
     */
    double normalizeMax(long long maxDisplacement, long long n) const {
        if (n < 2) return 0.0;
        return static_cast<double>(
            static_cast<long double>(maxDisplacement) / (static_cast<long double>(n) - 1.0L)
        );
    }

private:
    MetricNormalizer() = default;
};