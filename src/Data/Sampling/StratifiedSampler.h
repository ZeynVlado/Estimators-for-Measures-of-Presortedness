#ifndef STRATIFIED_SAMPLER_H
#define STRATIFIED_SAMPLER_H
//#define STRATIFIED_SAMPLER_DEBUG

#include <vector>
#include <algorithm>
#include <cassert>
#include <random>
#include <unordered_set>
#include "Sampler.h"
#include "StratBasedSampler.h"
#ifdef STRATIFIED_SAMPLER_DEBUG
#include <iostream>
#endif

/**
 * \class StratifiedSampler
 * \brief Stratified sampling over a 1D array with unbiased quota allocation.
 *
 * \section Description
 * The array is partitioned into consecutive, non-overlapping strata of fixed size
 * (the last stratum may be shorter). The target sample length is distributed across
 * strata so that each stratum contributes at least one element; the remaining quota
 * is then allocated in unbiased random rounds among strata that still have capacity.
 * Inside each stratum, a specified number of unique positions is drawn uniformly
 * at random. Concatenation of per-stratum picks (already increasing by construction)
 * forms the final sample.
 *
 */
class StratifiedSampler : public StratBasedSampler {
public:
    /**
     * \brief Construct a stratified sampler with a target sample length and stratum size.
     *
     * \section Description
     * Stores the parameters needed for stratified sampling. Array partitioning,
     * quota allocation, and element selection are performed in createSample.
     *
     * @param sampleLength Target sample length.
     * @param stratLength  Stratum length used to build fixed strata (last may be shorter).
     *
     */
    StratifiedSampler(int sampleLength, int stratLength)
        : StratBasedSampler(sampleLength), stratLength_(stratLength) {}

    /**
     * \brief Create a stratified sample from @param arr .
     *
     * \section Description
     * Builds fixed strata, assigns unbiased per-stratum quotas (at least one per
     * stratum), selects the required number of unique indices within each stratum,
     * and writes the corresponding values to the output sample. Global index order
     * is increasing by construction.
     *
    * \section Algorithm
    *  1. Partition [0, n) into strata of length \c stratLength_ via StratBasedSampler::buildStrata.
    *  2. Let \c K be the number of strata; assert \c sampleLength_ >= \c K and \c sampleLength_ <= \c n.
    *    3. Compute \c extra = sampleLength_ - K.
    *  4. Distribute per-stratum quotas using
    *     StratBasedSampler::assignExtrasEvenRandomUnique(strata, extra, gen):
    *       - initialize quota_i = 1 for all i,
    *       - while extra > 0, randomly (without replacement per round) add +1 to
    *         distinct active strata until \c extra reaches 0.
    *  5. For each stratum i:
    *       - if quota_i == len_i, take the entire stratum;
    *       - otherwise draw quota_i unique indices uniformly from the stratum
    *         and sort them locally (inside the stratum).
    *  6. Append indices from strata in ascending stratum order; global order is already sorted.
    *  7. Push corresponding values into the sample.
     *
     * \section Complexity
     * - Quota allocation: O(R * K) in the worst case, where R is the number of rounds.
     * - Index selection: for stratum i, O(quota_i) expected with Floyd’s algorithm.
     * - Total: O(R * K + sampleLength_).
     * - Space: O(K) for strata + O(sampleLength_) for output.
     * @param arr Input array.
     *
     */
    void createSample(const std::vector<int>& arr) override {
        sample_.clear();
        sampleIndices_.clear();

        const int n = static_cast<int>(arr.size());
        assert(n > 0);
        assert(sampleLength_ > 0);
        assert(stratLength_ > 0);
        assert(sampleLength_ <= n);

        // Fixed strata (last may be shorter)
        std::vector<Stratum> strata = buildStrata(n, stratLength_);
        const int K = static_cast<int>(strata.size());
        assert(K >= 1);
        assert(sampleLength_ >= K && "sampleLength must be >= number of strata (pick >=1 per stratum)");

        // Unbiased per-round random allocation of the remaining quota
        std::random_device rd;
        std::mt19937 gen(rd());
        const int extra = sampleLength_ - K;
        assignExtrasEvenRandomUnique(strata, extra, gen); // from StratBasedSampler

        // Per-stratum random indices (globally sorted by construction)
        std::vector<int> chosen = selectIndicesByQuotas(strata, gen);

        // Assemble sample (already increasing, no global sort required)
        sampleIndices_ = chosen;
        sample_.reserve(sampleLength_);
        for (int idx : chosen) sample_.push_back(arr[idx]);

#ifdef STRATIFIED_SAMPLER_DEBUG
        debugPrint(n, strata, chosen);
#endif

        assert(static_cast<int>(sampleIndices_.size()) == sampleLength_);
        assert(static_cast<int>(sample_.size())        == sampleLength_);
    }

private:
    int stratLength_;

    /**
     * \brief Select exactly st.quota(quota of stratum) unique indices inside each stratum.
     *
     * \section Description
     * Given a partition of the array into disjoint strata and an integer quota for
     * each stratum (with the quotas summing to the target sample size), this routine
     * constructs a globally ordered index set by sampling, within every stratum,
     * a subset of positions of the prescribed cardinality. Selection within a stratum
     * is uniform over all subsets of the required size; when the quota equals the stratum
     * length, the full stratum domain is admitted. Concatenation across strata follows
     * increasing stratum order, which yields a globally sorted index list without a
     * separate sorting step.
     *
     * \section Algorithm
     * - For each stratum:
     *   - If quota == len: push all indices in [start, end).
     *   - Else: sampleUniqueFromRangeFloyd(len, quota), add start, sort locally,
     *           append to output.
     *
     * \section Notes
     * - Floyd’s algorithm draws m unique integers from  [0..len) in O(m) time on average.
     *
     * \section Complexity
     * - Time O(sampleLength_)
     * - Space proportional to output.
     *
     * @param strata Strata with assigned quotas (sum(quota) == sampleLength_).
     * @param gen    RNG for per-stratum random draws.
     * @return Sorted list of chosen global indices of length sampleLength_.
     *
     */
    static std::vector<int> selectIndicesByQuotas(
        const std::vector<StratBasedSampler::Stratum>& strata,
        std::mt19937& gen
    ) {
        std::vector<int> out;
        int need = 0;
        for (const auto& st : strata) need += st.quota;
        out.reserve(need);

        for (const auto& st : strata) {
            const int m = st.quota;
            if (m == 0) continue;
            if (m == st.len) {
                // Take whole stratum
                for (int x = st.start; x < st.end; ++x) out.push_back(x);
            } else {
                // m unique in [0..len) -> shift by start -> local sort
                auto picks = sampleUniqueFromRangeFloyd(st.len, m, gen);
                for (int& v : picks) v += st.start;
                std::sort(picks.begin(), picks.end());
                out.insert(out.end(), picks.begin(), picks.end());
            }
        }
        return out; // already increasing globally
    }

    /**
     * \brief Draw m unique integers uniformly from [0, len) using Floyd’s algorithm.
     *
     * \section Description
     * Incremental method that iterates j from len - m to len - 1, each time
     * drawing t ~ Uniform{0..j}; if t is new, insert t; otherwise insert j.
     * Produces a uniform m-subset without auxiliary arrays of size len.
     *
     * \section Algorithm
     * - For j = len - m .. len - 1:
     *     - draw t from [0..j];
     *     - if t not yet used -> insert t; else insert j.
     *
     * \section Complexity
     * - Time: O(m).
     * - Space: O(m).
     *
     * @param len Population size.
     * @param m   Sample size (0 <= m <= len).
     * @param gen RNG.
     * @return Unordered set of m distinct integers in [0, len); caller may sort if needed.
     *
     */
    static std::vector<int> sampleUniqueFromRangeFloyd(int len, int m, std::mt19937& gen) {
        assert(m >= 0 && m <= len);
        std::vector<int> out; out.reserve(m);
        std::unordered_set<int> used; used.reserve(m * 2 + 1);
        for (int j = len - m; j < len; ++j) {
            std::uniform_int_distribution<int> d(0, j);
            int t = d(gen);
            if (!used.insert(t).second) used.insert(j);
        }
        out.insert(out.end(), used.begin(), used.end());
        return out;
    }

#ifdef STRATIFIED_SAMPLER_DEBUG
    /**
     * \brief Debug print of stratum metadata and chosen indices.
     */
    static void debugPrint(int n,
                           const std::vector<StratBasedSampler::Stratum>& strata,
                           const std::vector<int>& chosen) {
        std::cerr << "[StratifiedSampler] n=" << n
                  << " strata=" << strata.size()
                  << " sampleLength=" << chosen.size() << "\n";
        for (int i = 0; i < static_cast<int>(strata.size()); ++i) {
            std::cerr << "  stratum#" << i
                      << " start=" << strata[i].start
                      << " end="   << strata[i].end
                      << " len="   << strata[i].len
                      << " quota=" << strata[i].quota << "\n";
        }
        std::cerr << "  chosen indices:";
        for (int v : chosen) std::cerr << ' ' << v;
        std::cerr << "\n";
    }
#endif
};

#endif // STRATIFIED_SAMPLER_H