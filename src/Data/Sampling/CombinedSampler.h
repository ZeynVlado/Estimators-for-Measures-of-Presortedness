

#ifndef COMBINED_SAMPLER_H
#define COMBINED_SAMPLER_H

// #define COMBINED_SAMPLER_DEBUG

#include <vector>
#include <cassert>
#include <random>
#include <utility>
#include "StratBasedSampler.h"


/**
 * \class CombinedSampler
 * \brief Stratified-and-clustered sampling: one contiguous cluster per stratum.
 *
 * \section Description
 * The input array is partitioned into consecutive, non-overlapping strata of fixed
 * size (the last stratum may be shorter). A target sample length is allocated as
 * integer quotas across strata with the constraint that each stratum contributes at
 * least one element. Within each stratum, exactly one contiguous cluster is selected
 * whose length equals the assigned quota for that stratum. Clusters from different
 * strata do not overlap by construction, and the concatenation of all clusters (in
 * stratum order) yields an increasing sequence of indices of the requested length.
 *
 * \section Preconditions / Throws
 * - Asserts (debug): n > 0, stratLength_ > 0, sampleLength_ > 0, sampleLength_ <= n,
 *   K >= 1, sampleLength_ >= K; per-stratum quotas never exceed stratum lengths;
 *   final sample size equals sampleLength_. No exceptions thrown.
 *
 */
class CombinedSampler : public StratBasedSampler {
public:
    /**
     * \brief Construct a combined sampler with a target sample length and stratum size.
     *
     * @param sampleLength Target sample length.
     * @param stratLength  Stratum length used to build fixed strata (last may be shorter).
     *
     * /section Preconditions / Throws
     * - No exceptions thrown. Preconditions are asserted in \ref createSample.
     */
    CombinedSampler(int sampleLength, int stratLength)
        : StratBasedSampler(sampleLength), stratLength_(stratLength) {}

    /**
     * \brief Create a combined stratified–clustered sample from \p arr.
     *
     * \section Description
     * Builds fixed strata, assigns unbiased per-stratum quotas (at least one per stratum),
     * and for each stratum selects exactly one contiguous cluster whose length equals
     * the quota of that stratum. Emitted indices are globally sorted by construction.
     *
     * \section Algorithm
     *  1. Partition [0, n) into K strata of length stratLength_ via
     *     StratBasedSampler::buildStrata.
     *  2. Check feasibility: sampleLength_ <= n and sampleLength_ >= K.
     *  3. Let extra = sampleLength_ - K. Distribute quotas across strata using
     *     StratBasedSampler::assignExtrasEvenRandomUnique(strata, extra, gen):
     *       - initialize quota_i = 1 for all i;
     *       - while extra > 0, randomly (without replacement per round) add +1
     *         to distinct active strata (those with remaining capacity).
     *  4. For each stratum with quota m:
     *       - if m == len, choose the whole stratum [start, end);
     *       - else choose a start offset uniformly from [0 .. len - m] and form a single
     *         cluster [s, s+m).
     *  5. Emit elements of all per-stratum clusters in stratum order; indices are already
     *     globally increasing.
     *
     * \section Complexity
     * - Quota assignment: O(R * K) in the worst case, where R is the number of rounds.
     * - Cluster picking per stratum: O(1) for offset selection + O(m) to output m elements.
     * - Total: O(R * K + sampleLength_).
     * - Space: O(K) for strata + O(sampleLength_) for output.
     *
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

        // 1) Strata (last may be shorter)
        std::vector<Stratum> strata = buildStrata(n, stratLength_);
        const int K = static_cast<int>(strata.size());
        assert(K >= 1);
        assert(sampleLength_ >= K && "sampleLength must be >= number of strata");

        // 2) Unbiased per-round random allocation of the remaining quota
        const int extra = sampleLength_ - K;
        std::random_device rd;
        std::mt19937 gen(rd());
        assignExtrasEvenRandomUnique(strata, extra, gen);

#ifdef COMBINED_SAMPLER_DEBUG
        std::vector<std::pair<int,int>> clustersChosen;
        clustersChosen.reserve(K);
#endif

        for (const auto& st : strata) {
            if (st.quota <= 0) continue;          // by construction, quota >= 1
            assert(st.quota <= st.len);
            const auto [s, e] = pickClusterInStratum(st, gen);

#ifdef COMBINED_SAMPLER_DEBUG
            clustersChosen.emplace_back(s, e);
#endif
            for (int i = s; i < e; ++i) {
                sampleIndices_.push_back(i);
                sample_.push_back(arr[i]);
            }
        }

#ifdef COMBINED_SAMPLER_DEBUG
        debugPrint(n, strata, clustersChosen, sampleIndices_);
#endif

        assert(static_cast<int>(sampleIndices_.size()) == sampleLength_);
        assert(static_cast<int>(sample_.size())        == sampleLength_);
    }

private:
    int stratLength_;

    /**
     * \brief Select one contiguous cluster of length equal to the stratum quota.
     *
     * \section Algorithm
     * - Let len = end - start and m = quota.
     * - If m == 0: return {start, start}.
     * - If m == len: return {start, end}.
     * - Else: draw off ~ Uniform{0..len-m}; set s = start + off; return {s, s + m}.
     *
     * \section Complexity
     * - Time: O(1). Space: O(1).
     *
     * @param st  Stratum (start, end, len, quota).
     * @param gen Pseudorandom generator for offset sampling.
     * @return A pair (s, e) representing the chosen cluster [s, e).
     *
     */
    static std::pair<int,int> pickClusterInStratum(const Stratum& st, std::mt19937& gen) {
        const int m = st.quota;
        assert(m >= 0 && m <= st.len);
        if (m == 0)      return {st.start, st.start};
        if (m == st.len) return {st.start, st.end};
        std::uniform_int_distribution<int> d(0, st.len - m);
        const int off = d(gen);
        const int s = st.start + off;
        return {s, s + m};
    }

#ifdef COMBINED_SAMPLER_DEBUG
    /**
     * \brief Debug print of per-stratum metadata, chosen clusters, and resulting indices.
     */
    static void debugPrint(
        int n,
        const std::vector<StratBasedSampler::Stratum>& strata,
        const std::vector<std::pair<int,int>>& clustersChosen,
        const std::vector<int>& chosen
    ){
        std::cerr << "[CombinedSampler] n=" << n
                  << " strata=" << strata.size()
                  << " sampleLength=" << chosen.size() << "\n";
        for (int i = 0; i < static_cast<int>(strata.size()); ++i) {
            const auto& st = strata[i];
            std::cerr << "  stratum#" << i
                      << " start=" << st.start
                      << " end="   << st.end
                      << " len="   << st.len
                      << " quota=" << st.quota;
            if (i < static_cast<int>(clustersChosen.size())) {
                std::cerr << "  cluster=[" << clustersChosen[i].first
                          << "," << clustersChosen[i].second << ")";
            }
            std::cerr << "\n";
        }
        std::cerr << "  chosen indices:";
        for (int v : chosen) std::cerr << ' ' << v;
        std::cerr << "\n";
    }
#endif
};

#endif // COMBINED_SAMPLER_H