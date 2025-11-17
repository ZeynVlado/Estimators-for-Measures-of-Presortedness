#ifndef STRAT_BASED_SAMPLER_H
#define STRAT_BASED_SAMPLER_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <cassert>
#include "Sampler.h"

/**
 * \class StratBasedSampler
 * \brief Base class for samplers that operate on contiguous strata (blocks).
 *
 * \section Description
 * Provides common primitives for stratified sampling schemes that partition
 * an array into consecutive strata and then allocate a per-stratum quota.
 * Includes: (1) fixed-size stratum construction, (2) quota initialization
 * and capacity accounting, and (3) unbiased  quota–distribution policy.
 *
 */
class StratBasedSampler : public Sampler {
public:
    /**
     * \brief Construct the base class with a requested sample length.
     *
     * @param sampleLength Requested sample length;
     */
    explicit StratBasedSampler(int sampleLength)
        : Sampler(sampleLength) {}

    ~StratBasedSampler() override = default;

protected:
    /**
     * \struct Stratum
     * \brief Half-open index range [start, end) with an assigned quota.
     *
     * \section Description
     * Represents a contiguous block of the array. The quota field indicates
     * how many indices should be sampled from this stratum.
     */
    struct Stratum {
        int start = 0;
        int end   = 0;   ///< exclusive
        int len   = 0;   ///< end - start
        int quota = 0;   ///< number of items to take from this stratum
    };

    /**
     * \brief Partition [0, n) into consecutive strata of size \p stratSize (last may be shorter).
     *
     * \section Description
     * Produces K = ceil(n / stratSize) strata that tile the index range without
     * overlap: the first K-1 strata have length \p stratSize, the last stratum
     * has length n % stratSize (or \p stratSize if divisible).
     *
     * \section Algorithm
     * 1. Compute K = (n + stratSize - 1) / stratSize.
     * 2. For i in [0..K-1]:
     *    - start = i * stratSize
     *    - end   = min(start + stratSize, n)
     *    - emit Stratum{start, end, end-start, quota=0}
     *
     * \section Notes
     * - Strata are returned in ascending order and are disjoint, covering [0, n).
     *
     * \section Complexity
     * - Time: O(K).
     * - Space: O(K).
     *
     * @param n         Total array length.
     * @param stratSize Target stratum size (> 0).
     * @return Vector of consecutive strata.
     *
     */
    static std::vector<Stratum> buildStrata(int n, int stratSize) {

        const int K = (n + stratSize - 1) / stratSize;
        std::vector<Stratum> out; out.reserve(K);
        for (int i = 0; i < K; ++i) {
            const int start = i * stratSize;
            const int end   = std::min(start + stratSize, n);
            out.push_back(Stratum{ start, end, end - start, 0 });
        }
        return out;
    }


    /**
     * \brief Unbiased random distribution of the extra quota (random, unique per round).
     *
     * \section Description
     * Initializes each stratum with quota = 1, then repeatedly selects a random
     * subset of the currently “active” strata (those with positive capacity
     * cap_i = len_i - 1). Each round assigns at most one increment per active stratum (no repeats
     * within the round), yielding unbiased allocation over multiple rounds.
     * This prevents positional bias.
     *
     * \section Algorithm
     * 1. Call initBaseQuotas(strata) -> quota_i = 1 for all i.
     * 2. Compute cap = computeCapacities(strata) where cap_i = max(0, len_i - 1).
     * 3. While extra > 0:
     *    3.1 active = collectActive(cap), assert nonempty.
     *    3.2 Shuffle active uniformly.
     *    3.3 t = min(extra, |active|); for j in [0..t-1]:
     *        - i = active[j]; assert cap_i > 0; quota_i += 1; cap_i -= 1.
     *    3.4 extra -= t.
     *
     *
     * \section Complexity
     * - Time: O(R * K) in the worst case (R = ceil(extra / |active|)), Space: O(K).
     *
     * \param strata In/out: list of strata; quotas are modified in place.
     * \param extra  Total additional units to distribute after giving each stratum 1.
     * \param gen    Random number generator used for shuffling.
     *
     */
    static void assignExtrasEvenRandomUnique(std::vector<Stratum>& strata,
                                             int extra,
                                             std::mt19937& gen) {
        initBaseQuotas(strata);                 // quota_i = 1
        if (extra <= 0) return;

        std::vector<int> cap = computeCapacities(strata); // cap_i = len_i - 1
        const long long totalCap = sum(cap);
        assert(totalCap >= extra && "Not enough capacity across strata");

        while (extra > 0) {
            std::vector<int> active = collectActive(cap);
            assert(!active.empty() && "No active strata to allocate extras");

            std::shuffle(active.begin(), active.end(), gen);
            const int t = std::min(extra, static_cast<int>(active.size()));

            for (int j = 0; j < t; ++j) {
                const int idx = active[j];
                assert(cap[idx] > 0);
                strata[idx].quota += 1;
                cap[idx]          -= 1;
            }
            extra -= t;
        }
    }

    /**
     * \brief Initialize per-stratum quotas to 1.
     *
     * \section Description
     * Sets quota_i = 1 for each stratum. This enforces the constraint that at
     * least one element must be taken from every stratum before distributing
     * any additional units of quota.
     *
     * \section Algorithm
     * - For each stratum, set quota = 1 (assert len > 0).
     *
     * \section Complexity
     * - Time: O(K). Space: O(1) besides the input vector.
     *
     * @param strata In/out: list of strata.
     *
     */
    static void initBaseQuotas(std::vector<Stratum>& strata) {
        for (auto& st : strata) {
            assert(st.len > 0);
            st.quota = 1;
        }
    }

    /**
     * \brief Compute per-stratum capacities after the base allocation (quota = 1).
     *
     * \section Description
     * Capacity for stratum i is the maximum extra quota it can accept beyond the
     * base 1, i.e., cap_i = max(0, len_i - 1). Returns the vector of capacities.
     *
     * \section Algorithm
     * - For each stratum, push back max(0, len - 1).
     *
     * \section Complexity
     * - Time: O(K).
     * - Space: O(K).
     *
     * @param strata List of strata (len must be known).
     * @return Vector cap of size K where cap[i] = max(0, len_i - 1).
     */
    static std::vector<int> computeCapacities(const std::vector<Stratum>& strata) {
        std::vector<int> cap;
        cap.reserve(strata.size());
        for (const auto& st : strata) cap.push_back(std::max(0, st.len - 1));
        return cap;
    }

    /**
     * \brief Sum over a vector of numeric values.
     */
    template <class T>
    static long long sum(const std::vector<T>& v) {
        long long s = 0;
        for (auto x : v) s += x;
        return s;
    }

    /**
     * \brief Collect indices of strata with positive remaining capacity.
     *
     * \section Description
     * Returns the list of indices i where cap_i > 0, i.e., strata that can
     * receive at least one more unit of quota.
     *
     * \section Algorithm
     * - Iterate over cap; push back i whenever cap[i] > 0.
     *
     * \section Complexity
     * - Time: O(K). Space: O(K) for the output index list.
     *
     * @param cap Vector of per-stratum capacities.
     * @return Indices of active strata.
     */
    static std::vector<int> collectActive(const std::vector<int>& cap) {
        std::vector<int> active; active.reserve(cap.size());
        for (int i = 0; i < static_cast<int>(cap.size()); ++i)
            if (cap[i] > 0) active.push_back(i);
        return active;
    }

    /**
     * \brief Lift all active strata by a common base increment (bounded by capacity).
     *
     * \section Description
     * Adds \p base to each active stratum’s quota, but not exceeding individual
     * capacity; returns the total amount actually used (≤ base * |active|).
     *
     * \section Algorithm
     * - For each index in active:
     *    - add = min(base, cap[i]); quota[i] += add; cap[i] -= add; used += add.
     *
     * \section Complexity
     * - Time: O(|active|). Space: O(1) besides inputs.
     *
     * @param strata In/out: strata whose quotas are updated.
     * @param cap    In/out: capacities decreased in place.
     * @param active Indices of strata with cap>0.
     * @param base   Common increment to attempt for all active strata.
     * @return Total increment consumed across active strata.
     */
    static int applyBaseToActive(std::vector<Stratum>& strata,
                                 std::vector<int>& cap,
                                 const std::vector<int>& active,
                                 int base) {
        int used = 0;
        for (int idx : active) {
            if (cap[idx] <= 0) continue;
            const int add = std::min(base, cap[idx]);
            strata[idx].quota += add;
            cap[idx]          -= add;
            used              += add;
        }
        return used;
    }

    /**
     * \brief Assign at most one additional unit to each active stratum (single pass).
     *
     * \section Description
     * Distributes up to \p extra units by giving +1 to active strata left-to-right,
     * skipping those whose capacity is already zero. Returns how many units were
     * actually used during this pass.
     *
     * \section Algorithm
     * - For i in active (in order):
     *    - if extra == 0, stop.
     *    - if cap[i] > 0: quota[i] += 1; cap[i] -= 1; ++used; --extra.
     *
     * \section Notes
     * - Used by the deterministic policy to allocate the per-round remainder.
     *
     * \section Complexity
     * - Time: O(|active|). Space: O(1) besides inputs.
     *
     * @param strata In/out: strata whose quotas are updated.
     * @param cap    In/out: capacities decreased in place.
     * @param active Indices of strata with cap>0.
     * @param extra  Upper bound on increments to be assigned in this pass.
     * @return Number of units actually assigned (0..extra).
     */
    static int applyRemainderOnce(std::vector<Stratum>& strata,
                                  std::vector<int>& cap,
                                  const std::vector<int>& active,
                                  int extra) {
        int used = 0;
        for (int idx : active) {
            if (extra == 0) break;
            if (cap[idx] > 0) {
                strata[idx].quota += 1;
                cap[idx]          -= 1;
                --extra;
                ++used;
            }
        }
        return used;
    }
};

#endif // STRAT_BASED_SAMPLER_H