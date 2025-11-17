
#ifndef RUNGENERATORNOISEBASED_H
#define RUNGENERATORNOISEBASED_H
#include <assert.h>
#include <stdexcept>


#include "PresortednessMetrics.h"

#include "RunArrGenerator.h"


/**
 * \class RunGeneratorNoiseBased
 * \brief Build a run-structured array by injecting strict drops (“noise”) into a sorted base.
 *
 * \section Description
 * Produces an integer array of length N that contains exactly a target number of
 * strictly ascending runs. The method starts from a strictly increasing base sequence
 * [1..N] (one run), then randomly selects non-adjacent positions and replaces those
 * elements with values small enough to create strict drop boundaries. Each injected
 * drop adds exactly one run; using k drops yields 1+k runs.
 *
 * This implementation avoids precomputing “candidates” explicitly. The feasible
 * positions form a simple contiguous index range that can be derived analytically
 * from \p minValue and the strictly increasing base. Random, unique, non-adjacent
 * positions are then drawn directly and applied in-place.
 */
class RunGeneratorNoiseBased : public RunArrGenerator {
public:
    using RunArrGenerator::RunArrGenerator;

    /**
     * \brief Main pipeline: build an array with exactly \p runsTarget strictly ascending runs.
     *
     * \section Description
     * Let the strictly increasing base be a[i] = i+1 for i in [0..N-1].
     * A strict drop at index i (1 <= i <= N-2) is feasible iff there exists a value
     * in [minValue .. min(a[i-1], a[i+1]) - 1]. For the strictly increasing base,
     * this reduces to a[i-1] > minValue, i.e. i > minValue. Hence the feasible
     * index range is:
     *   i ∈ [ firstIdx .. N-2 ], where firstIdx = max(1, minValue+1)   (0-based).
     *
     * The method:
     *  - builds the strictly increasing base,
     *  - computes the feasible non-adjacent capacity over that contiguous range,
     *  - draws exactly (runsTarget-1) non-adjacent indices at random (unique, no repeats),
     *  - assigns at each chosen index a value uniformly from [minValue .. min(a[i-1],a[i+1])-1],
     *  - returns the modified array which now has exactly \p runsTarget runs.
     *
     * \section Algorithm
     * 1) Handle trivial cases:
     *    - if N == 0: require runsTarget == 0, else throw.
     *    - if N > 0: require runsTarget in [1..maxRuns], where maxRuns is computed below.
     * 2) data <- generateSortedUnique(N)  // 1..N
     * 3) dropCount <- runsTarget - 1; if dropCount == 0: return data.
     * 4) firstIdx <- max(1, minValue+1); lastIdx <- N-2.
     * 5) capacity <- non-adjacent capacity of [firstIdx..lastIdx], equals ceil(L/2) where L is the range length.
     *    - If dropCount > capacity: throw std::runtime_error.
     * 6) drops <- chooseNonAdjacentFromRangeExact(firstIdx, lastIdx, dropCount, rng)
     *    (random, unique, non-adjacent; no sorting required).
     * 7) For each i in drops:
     *      left  = data[i-1], right = data[i+1], hi = min(left, right) - 1 (== left-1 for strictly increasing base).
     *      draw v uniformly from [minValue .. hi]; data[i] = v.
     * 8) Return data.
     *
     * \section Notes
     * - Non-adjacent constraint guarantees that drops do not interfere with each other’s
     *   feasibility and that each adds exactly one run.
     * - Drops are selected randomly and uniquely per construction; no post-sorting is needed.
     * - Duplicates across different runs are expected and harmless; within each run values remain strictly increasing.
     *
     * \section Complexity
     * - Time: O(N) to build the base + O(dropCount) to apply the drops.
     * - Space: O(N) for the result; O(dropCount) temporary storage.
     *
     *
     * @param N           Total array length.
     * @param runsTarget  Target number of strictly ascending runs (R).
     * @param minValue    Inclusive lower bound for injected values.
     * @param maxValue    (Unused here) Upper bound of the domain; the base is [1..N].
     */
    std::vector<int> generateRunsFromSortedArrayNoise(size_t n,
                                                      int runsTarget,
                                                      int minValue,
                                                      int /*maxValue*/)
    {
        assert(n >= runsTarget);
        assert(runsTarget > 0);

        // Base: strictly increasing (one run)
        std::vector<int> data = generateSortedUnique(n);

        const int dropCount = runsTarget - 1;
        if (dropCount == 0) return data; // already 1 run

        // Feasible contiguous range for drops: i in [firstIdx..lastIdx] (0-based)
        // Condition for feasibility on strictly increasing base: a[i-1] > minValue -> i > minValue.
        const int firstIdx = std::max<int>(1, (minValue + 1));
        if (n < 3 || firstIdx > n - 2) {
            // No internal feasible positions exist
            throw std::runtime_error("Not enough feasible positions to place any drop");
        }
        const int lastIdx = n - 2;
        const int rangeLen = (lastIdx - firstIdx + 1);

        // Non-adjacent capacity on a linear block of length L is ceil(L/2)
        const size_t capacity = (rangeLen + 1) / 2;
        assert(capacity >= dropCount);

        // Draw exactly dropCount non-adjacent indices in [firstIdx..lastIdx]
        std::vector<size_t> drops = chooseNonAdjacentFromRangeExact(firstIdx, lastIdx, dropCount, this->rng);

        // Apply drops in-place (no need to sort)
        applyNoiseInsertion(data, drops, minValue);

#ifdef RUNGENERATORNOISEBASED_DEBUG
        DisorderMetrics dbg;
        assert(dbg.calculateRuns(data) == runsTarget);
#endif
        return data;
    }

private:
    /**
     * \brief Choose exactly \p k non-adjacent indices from [lo..hi] (inclusive), randomly and uniquely.
     *
     * \section Description
     * Builds a maximal non-adjacent set on [lo..hi] by taking every other index
     * starting at \p lo, which yields capacity = ceil(L/2) where L = hi-lo+1.
     * If capacity > k, randomly remove (capacity - k) distinct entries from this
     * set (preserving order of the remaining ones). The result is non-adjacent
     * by construction and has exactly \p k elements.
     *
     * \section Algorithm
     * - Build vector S = { lo, lo+2, lo+4, ... } up to <= hi.
     * - If |S| == k: return S.
     * - Else:
     *    - Generate indices [0..|S|-1], shuffle, take first (|S|-k) as removal positions,
     *      sort removal positions descending and erase from S to keep remaining order.
     *    - Return S (size == k).
     *
     * \section Complexity
     * - Time: O(L) to build S + O(|S|) to remove extras.
     * - Space: O(|S|).
     *
     */
    static std::vector<size_t>
    chooseNonAdjacentFromRangeExact(int lo, int hi, int k, std::mt19937& rng)
    {

        assert(lo <= hi);
        const int L = hi - lo + 1;
        const int cap = (L + 1) / 2;
        assert(cap >= k);

        // Maximal non-adjacent set: take every other index starting at lo
        std::vector<size_t> S; S.reserve(cap);
        for (size_t x = lo; x <= hi; x += 2) S.push_back(x);

        if (S.size() == k) return S;

        const size_t toRemove = S.size() - k;
        std::vector<size_t> idx(S.size());
        std::iota(idx.begin(), idx.end(), size_t{0});
        std::shuffle(idx.begin(), idx.end(), rng);
        idx.resize(toRemove);
        std::sort(idx.begin(), idx.end(), std::greater<size_t>());
        for (size_t j : idx) {
            S.erase(S.begin() + static_cast<std::ptrdiff_t>(j));
        }
        return S;
    }

    /**
     * \brief Insert strict drops at the specified positions.
     *
     * \section Description
     * For each i in \p positions, replace data[i] with a value drawn uniformly from
     * [minValue .. min(data[i-1], data[i+1]) - 1]. On a strictly increasing base,
     * this simplifies to [minValue .. data[i-1]-1], ensuring data[i] < data[i-1]
     * while also data[i] < data[i+1] holds automatically.
     *
     * \section Complexity
     * - Time: O(|positions|). Space: O(1).
     *
     */
    void applyNoiseInsertion(std::vector<int>& data,
                             const std::vector<size_t>& positions,
                             int minValue)
    {
        for (size_t i : positions) {
            int left = data[i-1], right = data[i+1];
            int hi = std::min(left, right) - 1;
            if (hi < minValue) throw std::runtime_error("");
            std::uniform_int_distribution<int> dist(minValue, hi);
            data[i] = dist(this->rng);
        }
    }

};
#endif // RUNGENERATORNOISEBASED_H