#ifndef RNDARRGENERATOR_H
#define RNDARRGENERATOR_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "DataGenerator.h"
/**
 * \class RndArrGenerator
 * \brief Utilities for random array generation (permutations, pools, multiset/i.i.d. models).
 *
 * \section Description
 * Provides four generators:
 *  (1) a uniformly random permutation of 1..N;
 *  (2) a random array of length N drawn from a k-element value pool that is chosen
 *      uniformly without replacement from [minValue, maxValue] (the pool has unique
 *      values; the output samples from the pool with replacement);
 *  (3) a multiset permutation: draw a uniform positive composition
 *      of N into u parts (stars & bars), build the multiset with labels base..base+u-1,
 *      and uniformly shuffle all elements;

 *
 * Additionally, a helper assembles and shuffles a multiset from an arbitrary profile.
 */
// simple gate: debugif(cond) << "msg" << '\n';
#define debugif(cond) if(cond) std::cerr

class RndArrGenerator: public DataGenerator {
public:
    using DataGenerator::DataGenerator; // default ctor with nondeterministic seed

    /**
    * \brief Generate a uniformly random permutation of integers 1..size.
    *
    * \section Description
    * Constructs the identity sequence [1,2,...,size] and shuffles it in place
    * with the member PRNG. Each permutation is produced with (approximately)
    * equal probability assuming a high-quality uniform PRNG.
    *
    * \section Algorithm
    *  - Allocate vector<int> of length \p size.
    *  - Fill with 1..size via std::iota.
    *  - Apply std::shuffle with the member RNG.
    *
    * \section Notes
    * - Values are 1-based by design. Change the offset
    *   in std::iota if 0-based is desired.
    *
    * \section Complexity
    * - Time:  O(size)
    * - Space: O(1) extra (besides the output vector)
    *
    * @param size Number of elements in the permutation.
    * @return A random permutation of {1, 2, ..., size}.
    *
    */
    std::vector<int> generatePermutation(size_t size) {
        if (size > static_cast<size_t>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument("generatePermutation: size exceeds INT_MAX");
        }
        std::vector<int> data(size);
        std::iota(data.begin(), data.end(), 1);
        std::shuffle(data.begin(), data.end(), rng);
        return data;
    }

    /**
     * \brief Generate a random multiset permutation with exactly \p universe distinct labels (Exact-profile model, §3.1).
     *
     * \section Description
     * Produces an array of length \p total that is a uniformly random linearization of a multiset with
     * exactly \p universe distinct labels (the labels are base, base+1, ..., base+u-1). First, a
     * uniformly random positive composition (profile) x = (x1,...,xu) of total into u parts (each xi ≥ 1)
     * is sampled via stars & bars; then the method builds the multiset
     * { base^x1, (base+1)^x2, ..., (base+u-1)^xu } and shuffles all elements uniformly.
     *
     * \section Algorithm
     *  - Validate inputs: u > 0 and total >= u (each part must get at least one item).
     *  - Draw a positive profile x with randomPositiveProfile(total, universe).
     *  - Return multisetPermutationFromProfile(x, base) which shuffles the whole multiset.
     *
     * \section Complexity
     * - Time:  O(total) overall (assemble O(total) + shuffle O(total)).
     * - Space: O(universe) for the profile plus the output of size total.
     *
     * @param total     Total length n of the output array (n >= u).
     * @param universe  Number of distinct labels u (u >= 1). Exactly u different labels will appear.
     * @param base      Starting label value (labels are base..base+u-1).
     * @return          A uniformly random permutation of the sampled multiset.
     *
     */
    std::vector<int> generateMultisetPermutation(size_t total,
                                                 int universe,
                                                 int base = 1)
    {
        // requires total >= universe so each label appears at least once
        if (universe <= 0) throw std::invalid_argument("generateMultisetPermutation: universe must be > 0");
        if (total < static_cast<size_t>(universe))
            throw std::invalid_argument("generateMultisetPermutation: total must be >= universe");

        std::vector<int> profile = randomPositiveProfile(total, universe);
        return multisetPermutationFromProfile(profile, base);
    }

    /**
     * \brief Generate a random integer array of length \p size using a k-ary value pool.
     *
     * \section Description
     * First builds a pool of exactly @param k unique values chosen uniformly without
     * replacement from the closed interval [@param minValue, @param maxValue]. Then emits
     * an array of length @param size by drawing each element independently and uniformly
     * from that pool (with replacement). The output may contain duplicates; the pool
     * itself contains unique values only.
     *
     * \section Algorithm
     *  - Validate inputs and clamp k to the range width if needed.
     *  - Choose k unique indices from [0..range) with Floyd’s algorithm
     *  - Map indices back to actual values: value = minValue + idx.
     *  - For i=0..size-1: draw a pool index uniformly and write the value.
     *
     * \section Complexity
     * - Time:  O(k) to build the pool + O(size) to generate the output.
     * - Space: O(k) additional for the pool.
     *
     * @param size      Length of the output array.
     * @param minValue  Inclusive lower bound of the value range.
     * @param maxValue  Inclusive upper bound of the value range.
     * @param k         Number of distinct values to include in the pool (k >= 1).
     * @return The generated array of length \p size.
     *
     */
    std::vector<int> generateRandom(size_t size,
                                    int minValue,
                                    int maxValue,
                                    int k) {
        if (size == 0) return {};

        if (minValue > maxValue)
            throw std::invalid_argument("generateRandom: minValue must be <= maxValue");
        if (k <= 0)
            throw std::invalid_argument("generateRandom: k must be positive");

        // Compute inclusive range width using 64-bit to avoid overflow.
        const long long range_ll =
            static_cast<long long>(maxValue) - static_cast<long long>(minValue) + 1LL;
        if (range_ll <= 0)
            throw std::invalid_argument("generateRandom: invalid range (overflow or empty)");
        int range = static_cast<int>(range_ll);

        if (k > range) k = range;

        // Choose k unique indices in [0 .. range), then map to actual values.
        std::vector<int> pickedIdx = floyd_unique_indices(range, k, rng);

        std::vector<int> pool;
        pool.reserve(pickedIdx.size());
        for (int idx : pickedIdx) {
            pool.push_back(minValue + idx);
        }

        // Fill the output by sampling uniformly from the pool with replacement.
        std::vector<int> data(size);
        std::uniform_int_distribution<int> choose(0, static_cast<int>(pool.size()) - 1);
        for (int& v : data) {
            v = pool[choose(rng)];
        }
        return data;
    }

private:

    /**
     * \brief Draw a uniformly random positive profile x of length \p universe summing to \p total (Exact-profile helper).
     *
     * \section Description
     * Implements stars & bars for positive compositions: choose (u-1) distinct cut positions
     * from {1..n-1} uniformly, then take consecutive differences to obtain x1,...,xu (each ≥ 1).
     *
     * \section Algorithm
     *  - Sample k = u-1 unique integers from [0..(n-2)] via Floyd’s algorithm, shift by +1 -> {1..n-1}, sort.
     *  - Let b1<...<bk be the cuts; set x1=b1, x2=b2-b1, ..., xu=n-bk.
     *
     * \section Complexity
     * - Time:  O(universe)
     * - Space: O(universe)
     *
     * @param total     Total sum n (n >= u).
     * @param universe  Number of positive parts u (u >= 1).
     * @return          Profile vector x of length u with sum(x) = n and each x_i >= 1.
     *
     */
    std::vector<int> randomPositiveProfile(size_t total, int universe) {
        if (universe <= 0) throw std::invalid_argument("randomPositiveProfile: universe must be > 0");
        if (total < static_cast<size_t>(universe))
            throw std::invalid_argument("randomPositiveProfile: total must be >= universe");

        if (universe == 1) return { static_cast<int>(total) };

        const int n = static_cast<int>(total);
        const int k = universe - 1;
        std::vector<int> bars0 = floyd_unique_indices(n - 1, k, rng);
        std::sort(bars0.begin(), bars0.end());

        std::vector<int> profile(universe, 0);
        int prev = 0;
        for (int i = 0; i < k; ++i) {
            int b = bars0[i] + 1;
            profile[i] = b - prev;
            prev = b;
        }
        profile[universe - 1] = n - prev;
        return profile;
    }

    /**
     * \brief Assemble and uniformly shuffle a multiset according to the given (nonnegative) profile.
     *
     * \section Description
     * Builds a flat array that contains exactly profile[i] copies of label (base+i), for i=0..u-1,
     * where u = profile.size(), and then shuffles the entire array with the member PRNG. This
     * yields a uniform random permutation of the multiset; all n!/∏(x_i!) linearizations are equiprobable.
     *
     * \section Algorithm
     *  - Validate the profile: non-empty, all counts >= 0; compute total n.
     *  - Reserve n and push label (base+i) exactly profile[i] times for each i.
     *  - Apply std::shuffle to the whole array.
     *
     * \section Complexity
     * - Time:  O(n) where n = sum(profile).
     * - Space: O(1) extra besides the output vector of size n.
     *
     * @param profile   Counts per label (length u). Zeros are allowed.
     * @param base      Starting label value; labels used are base..base+u-1.
     * @return          A uniformly shuffled array realizing the given multiset.
     *
     */
    std::vector<int> multisetPermutationFromProfile(const std::vector<int>& profile, int base) {
        if (profile.empty()) return {};

        long long total_ll = 0;
        for (int c : profile) {
            if (c < 0) throw std::invalid_argument("multisetPermutationFromProfile: negative count");
            total_ll += c;
        }
        if (total_ll == 0) return {};

        std::vector<int> data;
        data.reserve(static_cast<size_t>(total_ll));
        for (size_t i = 0; i < profile.size(); ++i) {
            int label = base + static_cast<int>(i);
            for (int c = 0; c < profile[i]; ++c) data.push_back(label);
        }
        std::shuffle(data.begin(), data.end(), rng);
        return data;
    }

    /**
     * \brief Sample k unique integers from [0, range_size) using Floyd’s algorithm.
     *
     * \section Description
     * Iterates r from range_size - k to range_size - 1; each step draws t uniformly
     * from [0..r]. If t is new, t is inserted; otherwise r is inserted. Produces a
     * uniform k-subset without constructing or shuffling the entire population.
     *
     * \section Algorithm
     * - For r = range_size - k .. range_size - 1:
     *     - Draw t ∈ [0..r].
     *     - If t not in chosen: insert t; else insert r.
     * - Return the collected set as a sorted vector.
     *
     * \section Notes
     * - Sorting the result is optional; it provides a stable order (deterministic mapping without extra sorting).
     *
     * \section Complexity
     * - Expected Time: O(k);
     * - Space: O(k).
     *
     * @param range_size Population size (must be >= k and >= 0).
     * @param k          Sample size (0 <= k <= range_size).
     * @param gen        PRNG.
     * @return Sorted vector of k distinct integers in [0, range_size).
     *
     */
    static std::vector<int> floyd_unique_indices(int range_size, int k, std::mt19937& gen) {
        std::vector<int> out;
        if (range_size <= 0 || k <= 0) return out;
        if (k > range_size) k = range_size;

        std::unordered_set<int> chosen;
        chosen.reserve(static_cast<size_t>(k) * 2);

        std::uniform_int_distribution<int> dis; // rebind per-iteration
        for (int r = range_size - k; r < range_size; ++r) {
            dis = std::uniform_int_distribution<int>(0, r);
            const int v = dis(gen);
            if (!chosen.insert(v).second) {
                chosen.insert(r);
            }
        }
        return out;
    }
};

#endif // RNDARRGENERATOR_H
