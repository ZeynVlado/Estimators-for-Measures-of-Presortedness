

#ifndef CLUSTERSAMPLER_H
#define CLUSTERSAMPLER_H

//#define CLUSTER_SAMPLER_DEBUG

#include <vector>
#include <algorithm>
#include <random>
#include <cassert>
#include <iostream>
#include <cmath>
#include "Sampler.h"
/**
 * \class ClusterSampler
 * \brief Cluster-based sampling over a 1D array.
 *
 * \section Description
 * This sampler builds a sample by covering the array with a small number of
 * non-overlapping contiguous segments (“clusters”) and taking elements only
 * from within those segments. First, we place equally sized base clusters
 * uniformly at random across the array without overlap. Then, if the requested
 * sample is not fully covered by those clusters, we enlarge some clusters to
 * the right by one element at a time, choosing which clusters to extend using
 * unbiased random rounds. Finally, we concatenate all covered elements (in
 * order) to obtain the sample of the desired size.
 *
 */
class ClusterSampler : public Sampler {
public:
    /**
     * \struct Cluster
     * \brief Half-open index range [start, end) representing a contiguous segment.
     */
    struct Cluster { int start; int end; }; // [start, end)

    /**
     * \brief Construct a ClusterSampler.
     *
     * \section Description
     * Stores the requested sample length and base cluster size. The base size
     * caps the initial (equal) length of all clusters before distributing any
     * remainder.
     *
     * @param sampleLength Desired sample length (upper-bounded by array size).
     * @param clusterSize  Desired base cluster size (must be <= sampleLength).
     *
     */
    ClusterSampler(int sampleLength, int clusterSize)
        : Sampler(sampleLength), clusterSize_(clusterSize) {}

    /**
     * \brief Build the sample from \p arr using the cluster-based strategy.
     *
     * \section Description
     * Constructs non-overlapping base clusters placed uniformly at random,
     * possibly extends them to match exactly the requested sample size, and
     * returns the concatenated elements in ascending index order.
     *
     * \section Algorithm
     * 1. Validate preconditions; clear output buffers.
     * 2. Compute n = arr.size(), targetSample = min(sampleLength_, n).
     * 3. Set clusterLength = min(clusterSize_, targetSample).
     * 4. Compute clusterCount = floor(targetSample / clusterLength) (>= 1).
     * 5. Place clusters via layout_random_stars_bars.
     * 6. Let baseCovered = clusterCount * clusterLength; reminder = targetSample - baseCovered.
     * 7. If reminder > 0:
     *      7.1 Compute total right-extension capacity via right_extension_capacity.
     *      7.2 Assert reminder <= capacity.
     *      7.3 Distribute with distribute_reminder_random (unbiased, unique per round).
     * 8. Append all cluster elements into sample_/sampleIndices_ in ascending order.
     * 9. Assert final size equals targetSample.
     *
     * \section Notes
     * - Stars&bars makes all admissible gap splits equally likely.
     * - Random unique-per-round extension avoids left-to-right bias.
     *
     * \section Complexity
     * - Time: O(k log k + targetSample), where k = clusterCount
     *   (stars&bars sorting is O(k log k); copying elements is O(targetSample)).
     * - Space: O(k) for clusters + O(targetSample) for outputs.
     *
     * @param arr Input array to sample from.
     */
    void createSample(const std::vector<int>& arr) override {
        sample_.clear();
        sampleIndices_.clear();

        const int n = static_cast<int>(arr.size());

        // Preconditions
        assert(n > 0 && "Array must be non-empty");
        assert(n > sampleLength());
        assert(clusterSize_ > 0 && "clusterSize must be > 0");
        assert(sampleLength_ > 0 && "sampleLength must be > 0");
        assert(clusterSize_ <= sampleLength_ && "clusterSize must be <= sampleLength");

        const int targetSample  = std::min(sampleLength_, n);
        const int clusterLength = std::min(clusterSize_, targetSample);

        // clusterCount = floor(targetSample / clusterLength)
        const int clusterCount = target_sample_cluster_count(targetSample, clusterLength);
        assert(clusterCount >= 1);

        std::random_device rd;
        std::mt19937 gen(rd());

        // Random non-overlapping clusters via Stars & Bars
        std::vector<Cluster> clusters =
            layout_random_stars_bars(n, clusterLength, clusterCount, gen);

        // Compute reminder
        const int baseCovered = clusterCount * clusterLength;
        assert(baseCovered <= targetSample);
        int reminder = targetSample - baseCovered;

        // Random, unbiased reminder distribution
        if (reminder > 0) {
            const int capacity = right_extension_capacity(clusters, n);
            assert(reminder <= capacity &&
                   "Not enough right-side gap to place reminder safely");
            distribute_reminder_random(reminder, n, clusters, gen);
        }

#ifdef CLUSTER_SAMPLER_DEBUG
        debug_print_clusters(clusters, targetSample);
#endif

        // Concatenate in ascending order (clusters are already sorted by start)
        for (const auto& cl : clusters) {
            for (int i = cl.start; i < cl.end; ++i) {
                sample_.push_back(arr[i]);
                sampleIndices_.push_back(i);
            }
        }

#ifdef CLUSTER_SAMPLER_DEBUG
        std::cerr << "[ClusterSampler] sample.size() = " << sample_.size() << '\n';
#endif

        // Final size assertion
        assert(static_cast<int>(sample_.size()) == targetSample &&
               "Sample size must equal requested targetSample");
    }

private:

    /**
     * \brief Compute how many full clusters fit into \p targetSample.
     *
     * \section Description
     * Returns the integer number of equal-length base clusters that can be
     * placed such that their total length does not exceed the requested sample.
     *
     * @param targetSample Effective requested sample length (<= n).
     * @param clusterLength Fixed base cluster length (>= 1).
     * @return floor(targetSample / clusterLength).
     *
     * \section Preconditions / Throws
     * - Asserts: clusterLength > 0.
     */
    static int target_sample_cluster_count(int targetSample, int clusterLength) {
        assert(clusterLength > 0);
        return targetSample / clusterLength; // floor
    }


    /**
     * \brief Stars & Bars layout: place k non-overlapping clusters of length L uniformly at random in [0, n).
     *
     * \section Description
     * All valid placements (positions of k equal-length, non-overlapping blocks) are equiprobable.
     * Model: arrange G = n - k*L “stars” and k “bars” uniformly among G+k slots (classical stars & bars).
     * Choosing k bar positions without replacement in [0 .. G+k-1] yields a uniform composition
     * (g_0,...,g_k) of G, where g_i is the gap size before/between/after clusters.
     *
     * \section Algorithm (uniform)
     * 1) G = n - k*L (totalGap), M = G + k.
     * 2) Sample k distinct bar positions B (sorted) from [0 .. M-1] uniformly.
     * 3) Convert to gaps:
     *      g_0 = B[0];
     *      g_i = B[i] - B[i-1] - 1  for i=1..k-1;
     *      g_k = (M-1) - B[k-1];
     *    (Note: sum g_i = G).
     * 4) Place clusters in increasing order:
     *      pos = g_0;
     *      for i in 0..k-1: cluster=[pos, pos+L), pos = pos+L + g_{i+1}.
     *
     * \section Complexity
     * - Time: O(k log k) due to sorting the k bar positions; building gaps/clusters is O(k).
     * - Space: O(k).
     */
    static std::vector<Cluster> layout_random_stars_bars(int n, int L, int k, std::mt19937& gen) {
        std::vector<Cluster> out;

        const long long totalLen = 1LL * k * L;
        assert(totalLen <= n && "k*L must fit into n (non-overlapping clusters)");

        const int G = n - static_cast<int>(totalLen); // total free gap
        const int M = G + k;                           // slots = stars + bars

        // Uniformly choose k distinct bar positions in [0..M-1]
        std::vector<int> bars = sample_k_unique_sorted(M, k, gen); // sorted

        // Convert bar positions to non-negative gaps summing to G
        std::vector<int> gaps; gaps.reserve(k + 1);
        if (k == 0) {
            gaps.push_back(G);
        } else {
            gaps.push_back(bars[0]); // g0
            for (int i = 1; i < k; ++i)
                gaps.push_back(bars[i] - bars[i-1] - 1);
            gaps.push_back((M - 1) - bars.back()); // g_k
        }

        // Build clusters in increasing order
        out.reserve(k);
        int pos = gaps[0];
        for (int i = 0; i < k; ++i) {
            const int s = pos;
            const int e = s + L;
            out.push_back({ s, e });
            pos = e + gaps[i + 1];
        }

        // checks
        for (int i = 1; i < k; ++i) {
            assert(out[i - 1].end <= out[i].start && "Clusters must not overlap");
        }
        assert(out.back().end <= n && "Last cluster must fit into array bounds");
        return out;
    }

    /**
     * \brief Helper: uniformly sample k distinct integers from [0..M-1] and return them sorted.
     * Uses Floyd's algorithm.
     */
    static std::vector<int> sample_k_unique_sorted(int M, int k, std::mt19937& gen) {
        std::unordered_set<int> chosen;
        chosen.reserve(static_cast<size_t>(k) * 2);

        std::uniform_int_distribution<int> dist; // will set range in-loop
        for (int j = M - k; j < M; ++j) {
            dist = std::uniform_int_distribution<int>(0, j);
            int t = dist(gen);
            if (!chosen.insert(t).second) {
                chosen.insert(j);
            }
        }
        std::vector<int> out;
        out.reserve(k);
        for (int v : chosen) out.push_back(v);
        std::sort(out.begin(), out.end());
        return out;
    }
    /**
     * \brief Sum of per-cluster right-side capacities to extend by +1 without overlap.
     *
     * \section Description
     * For each cluster, compute how many single-step rightward extensions are
     * possible before hitting the next cluster’s start (or the array end for
     * the last cluster). Returns the sum over all clusters.
     *
     * \section Algorithm
     * - For each cluster i, define limit_i as:
     *      - limit_i = clusters[i+1].start for i < k-1
     *      - limit_{k-1} = n
     * - Capacity_i = max(0, limit_i - clusters[i].end).
     * - Return sum_i Capacity_i.
     *
     * \section Notes
     * - Includes internal gaps and the right tail after the last cluster.
     *
     * \section Complexity
     * - Time: O(k), Space: O(1).
     *
     * @param cs Sorted non-overlapping clusters.
     * @param n  Array length.
     * @return Total number of "+1" extensions possible to the right.
     */
    static int right_extension_capacity(const std::vector<Cluster>& cs, int n) {
        const int k = static_cast<int>(cs.size());
        if (k == 0) return 0;
        int cap = n - cs.back().end;                    // right tail gap g_k
        for (int i = 0; i + 1 < k; ++i)                 // internal gaps g_1..g_{k-1}
            cap += cs[i + 1].start - cs[i].end;
        return cap;
    }

    /**
     * \brief Unbiased reminder distribution.
     *
     * \section Description
     * In each round, we collect all clusters that can grow by one index to the
     * right without crossing a neighbor. We then randomize their order and
     * extend as many distinct clusters by +1 as the remaining reminder allows.
     * Repeating this process eliminates positional bias when capacities differ.
     *
     * \section Algorithm
     * 1. While reminder > 0:
     *    1.1 Build a list of eligible cluster indices i where clusters[i].end < limit_i
     *        (limit_i = next.start for non-last, else n).
     *    1.2 Shuffle this list uniformly (std::shuffle with \p gen).
     *    1.3 Let take = min(reminder, eligible.size()).
     *    1.4 For the first \c take indices in the shuffled list: do clusters[i].end += 1.
     *    1.5 Decrease reminder by \c take and repeat.
     *
     * \section Notes
     * - “Random & unique per round” means no cluster gets two increments within
     *   the same round; which clusters are extended changes randomly each round.
     *
     * \section Complexity
     * - Time: O(R * k) worst-case, where R ~= ceil(reminder / avgEligiblePerRound).
     * - Space: O(k) for the temporary eligible list.
     *
     * @param reminder In/out: how many +1 increments still need to be allocated.
     * @param n        Array length (rightmost bound).
     * @param clusters In/out: clusters to extend.
     * @param gen      RNG used to shuffle per round.
     *
     * \section Preconditions / Throws
     * - Asserts (debug): reminder <= right_extension_capacity(clusters, n);
     *   at exit, reminder == 0. No exceptions thrown.
     */
    static void distribute_reminder_random(int& reminder,
                                           int n,
                                           std::vector<Cluster>& clusters,
                                           std::mt19937& gen) {
        if (reminder <= 0 || clusters.empty()) return;

        while (reminder > 0) {
            // Gather clusters that can grow by +1 to the right
            std::vector<int> growable;
            growable.reserve(clusters.size());
            for (int i = 0; i < static_cast<int>(clusters.size()); ++i) {
                const int limit = (i + 1 < static_cast<int>(clusters.size()))
                                  ? clusters[i + 1].start
                                  : n;
                if (clusters[i].end < limit) {
                    growable.push_back(i);
                }
            }

            if (growable.empty()) break; // safety for release builds

            // Random, no repeats per round
            std::shuffle(growable.begin(), growable.end(), gen);

            const int take = std::min(reminder, static_cast<int>(growable.size()));
            for (int j = 0; j < take; ++j) {
                clusters[growable[j]].end += 1;
            }
            reminder -= take;
        }

        assert(reminder == 0 && "Not enough right-side gap to place reminder safely");
    }


#ifdef CLUSTER_SAMPLER_DEBUG
    /**
     * \brief Debug print of final clusters and their lengths.

     */
    static void debug_print_clusters(const std::vector<Cluster>& clusters, int targetSample) {
        const int k = static_cast<int>(clusters.size());
        std::cerr << "[ClusterSampler] clusters.count = " << k << '\n';
        long long sumL = 0;
        for (int i = 0; i < k; ++i) {
            const int l = clusters[i].end - clusters[i].start;
            sumL += l;
            std::cerr << "  #" << i
                      << " start=" << clusters[i].start
                      << " end="   << clusters[i].end
                      << " len="   << l << '\n';
        }
        std::cerr << "[ClusterSampler] sum(len_i) = " << sumL
                  << " (targetSample=" << targetSample << ")\n";
    }
#endif

private:
    int clusterSize_;
};

#endif // CLUSTERSAMPLER_H