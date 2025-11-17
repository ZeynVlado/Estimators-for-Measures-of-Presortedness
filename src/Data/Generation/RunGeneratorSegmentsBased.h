

#ifndef RUNGENERATORSEGMENTSBASED_H
#define RUNGENERATORSEGMENTSBASED_H
//#define RUNGENERATORSEGMENTSBASED_DEBUG
#include <assert.h>

#include "PresortednessMetrics.h"
#include "RndArrGenerator.h"
#include "RunArrGenerator.h"


/**
 * \class RunGeneratorSegmentsBased
 * \brief Segment-based run synthesis utility.
 *
 * \section Description
 *
 * Provides utilities to synthesize integer arrays that contain exactly a requested
 * number of strictly ascending runs. The class exposes two public pipelines:
 *  - start from a uniformly random permutation of {1..N};
 *  - start from a strictly increasing base [1..N].
 *
 * Both pipelines:
 *  - partition the index range into \c runsCount positive-length segments,
 *  - enforce local ascending order inside each segment,
 *  - reorder segments by strictly descending segment maximum,
 *  - concatenate segments to produce the final array.
 *
 */
class RunGeneratorSegmentsBased : public RunArrGenerator {
public:
    using RunArrGenerator::RunArrGenerator; // inherit default ctor (rng seeded via random_device)

    /**
     * \brief Main pipeline: exactly \p runsCount ascending runs starting from a random permutation.
     *
     * \section Description
     * Builds a random permutation of 1..size, splits it into \p runsCount segments
     * (positive lengths summing to size), sorts each segment ascending, reorders
     * segments by strictly descending segment maximum, and concatenates them.
     * This produces exactly \p runsCount strictly ascending runs.
     *
     * \section Algorithm
     * 1. perm <- uniformly random permutation of {1..size} using the member RNG.
     * 2. (runLengths, runStarts) <- positive segment lengths that sum to size,
     *    and their prefix-starts.
     * 3. Sort each segment ascending in-place.
     * 4. Compute per-segment metadata (start, len, min, max).
     * 5. Reorder segments by strictly descending max to force a boundary drop.
     * 6. Concatenate segments in the new order .
     *
     * \section Complexity
     * - Time: O(size log(size / runsCount)) amortized for local sorts + O(size) concatenation.
     * - Space: O(runsCount) for metadata + O(size) for output.
     *
     */
    std::vector<int> generateRunsFromPermutationSortedSegments(size_t size,
                                                               size_t runsCount)
    {
        if (size == 0) return {};

        assert(static_cast<int>(size) >= static_cast<int>(runsCount));
        assert(runsCount > 0);

        // Random permutation of 1..size
        std::vector<int> perm = generatePermutation(size);

        // Positive lengths summing to size (no modulo bias).
        std::vector<size_t> runLengths = buildRunLengths(size, runsCount, this->rng);
        std::vector<size_t> runStarts  = buildRunStarts(runLengths);

        assert(runLengths.size() == runsCount);
        assert(runStarts.size()  == runsCount);

        // Sum of lengths equals size
        const size_t totalLen = std::accumulate(runLengths.begin(), runLengths.end(), size_t{0});
        assert(totalLen == size && "Sum(runLengths) must equal size");

        // Sort each segment in strictly ascending order.
        sortEachSegment(perm, runStarts, runLengths);

        // Compute metadata.
        std::vector<Segment> segments = computeSegmentsMetaData(perm, runStarts, runLengths);

        // Order segments by strictly descending max.
        orderSegmentsByMaxDesc(segments);

        // Concatenate.
        std::vector<int> out = concatBySegments(perm, segments);

        assert(out.size() == perm.size());
        assert(PresortednessMetrics::instance().calculateRuns(out)  == static_cast<long long>(runsCount));

#ifdef RUNGENERATORSEGMENTSBASED_DEBUG
        debugDumpRuns(out, "[PermutationSortedSegments]");
#endif
        return out;
    }

    /**
     * \brief Alternative pipeline: exactly \p runsCount runs starting from a strictly increasing base.
     *
     * \section Description
     * Uses [1..size] as the base array, then applies the same segmentation,
     * metadata building, and segment reordering procedure to obtain exactly
     * \p runsCount strictly ascending runs in the result.
     *
     * \section Algorithm
     * 1) base <- [1..size].
     * 2) Build lengths and starts as above.
     * 3) Compute segment metadata; reorder segments by max descending.
     * 4) Concatenate segments → output.
     *
     * \section Complexity
     * - Time: O(size) (local sorts are trivial on a strictly increasing base)
     *          + O(size) concatenation.
     * - Space: O(runsCount) for metadata + O(size) for output.
     *
     */
    std::vector<int> generateRunsFromSortedArrayPermutatedSegments(size_t size,
                                                                   size_t runsCount,
                                                                   int /*minValue*/,
                                                                   int /*maxValue*/)
    {
        assert(static_cast<int>(size) >= static_cast<int>(runsCount));
        assert(runsCount > 0);

        std::vector<int> base = generateSortedUnique(size);

        std::vector<size_t> runLengths = buildRunLengths(size, runsCount, this->rng);
        std::vector<size_t> runStarts  = buildRunStarts(runLengths);

        assert(runLengths.size() == runsCount);
        assert(runStarts.size()  == runsCount);

        const size_t totalLen = std::accumulate(runLengths.begin(), runLengths.end(), size_t{0});
        assert(totalLen == size && "Sum(runLengths) must equal size");

        std::vector<Segment> segs = computeSegmentsMetaData(base, runStarts, runLengths);
        orderSegmentsByMaxDesc(segs);

        std::vector<int> out = concatBySegments(base, segs);

       ;
        assert(out.size() == base.size());
        assert( PresortednessMetrics::instance().calculateRuns(out) == static_cast<long long>(runsCount));

#ifdef RUNGENERATORSEGMENTSBASED_DEBUG
        debugDumpRuns(out, "[SortedArrayPermutatedSegments]");
#endif
        return out;
    }

private:
    struct Segment {
        size_t start;
        size_t len;
        int    minv;
        int    maxv;
    };

    /**
     * \brief Sort each segment ascending in-place.
     *
     * \section Description
     * For each (start, length) pair, sorts the slice [start, start+length) of \p a
     * using std::sort. Assumes segments are disjoint and lie within array bounds.
     *
     * \section Algorithm
     * - For i = 0..(#segments-1):
     *   - s <- runStarts[i], L <- runLengths[i]; if L == 0 skip (defensive).
     *   - sort(a.begin()+s, a.begin()+s+L) with std::sort.
     *
     * \section Complexity
     * - Time: sum over segments of O(L_i log L_i); amortized ~ O(size log(size / runsCount)).
     * - Space: O(1) extra (besides call stack).
     */
    void sortEachSegment(std::vector<int>& a,
                         const std::vector<size_t>& runStarts,
                         const std::vector<size_t>& runLengths)
    {
        const size_t n = runStarts.size();
        for (size_t i = 0; i < n; ++i) {
            const size_t L = runLengths[i];
            const size_t s = runStarts[i];
            if (L == 0) continue;
            const auto ds = static_cast<std::ptrdiff_t>(s);
            const auto de = static_cast<std::ptrdiff_t>(s + L);
            std::sort(a.begin() + ds, a.begin() + de);
        }
    }


    /**
     * \brief Compute per-segment metadata assuming each segment is already sorted ascending.
     *
     * \section Description
     * For segment i = [start, start+len), compute:
     *  - min = a[start],
     *  - max = a[start+len-1],
     *  - and record (start, len, min, max).
     *
     * \section Algorithm
     * - For i over segments:
     *   - s <- runStarts[i], L <- runLengths[i]; if L == 0 continue.
     *   - minv <- a[s]; maxv <- a[s+L-1];
     *   - push_back {s, L, minv, maxv}.
     *
     * \section Complexity
     * - Time: O(#segments).
     * - Space: O(#segments) for the returned vector.
     *
     */
    std::vector<Segment> computeSegmentsMetaData(const std::vector<int>& a,
                                                 const std::vector<size_t>& runStarts,
                                                 const std::vector<size_t>& runLengths)
    {
        std::vector<Segment> segments;
        segments.reserve(runStarts.size());
        for (size_t i = 0; i < runStarts.size(); ++i) {
            const size_t L = runLengths[i];
            const size_t s = runStarts[i];

            assert(L > 0 && "Each segment length must be positive");
            assert(s <= a.size() && "Start must be within [0..n]");
            assert(s + L <= a.size() && "End must not exceed n");

            if (L == 0) continue;
            const int minv = a[s];
            const int maxv = a[s + L - 1];
            segments.push_back(Segment{ s, L, minv, maxv });
        }
        return segments;
    }

    /**
     * \brief Reorder segments by strictly descending max value.
     *
     * \section Description
     * Sorts the segment descriptors so that segment maxima strictly decrease.
     * This guarantees a strict drop at every boundary after concatenation,
     * producing exactly one run boundary per adjacent pair.
     *
     * \section Complexity
     * - Time: O(#segments * log #segments).
     * - Space: O(1) extra.
     */
    void orderSegmentsByMaxDesc(std::vector<Segment>& segs)
    {
        std::sort(segs.begin(), segs.end(),
                  [](const Segment& a, const Segment& b) {
                      return a.maxv > b.maxv;
                  });
    }

    /**
     * \brief Concatenate segments in their current order.
     *
     * \section Description
     * Appends slices [start, start+len) from \p src in the order defined by \p segs.
     * The resulting array has total length equal to the sum of segment lengths.
     *
     * \section Algorithm
     * - out.reserve(sum(len_i)).
     * - For each seg in segs: out.insert(out.end(), src.begin()+start, src.begin()+start+len).
     *
     * \section Complexity
     * - Time: O(size) (single pass copy).
     * - Space: O(size) for the output vector.
     *
     */
    std::vector<int> concatBySegments(const std::vector<int>& src,
                                      const std::vector<Segment>& segs)
    {
        size_t total = 0;
        for (const auto& s : segs) total += s.len;

        std::vector<int> out;
        out.reserve(total);

        for (const auto& seg : segs) {
            const auto ds = static_cast<std::ptrdiff_t>(seg.start);
            const auto de = static_cast<std::ptrdiff_t>(seg.start + seg.len);
            out.insert(out.end(), src.begin() + ds, src.begin() + de);
        }
        return out;
    }

    // Small local helper: 1..size permutation using this->rng
    std::vector<int> generatePermutation(size_t size) {
        if (size == 0) return {};
        std::vector<int> v(size);
        std::iota(v.begin(), v.end(), 1);
        std::shuffle(v.begin(), v.end(), rng);
        return v;
    }

#ifdef RUNGENERATORSEGMENTSBASED_DEBUG
    /**
     * \brief Debug dump of all runs: indices, lengths, and elements.
     */
    static void debugDumpRuns(const std::vector<int>& a, const char* tag) {
        std::cerr << "\n[RunGeneratorSegmentsBased::debugDumpRuns] " << tag << "\n";
        const size_t n = a.size();
        if (n == 0) {
            std::cerr << "  (empty array)\n";
            return;
        }

        std::vector<std::pair<size_t,size_t>> runs;
        size_t start = 0;
        for (size_t i = 1; i < n; ++i) {
            if (a[i] < a[i-1]) {
                runs.emplace_back(start, i);
                start = i;
            }
        }
        runs.emplace_back(start, n);

        for (size_t r = 0; r < runs.size(); ++r) {
            const auto [s,e] = runs[r];
            std::cerr << "  run #" << r
                      << " start=" << s
                      << " end="   << e
                      << " len="   << (e - s)
                      << "  elements:";
            for (size_t i = s; i < e; ++i) std::cerr << ' ' << a[i];
            std::cerr << "\n";
        }

        const long long runsByMetric = DisorderMetrics::instance().calculateRuns(a) ;
        std::cerr << "  final length = " << n
                  << " | runs (scan) = " << runs.size()
                  << " | runs (DisorderMetrics) = " << runsByMetric
                  << "\n\n";
    }
#endif
};

#endif // RUNGENERATORSEGMENTSBASED_H