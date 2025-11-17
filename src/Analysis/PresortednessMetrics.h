
#ifndef DISORDERMETRICS_H
#define DISORDERMETRICS_H
#include <algorithm>
#include <cmath>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <vector>
/// \class PresortednessMetrics
/// \brief A collection of presortedness metrics (singleton).
///
/// Stateless utility that computes common presortedness metrics on integer arrays:
/// runs count, run-length entropy, distinct count, value-frequency entropy, Hartley entropy,
/// inversions, removals (n - LIS), oscillations, Spearman footrule sum, Hamming distance.
class PresortednessMetrics {
public:
    /// \brief Access the singleton instance.
    static PresortednessMetrics& instance() {
        static PresortednessMetrics inst;
        return inst;
    }

    PresortednessMetrics(const PresortednessMetrics&) = delete;
    PresortednessMetrics& operator=(const PresortednessMetrics&) = delete;
    PresortednessMetrics(PresortednessMetrics&&) = delete;
    PresortednessMetrics& operator=(PresortednessMetrics&&) = delete;

    /**
     * \brief Count strictly increasing runs.
     *
     * \section Description
     * The number of runs is the count of maximal contiguous increasing blocks
     * in the array. A new run starts at position i if a[i] < a[i-1].
     *
     * \section Algorithm
     * - If array is empty: return 0.
     * - Initialize runs = 1.
     * - For i = 1..n-1: if a[i] < a[i-1], ++runs.
     *
     * \section Complexity
     * - Time: O(n). Space: O(1).
     */
    long long calculateRuns(const std::vector<int>& arr) const {
        if (arr.empty()) return 0;
        long long runsCount = 1;
        for (size_t i = 1; i < arr.size(); ++i) {
            if (arr[i] < arr[i - 1]) ++runsCount;
        }
        return runsCount;
    }

    /**
     * \brief Shannon entropy of run-length distribution.
     *
     * \section Description
     * Let r_1, …, r_k be lengths of increasing runs with sum n.
     * Define probabilities p_i = r_i / n. The run-length entropy is
     * H = -∑ p_i log p_i (natural logarithm).
     *
     * \section Algorithm
     * - Scan array once to produce run lengths r_i.
     * - For each r_i, accumulate H -= (r_i/n) * log(r_i/n).
     *
     * \section Complexity
     * - Time: O(n). Space: O(k) for run lengths (k = number of runs).
     */
    double calculateRunsEntropy(const std::vector<int>& arr) const {
        const size_t n = arr.size();
        if (n == 0) return 0.0;

        std::vector<long long> runLengths;
        runLengths.reserve(n);
        long long curLen = 1;

        for (size_t i = 1; i < n; ++i) {
            if (arr[i] < arr[i - 1]) {
                runLengths.push_back(curLen);
                curLen = 1;
            } else {
                ++curLen;
            }
        }
        runLengths.push_back(curLen);

        long double H = 0.0L;
        const long double N = static_cast<long double>(n);
        for (long long r : runLengths) {
            long double p = static_cast<long double>(r) / N;
            H -= p * std::log(p);
        }
        return static_cast<double>(H);
    }

    /**
     * \brief Number of distinct values.
     *
     * \section Description
     * The count of unique elements present in the array.
     *
     * \section Algorithm
     * - Insert all elements into an unordered_set; return its size.
     *
     * \section Complexity
     * - Time: O(n) expected. Space: O(U) where U is number of distinct values.
     */
    long long calculateDistinctElements(const std::vector<int>& arr) const {
        std::unordered_set<int> uniq;
        uniq.reserve(arr.size());
        for (int x : arr) uniq.insert(x);
        return static_cast<long long>(uniq.size());
    }

    /**
     * \brief Shannon entropy of the empirical value-frequency distribution.
     *
     * \section Description
     * Let f(v) be the frequency of value v in the array, and n the array length.
     * Define p(v) = f(v)/n. The value-frequency entropy is
     * H = -∑_v p(v) log p(v) (natural logarithm).
     *
     * \section Algorithm
     * - Tally frequencies into a hash map.
     * - For each frequency f, accumulate H -= (f/n) * log(f/n).
     *
     * \section Complexity
     * - Time: O(n) expected. Space: O(U) where U is number of distinct values.
     */
    double calculateValueFrequencyEntropy(const std::vector<int>& arr) const {
        const long long n = static_cast<long long>(arr.size());
        if (n == 0) return 0.0;

        std::unordered_map<long long, long long> freq;
        freq.reserve(arr.size());
        for (long long x : arr) ++freq[x];

        const long double N = static_cast<long double>(n);
        long double H = 0.0L;
        for (const auto& kv : freq) {
            const long double p = static_cast<long double>(kv.second) / N;
            H -= p * std::log(p); // natural log
        }
        return static_cast<double>(H);
    }

    /**
     * \brief Hartley entropy (H0) of the set of observed values.
     *
     * \section Description
     * Let U be the number of distinct values. The Hartley entropy is
     * H0 = ln(U) (natural logarithm).
     *
     * \section Algorithm
     * - Compute U = number of distinct elements.
     * - If U <= 1 → return 0; else return ln(U).
     *
     * \section Complexity
     * - Time: O(n) expected. Space: O(U).
     */
    double calculateDistinctElementEntropy(const std::vector<int>& arr) const {
        const long long U = calculateDistinctElements(arr);
        if (U <= 1) return 0.0;
        return static_cast<double>(std::log(static_cast<long double>(U)));
    }

    /**
     * \brief Number of inversions.
     *
     * \section Description
     * An inversion is a pair (i,j) with i<j and a[i] > a[j]. The metric equals
     * the total number of such pairs.
     *
     * \section Algorithm
     * - Use mergesort-based counting:
     *   during merge of two sorted halves, when taking an element from the right half
     *   over the left, add (remaining in left) to the inversion count.
     *
     * \section Complexity
     * - Time: O(n log n). Space: O(n).
     */
    long long calculateInversions(const std::vector<int>& arr) const {
        if (arr.size() < 2) return 0;
        std::vector<int> tmp = arr;
        return mergeSortAndCollect(tmp, 0, static_cast<int>(tmp.size()) - 1);
    }

    /**
     * \brief Removals to sort: n − LIS length.
     *
     * \section Description
     * The minimal number of elements to remove to make the array non-decreasing,
     * equal to n − LIS, where LIS is the length of the longest non-decreasing subsequence
     * (here LIS computed in non-decreasing sense via upper_bound).
     *
     * \section Algorithm
     * - Patience sorting approach:
     *   maintain tails[]; for each x use upper_bound(tails, x) to place/replace;
     *   LIS = tails.size(); return n − LIS.
     *
     * \section Complexity
     * - Time: O(n log n). Space: O(n) in the worst case.
     *
     * \section Notes
     * - Using upper_bound makes the LIS non-decreasing (ties allowed).
     */
    long long calculateRem(const std::vector<int>& arr) const {
        if (arr.empty()) return 0;
        std::vector<int> tail;
        tail.reserve(arr.size());
        for (int x : arr) {
            auto it = std::upper_bound(tail.begin(), tail.end(), x);
            if (it == tail.end()) tail.push_back(x);
            else *it = x;
        }
        long long lis = static_cast<long long>(tail.size());
        return static_cast<long long>(arr.size()) - lis;
    }

    /**
     * \brief Oscillation count (peaks + valleys).
     *
     * \section Description
     * The number of interior indices i (1 <= i <= n−2) that are either a strict peak
     * (a[i-1] < a[i] > a[i+1]) or a strict valley (a[i-1] > a[i] < a[i+1]).
     *
     * \section Algorithm
     * - For i = 1..n-2:
     *   - if (a[i-1] < a[i] && a[i] > a[i+1]) ++count;
     *   - else if (a[i-1] > a[i] && a[i] < a[i+1]) ++count.
     *
     * \section Complexity
     * - Time: O(n). Space: O(1).
     */
    long long calculateOsc(const std::vector<int>& arr) const {
        if (arr.size() < 3) return 0;
        long long count = 0;
        for (size_t i = 1; i + 1 < arr.size(); ++i) {
            const bool isPeak   = (arr[i - 1] < arr[i] && arr[i] > arr[i + 1]);
            const bool isValley = (arr[i - 1] > arr[i] && arr[i] < arr[i + 1]);
            if (isPeak || isValley) ++count;
        }
        return count;
    }

    /**
     * \brief Spearman footrule distance to the stably sorted order.
     *
     * \section Description
     * Let π be the stable sort permutation: positions of each value if the array
     * were stably sorted non-decreasing. The footrule distance is
     * ∑_i |i − π(i)|.
     *
     * \section Algorithm
     * - Make a stably sorted copy of the array.
     * - Build for each value v a queue of its target positions in the sorted copy.
     * - Scan original array left-to-right; for current value v take and pop the leftmost
     *   target position t from v’s queue; add |i − t| to the sum.
     *
     * \section Complexity
     * - Time: O(n log n) due to sorting. Space: O(n).
     */
    long long calculateDis(const std::vector<int>& arr) const {
        std::vector<int> sorted = arr;
        std::stable_sort(sorted.begin(), sorted.end());

        std::unordered_map<int, std::deque<int>> pos;
        pos.reserve(sorted.size());
        for (int i = 0; i < static_cast<int>(sorted.size()); ++i) {
            pos[sorted[i]].push_back(i);
        }

        long long sum = 0;
        for (int i = 0; i < static_cast<int>(arr.size()); ++i) {
            int idx = pos[arr[i]].front();
            pos[arr[i]].pop_front();
            sum += std::llabs(static_cast<long long>(idx) - i);
        }
        return sum;
    }

    /**
     * \brief Hamming distance to the stably sorted array.
     *
     * \section Description
     * The number of positions i where the current value differs from the value at i
     * in the stably sorted non-decreasing version of the array.
     *
     * \section Algorithm
     * - Make a stably sorted copy s of the array.
     * - Count positions i with arr[i] != s[i].
     *
     * \section Complexity
     * - Time: O(n log n) due to sorting. Space: O(n).
     */
    long long calculateHam(const std::vector<int>& arr) const {
        std::vector<int> sorted = arr;
        std::stable_sort(sorted.begin(), sorted.end());
        long long count = 0;
        for (size_t i = 0; i < arr.size(); ++i) {
            if (arr[i] != sorted[i]) ++count;
        }
        return count;
    }


private:
    PresortednessMetrics() = default;

    /**
     * \brief Merge step that also counts cross-inversions.
     *
     * \section Description
     * Counts pairs (i,j) with i in left half and j in right half such that L[i] > R[j].
     *
     * \section Algorithm
     * - Standard merge; when taking R[j] over L[i], add (leftSize − i) to the count.
     *
     * \section Complexity
     * - Time: O(n). Space: O(n) for temporary halves.
     */
    long long mergeAndCollect(std::vector<int>& arr, int left, int mid, int right) const {
        const int leftSize  = mid - left + 1;
        const int rightSize = right - mid;

        std::vector<int> L(leftSize), R(rightSize);
        for (int i = 0; i < leftSize;  ++i) L[i] = arr[left + i];
        for (int j = 0; j < rightSize; ++j) R[j] = arr[mid + 1 + j];

        int i = 0, j = 0, k = left;
        long long invCount = 0;

        while (i < leftSize && j < rightSize) {
            if (L[i] <= R[j]) {
                arr[k++] = L[i++];
            } else {
                invCount += (leftSize - i);
                arr[k++] = R[j++];
            }
        }
        while (i < leftSize)  arr[k++] = L[i++];
        while (j < rightSize) arr[k++] = R[j++];
        return invCount;
    }

    /**
     * \brief Recursive mergesort that returns the inversion count.
     *
     * \section Algorithm
     * - Recurse on halves; sum their inversion counts; add cross-inversions from merge step.
     *
     * \section Complexity
     * - Time: O(n log n). Space: O(n).
     */
    long long mergeSortAndCollect(std::vector<int>& arr, int left, int right) const {
        long long invCount = 0;
        if (left < right) {
            const int mid = left + (right - left) / 2;
            invCount += mergeSortAndCollect(arr, left, mid);
            invCount += mergeSortAndCollect(arr, mid + 1, right);
            invCount += mergeAndCollect(arr, left, mid, right);
        }
        return invCount;
    }
};

#endif // DISORDER_METRICS_SINGLETON_H