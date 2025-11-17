//
// EstimatorsForMeasuresOfPresortednessTest.h
//

#ifndef PresortednessMeasuresTEST_H
#define PresortednessMeasuresTEST_H

//#include ".../EstimatorsForMeasuresOfPresortedness/cmake-build-debug/_deps/googletest-src/googletest/include/gtest/gtest.h"

#include <vector>
#include <algorithm>
#include <cmath>

#include "PresortednessMetrics.h"


class PresortednessMeasures : public ::testing::Test { };

// BASIC EDGE CASES

TEST_F(PresortednessMeasures, EmptyAndSingle) {
    auto& dm = PresortednessMetrics::instance();

    std::vector<int> empty;
    EXPECT_EQ(dm.calculateRuns(empty), 0);
    EXPECT_DOUBLE_EQ(dm.calculateRunsEntropy(empty), 0.0);
    EXPECT_EQ(dm.calculateDistinctElements(empty), 0);
    EXPECT_DOUBLE_EQ(dm.calculateValueFrequencyEntropy(empty), 0.0);
    EXPECT_DOUBLE_EQ(dm.calculateDistinctElementEntropy(empty), 0.0);
    EXPECT_EQ(dm.calculateInversions(empty), 0);
    EXPECT_EQ(dm.calculateRem(empty), 0);
    EXPECT_EQ(dm.calculateOsc(empty), 0);
    EXPECT_EQ(dm.calculateDis(empty), 0);
    EXPECT_EQ(dm.calculateHam(empty), 0);

    std::vector<int> one{5};
    EXPECT_EQ(dm.calculateRuns(one), 1);
    EXPECT_NEAR(dm.calculateRunsEntropy(one), 0.0, 1e-15);
    EXPECT_EQ(dm.calculateDistinctElements(one), 1);
    EXPECT_NEAR(dm.calculateValueFrequencyEntropy(one), 0.0, 1e-15);
    EXPECT_NEAR(dm.calculateDistinctElementEntropy(one), 0.0, 1e-15);
    EXPECT_EQ(dm.calculateInversions(one), 0);
    EXPECT_EQ(dm.calculateRem(one), 0);
    EXPECT_EQ(dm.calculateOsc(one), 0);
    EXPECT_EQ(dm.calculateDis(one), 0);
    EXPECT_EQ(dm.calculateHam(one), 0);
}

// RUNS

TEST_F(PresortednessMeasures, NumberOfRunsCoreCases) {
    auto& dm = PresortednessMetrics::instance();

    // Strictly increasing -> 1 run
    std::vector<int> inc{1,2,3,4};
    EXPECT_EQ(dm.calculateRuns(inc), 1);

    // Non-decreasing with duplicates -> 1 run
    std::vector<int> nondec{1,1,2,2,3,3};
    EXPECT_EQ(dm.calculateRuns(nondec), 1);

    // Strictly decreasing -> n runs
    std::vector<int> dec{5,4,3,2,1};
    EXPECT_EQ(dm.calculateRuns(dec), 5);

    // Mixed with strict drops at 3->2 and 5->1 -> 3 runs
    std::vector<int> mix{1,3,2,2,5,1,1,4};
    EXPECT_EQ(dm.calculateRuns(mix), 3);

    // All equal -> 1 run
    std::vector<int> equal{7,7,7,7};
    EXPECT_EQ(dm.calculateRuns(equal), 1);
}

// SORTED VS REVERSE

TEST_F(PresortednessMeasures, SortedAndReverse) {
    auto& dm = PresortednessMetrics::instance();

    // Strictly increasing permutation
    std::vector<int> inc{1,2,3,4};
    // Runs: one block
    EXPECT_EQ(dm.calculateRuns(inc), 1);
    // Run entropy: one run of length 4 -> H = 0
    EXPECT_NEAR(dm.calculateRunsEntropy(inc), 0.0, 1e-12);
    // Inversions: already sorted -> 0
    EXPECT_EQ(dm.calculateInversions(inc), 0);
    // Removals: LIS length = 4 -> 0
    EXPECT_EQ(dm.calculateRem(inc), 0);
    // Oscillations: monotone -> 0
    EXPECT_EQ(dm.calculateOsc(inc), 0);
    // Dis/Ham/Max: all at correct positions
    EXPECT_EQ(dm.calculateDis(inc), 0);
    EXPECT_EQ(dm.calculateHam(inc), 0);
    // Distinct / entropies
    EXPECT_EQ(dm.calculateDistinctElements(inc), 4);
    EXPECT_NEAR(dm.calculateValueFrequencyEntropy(inc), std::log(4.0), 1e-12);
    EXPECT_NEAR(dm.calculateDistinctElementEntropy(inc), std::log(4.0), 1e-12);

    // Strictly decreasing permutation
    std::vector<int> dec{4,3,2,1};
    // Runs: each element starts a new run
    EXPECT_EQ(dm.calculateRuns(dec), 4);
    // Run entropy: 4 runs of length 1 -> H = ln 4
    EXPECT_NEAR(dm.calculateRunsEntropy(dec), std::log(4.0), 1e-12);
    // Inversions: C(4,2) = 6
    EXPECT_EQ(dm.calculateInversions(dec), 6);
    // LIS length = 1 -> removals = 3
    EXPECT_EQ(dm.calculateRem(dec), 3);
    // Oscillations: monotone decreasing -> 0
    EXPECT_EQ(dm.calculateOsc(dec), 0);
    // Distinct / entropies
    EXPECT_EQ(dm.calculateDistinctElements(dec), 4);
    EXPECT_NEAR(dm.calculateValueFrequencyEntropy(dec), std::log(4.0), 1e-12);
    EXPECT_NEAR(dm.calculateDistinctElementEntropy(dec), std::log(4.0), 1e-12);
    // Displacement / Hamming / Max displacement:
    // mapping to sorted [1,2,3,4] gives dis = 8, ham = 4, max = 3
    EXPECT_EQ(dm.calculateDis(dec), 8);
    EXPECT_EQ(dm.calculateHam(dec), 4);
}

// DUPLICATES & STABILITY

TEST_F(PresortednessMeasures, DuplicatesAndStability) {
    auto& dm = PresortednessMetrics::instance();

    // Array with two values, each repeated twice
    std::vector<int> a{2,2,1,1};
    // Runs: [2,2] and [1,1] -> 2 runs
    EXPECT_EQ(dm.calculateRuns(a), 2);
    // Inversions: (2,1) pairs -> 4
    EXPECT_EQ(dm.calculateInversions(a), 4);
    // LIS non-decreasing: length 2 -> removals = 2
    EXPECT_EQ(dm.calculateRem(a), 2);
    // Oscillations: no sign changes -> 0
    EXPECT_EQ(dm.calculateOsc(a), 0);

    // Displacement / Hamming / Max displacement:
    // sorted = [1,1,2,2]:
    // positions (0->2,1->3,2->0,3->1) -> dis = 8, ham = 4, max = 2
    EXPECT_EQ(dm.calculateDis(a), 8);
    EXPECT_EQ(dm.calculateHam(a), 4);

    // Value-frequency entropy: two values with p=0.5 each -> ln 2
    EXPECT_NEAR(dm.calculateValueFrequencyEntropy(a), std::log(2.0), 1e-12);
    // Hartley entropy: U=2 -> ln 2
    EXPECT_NEAR(dm.calculateDistinctElementEntropy(a), std::log(2.0), 1e-12);
    // Run entropy: run lengths 2,2 (n=4):
    // p=0.5,0.5 -> H = ln 2
    EXPECT_NEAR(dm.calculateRunsEntropy(a), std::log(2.0), 1e-12);
}

// OSCILLATIONS

TEST_F(PresortednessMeasures, OscillationPeaksValleys) {
    auto& dm = PresortednessMetrics::instance();

    // 1,3,2,4,3 has pattern peak-valley-peak: 3 oscillations
    std::vector<int> b{1,3,2,4,3};
    EXPECT_EQ(dm.calculateOsc(b), 3);

    // Flat array: no direction changes
    std::vector<int> flat{2,2,2,2};
    EXPECT_EQ(dm.calculateOsc(flat), 0);

    // Perfect zigzag: 1,3,1,3,1,3 -> 4 interior direction changes
    std::vector<int> zigzag{1,3,1,3,1,3};
    EXPECT_EQ(dm.calculateOsc(zigzag), 4);
}

// RUN ENTROPY EXAMPLES

TEST_F(PresortednessMeasures, RunsEntropyExamples) {
    auto& dm = PresortednessMetrics::instance();

    // b has run lengths 2,2,1 (n=5)
    std::vector<int> b{1,3,2,4,3};
    {
        const double n = 5.0;
        const double p1 = 2.0 / n;
        const double p2 = 2.0 / n;
        const double p3 = 1.0 / n;
        const double H =
            - (p1*std::log(p1) + p2*std::log(p2) + p3*std::log(p3));
        EXPECT_NEAR(dm.calculateRunsEntropy(b), H, 1e-12);
    }

    // Single full run -> H = 0
    std::vector<int> c{1,2,3,3,4,5};
    EXPECT_NEAR(dm.calculateRunsEntropy(c), 0.0, 1e-15);
}

// VALUE-FREQ & HARTLEY ENTROPY

TEST_F(PresortednessMeasures, ValueAndHartleyEntropy) {
    auto& dm = PresortednessMetrics::instance();

    // c: values 1,1,2,3 -> frequencies 2,1,1
    std::vector<int> c{1,1,2,3};
    {
        const double n = 4.0;
        const double p1 = 2.0 / n;
        const double p2 = 1.0 / n;
        const double p3 = 1.0 / n;
        const double H =
            - (p1*std::log(p1) + p2*std::log(p2) + p3*std::log(p3));
        EXPECT_NEAR(dm.calculateValueFrequencyEntropy(c), H, 1e-12);
    }
    // U = 3 -> Hartley entropy ln 3
    EXPECT_NEAR(dm.calculateDistinctElementEntropy(c), std::log(3.0), 1e-12);

    // All same value
    std::vector<int> allSame{7,7,7};
    EXPECT_NEAR(dm.calculateValueFrequencyEntropy(allSame), 0.0, 1e-15);
    EXPECT_NEAR(dm.calculateDistinctElementEntropy(allSame), 0.0, 1e-15);
    EXPECT_EQ(dm.calculateDistinctElements(allSame), 1);
}

// DIS / HAM /

TEST_F(PresortednessMeasures, DisHamMaxExamples) {
    auto& dm = PresortednessMetrics::instance();

    // Already sorted with duplicates -> all displacements zero
    std::vector<int> s{1,1,2,2,3};
    EXPECT_EQ(dm.calculateDis(s), 0);
    EXPECT_EQ(dm.calculateHam(s), 0);

    // Swap extremes to create nontrivial displacement
    std::vector<int> t{3,1,2,2,1};
    // Sorted stable: [1 (idx1),1 (idx4),2 (idx2),2(idx3),3(idx0)]
    // Distances: 0->4,1->0,2->2,3->3,4->1  (sum = 10, max = 4, ham = 4)
    EXPECT_EQ(dm.calculateDis(t), 8);
    EXPECT_EQ(dm.calculateHam(t), 2);
}

// REM

TEST_F(PresortednessMeasures, RemMatchesKnownLIS) {
    auto& dm = PresortednessMetrics::instance();

    // a: LIS non-decreasing has length 4 (e.g., 1,2,2,4) -> removals = 2
    std::vector<int> a{3,1,2,2,4,3};
    EXPECT_EQ(dm.calculateRem(a), 2);

    // b: best non-decreasing subsequence of length 3 (e.g., 1,2,2) -> removals = 3
    std::vector<int> b{4,4,4,1,2,2};
    EXPECT_EQ(dm.calculateRem(b), 3);
}

#endif // PresortednessMeasuresTEST_H
