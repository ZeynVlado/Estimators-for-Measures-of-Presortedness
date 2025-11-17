#include <gtest/gtest.h>

#include "PresortednessMetrics.h"

#include <random>

#include <vector>
#include <stdexcept>

#include "Generation/RunArrGenerator.h"


#ifndef RUNGENERATORS_TEST_H
#define RUNGENERATORS_TEST_H

#include "Generation/RunGeneratorNoiseBased.h"
#include "Generation/RunGeneratorSegmentsBased.h"
#include "Generation/RndArrGenerator.h"
#include "Generation/RunGeneratorNoiseBased.h"


//
class DataGeneratorTest : public ::testing::Test {
protected:
    PresortednessMetrics& dm = PresortednessMetrics::instance();
};

// Tests for RunGeneratorNoiseBased
TEST_F(DataGeneratorTest, NoiseBased_GeneratesExactNumberOfRuns) {
    RunGeneratorNoiseBased gen;

    struct Case { size_t n; int runs; int minV; int maxV; };
    std::vector<Case> cases = {
        {10, 1, 0, 100},
        {10, 2, 0, 100},
        {25, 5, 0, 100},
        {50, 10, 0, 100},
        {50, 20, 0, 100},
    };

    for (const auto& c : cases) {
        SCOPED_TRACE(testing::Message()
                     << "n=" << c.n << " runs=" << c.runs
                     << " minV=" << c.minV << " maxV=" << c.maxV);
        auto a = gen.generateRunsFromSortedArrayNoise(c.n, c.runs, c.minV, c.maxV);
        ASSERT_EQ(a.size(), c.n);
        const long long r = dm.calculateRuns(a);
        EXPECT_EQ(r, c.runs);
    }
}

// Tests for RunGeneratorSegmentsBased
TEST_F(DataGeneratorTest, SegmentsBased_FromPermutation_GeneratesExactNumberOfRuns) {
    RunGeneratorSegmentsBased gen;

    std::vector<std::pair<size_t,size_t>> cases = {
        {1, 1},
        {8, 1},
        {8, 2},
        {16, 4},
        {32, 8},
        {64, 16},
    };

    for (auto [n, runs] : cases) {
        SCOPED_TRACE(testing::Message() << "n=" << n << " runs=" << runs);
        auto a = gen.generateRunsFromPermutationSortedSegments(n, runs);
        ASSERT_EQ(a.size(), n);
        const long long r = dm.calculateRuns(a);
        EXPECT_EQ(r, static_cast<long long>(runs));
    }
}

TEST_F(DataGeneratorTest, SegmentsBased_FromSorted_GeneratesExactNumberOfRuns) {
    RunGeneratorSegmentsBased gen;

    struct Case { size_t n; size_t runs; int minV; int maxV; };
    std::vector<Case> cases = {
        {5,  1, 0, 100},
        {5,  2, 0, 100},
        {12, 3, 0, 100},
        {30, 6, 0, 100},
    };

    for (const auto& c : cases) {
        SCOPED_TRACE(testing::Message() << "n=" << c.n << " runs=" << c.runs);
        auto a = gen.generateRunsFromSortedArrayPermutatedSegments(c.n, c.runs, c.minV, c.maxV);
        ASSERT_EQ(a.size(), c.n);
        const long long r = dm.calculateRuns(a);
        EXPECT_EQ(r, static_cast<long long>(c.runs));
    }
}


TEST_F(DataGeneratorTest, MixedGenerators_RandomizedSmoke_ExactRuns) {
    std::vector<size_t> sizes = {10, 15, 20, 30};
    for (size_t n : sizes) {

        std::vector<int> targets = {1, 2, static_cast<int>(n/5), static_cast<int>(n/3)};
        for (int R : targets) {
            if (R <= 0 || static_cast<size_t>(R) > n) continue;

            // NoiseBased
            {
                RunGeneratorNoiseBased g;
                auto a = g.generateRunsFromSortedArrayNoise(
                    n, R, /*minValue*/0, /*maxValue*/100);
                EXPECT_EQ(dm.calculateRuns(a), R)
                    << "NoiseBased failed for n="<<n<<" R="<<R;
            }

            // SegmentsBased from permutation
            {
                RunGeneratorSegmentsBased g;
                auto a = g.generateRunsFromPermutationSortedSegments(n, static_cast<size_t>(R));
                EXPECT_EQ(dm.calculateRuns(a), R)
                    << "Segments/Perm failed for n="<<n<<" R="<<R;
            }

            // SegmentsBased from sorted base
            {
                RunGeneratorSegmentsBased g;
                auto a = g.generateRunsFromSortedArrayPermutatedSegments(
                    n, static_cast<size_t>(R), 0, 100);
                EXPECT_EQ(dm.calculateRuns(a), R)
                    << "Segments/Sorted failed for n="<<n<<" R="<<R;
            }
        }
    }
}

// Tests for RndArrGenerator::generateMultisetPermutation
TEST_F(DataGeneratorTest, MultisetPermutation_HasExactUniverseOfDistinctValues) {
    RndArrGenerator gen;
    struct Case { size_t total; int universe; int base; };
    std::vector<Case> cases = {
        {1, 1, 5},
        {5, 2, 1},
        {10, 3, 10},
        {20, 5, -7},
        {50, 7, 0},
    };

    for (const auto& c : cases) {
        SCOPED_TRACE(testing::Message()
                     << "total=" << c.total
                     << " universe=" << c.universe
                     << " base=" << c.base);

        auto a = gen.generateMultisetPermutation(c.total, c.universe, c.base);

        // Length matches
        ASSERT_EQ(a.size(), c.total);

        const int lo = c.base;
        const int hi = c.base + c.universe - 1;
        for (int v : a) {
            EXPECT_LE(lo, v);
            EXPECT_LE(v, hi);
        }

        std::vector<int> freq(c.universe, 0);
        for (int v : a) ++freq[v - c.base];

        int nonZero = 0;
        for (int cnt : freq) {
            EXPECT_GT(cnt, 0) << "Each label must appear at least once";
            if (cnt > 0) ++nonZero;
        }
        EXPECT_EQ(nonZero, c.universe);

        EXPECT_EQ(dm.calculateDistinctElements(a), c.universe);
    }
}

#endif // RUNGENERATORS_TEST_H








