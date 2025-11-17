#ifndef SAMPLERTEST_H
#define SAMPLERTEST_H

#define private public
#include "Sampling/Sampler.h"
#undef private


#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <vector>

#include <numeric>

#include "Sampling/ClusterSampler.h"
#include "Sampling/CombinedSampler.h"
#include "Sampling/StratifiedSampler.h"

namespace {
constexpr int N      = 10000;
constexpr int SAMPLE = 100;

std::vector<int> makeArray(int n) {
    std::vector<int> a(n);
    std::iota(a.begin(), a.end(), 0);
    return a;
}

bool strictlyIncreasing(const std::vector<int>& v) {
    for (size_t i = 1; i < v.size(); ++i) if (!(v[i] > v[i-1])) return false;
    return true;
}


} // namespace

// STRATIFIED

TEST(StratifiedSamplerTest, BasicPropertiesAndStrataCoverage) {
    auto arr = makeArray(N);

    std::vector<int> stratLens = {
        (int)std::sqrt((double)N), (int)(2*std::sqrt((double)N)), N/4, N/10, N/20
    };

    for (int stratLen : stratLens) {
        StratifiedSampler s(SAMPLE, stratLen);
        s.createSample(arr);

        const auto& idx = s.sampleIndices();
        const auto& smp = s.sample();

        ASSERT_EQ((int)idx.size(), SAMPLE);
        ASSERT_EQ((int)smp.size(), SAMPLE);

        ASSERT_TRUE(strictlyIncreasing(idx));

        ASSERT_EQ((int)arr.size(), N);
        for (size_t k = 0; k < idx.size(); ++k)
            EXPECT_EQ(smp[k], arr[idx[k]]);


        const int K = (N + stratLen - 1) / stratLen;
        std::vector<int> perStratum(K, 0);
        for (int id : idx) {
            int k = std::min(id / stratLen, K - 1);
            ++perStratum[k];
        }
        ASSERT_EQ((int)perStratum.size(), K);
        for (int k = 0; k < K; ++k)
            EXPECT_GE(perStratum[k], 1) << "stratum=" << k << " stratLen=" << stratLen;
    }
}

// CLUSTER

TEST(ClusterSamplerTest, RunsMatchBaseClusterCountAndBasics) {
    auto arr = makeArray(N);

    std::vector<int> clusterSizes = {
        (int)std::sqrt((double)SAMPLE),                 // ~10
        (int)(2*std::sqrt((double)SAMPLE)),             // ~20
        (int)std::log2((double)SAMPLE),                 // ~6
        (int)std::log2((double)SAMPLE * SAMPLE)         // ~13
    };

    for (int cs : clusterSizes) {
        ClusterSampler s(SAMPLE, cs);
        s.createSample(arr);

        const auto& idx = s.sampleIndices();
        const auto& smp = s.sample();

        ASSERT_EQ((int)idx.size(), SAMPLE);
        ASSERT_EQ((int)smp.size(), SAMPLE);
        ASSERT_TRUE(strictlyIncreasing(idx));
        ASSERT_EQ((int)arr.size(), N);
        for (size_t k = 0; k < idx.size(); ++k)
            EXPECT_EQ(smp[k], arr[idx[k]]);


    }
}

// COMBINED

TEST(CombinedSamplerTest, PerStratumAtLeastOneAndContiguousInsideStratum) {
    auto arr = makeArray(N);

    std::vector<int> stratLens = {
        (int)std::sqrt((double)N), (int)(2*std::sqrt((double)N)), N/4, N/10, N/20
    };

    for (int stratLen : stratLens) {
        CombinedSampler s(SAMPLE, stratLen);
        s.createSample(arr);

        const auto& idx = s.sampleIndices();
        const auto& smp = s.sample();

        ASSERT_EQ((int)idx.size(), SAMPLE);
        ASSERT_EQ((int)smp.size(), SAMPLE);
        ASSERT_TRUE(strictlyIncreasing(idx));
        ASSERT_EQ((int)arr.size(), N);
        for (size_t k = 0; k < idx.size(); ++k)
            EXPECT_EQ(smp[k], arr[idx[k]]);


        const int K = (N + stratLen - 1) / stratLen;
        std::vector<std::vector<int>> byStratum(K);
        for (int id : idx) {
            int k = std::min(id / stratLen, K - 1);
            byStratum[k].push_back(id);
        }

        for (int k = 0; k < K; ++k) {
            ASSERT_GE((int)byStratum[k].size(), 1) << "stratum=" << k << " stratLen=" << stratLen;

            const auto& v = byStratum[k];
            const int mn = v.front();
            const int mx = v.back();
            bool contiguous = (int)v.size() == (mx - mn + 1);
            EXPECT_TRUE(contiguous) << "stratum=" << k << " size=" << v.size()
                                    << " mn=" << mn << " mx=" << mx
                                    << " stratLen=" << stratLen;
        }
    }
}


#endif // SAMPLERTEST_H


