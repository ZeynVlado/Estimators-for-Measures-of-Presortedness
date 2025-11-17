
#ifndef SAMPLER_H
#define SAMPLER_H

#include <vector>
#include <random>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "SamplingStrategy.h"

/**
 * \class Sampler
 * \brief Abstract base for sampling strategies over a 1D array.
 *
 * \section Description
 * Encapsulates the minimal interface and storage for sampling methods:
 * a target sample length, output buffers for sampled values and their
 * source indices, and a single virtual method \ref createSample that
 * materializes the sample given an input array. Concrete subclasses
 * define the sampling policy (e.g., stratified, clustered).
 */
class Sampler {
public:
    explicit Sampler(int sampleLength)
        : sampleLength_(sampleLength) {}

    virtual ~Sampler() = default;

    /// Build a sample from \p arr according to the concrete policy.
    virtual void createSample(const std::vector<int>& arr) = 0;

    /// Sampled values, in the order produced by the policy.
    const std::vector<int>& sample() const        { return sample_; }
    /// Source indices corresponding to \ref sample().
    const std::vector<int>& sampleIndices() const { return sampleIndices_; }

    /// Requested sample length (may be clamped by subclasses to array size).
    int  sampleLength() const { return sampleLength_; }
    /// Set requested sample length (no validation here; subclasses enforce feasibility).
    void setSampleLength(int v) { sampleLength_ = v; }

protected:
    std::vector<int> sample_;         ///< Output sample values.
    std::vector<int> sampleIndices_;  ///< Indices of chosen elements in the input array.
    int sampleLength_;                ///< Target sample length requested by the user.
};

#endif // SAMPLER_H
