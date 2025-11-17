

#ifndef SAMPLERFACTORY_H
#define SAMPLERFACTORY_H
#include "ClusterSampler.h"
#include "CombinedSampler.h"
#include "StratifiedSampler.h"

/**
 * Factory that returns value objects of concrete samplers.
 */
class SamplerFactory {
public:
    static StratifiedSampler makeStratified(int sampleLength, int stratLength) {
        return StratifiedSampler(sampleLength, stratLength);
    }
    static ClusterSampler makeCluster(int sampleLength, int clusterSize) {
        return ClusterSampler(sampleLength, clusterSize);
    }
    static CombinedSampler makeCombined(int sampleLength, int stratSize) {
        return CombinedSampler(sampleLength, stratSize);
    }
};

#endif // SAMPLERFACTORY_H
