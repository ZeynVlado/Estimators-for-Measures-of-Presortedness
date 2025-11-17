#include <iostream>

#include "Generation/RunGeneratorNoiseBased.h"
#include "Generation/RunGeneratorSegmentsBased.h"
#include "Sampling/ClusterSampler.h"
#include "Sampling/CombinedSampler.h"
#include "Sampling/StratifiedSampler.h"

#include "src/Analysis/Evaluator.h"
#include "src/Data/Generation/DataGenerator.h"
#include "src/Data/Sampling/Sampler.h"
#include "src/Analysis/PresortednessMetrics.h"
#include "src/Data/ExperimentConfigurator.h"

namespace fs = std::filesystem;
int main() {
    // const int n        = 10000;
    // const int minValue = 1;
    // const int maxValue = 10000;
    // const std::string root =
    //     "/src/experiment_data_input";
    //
    // const int countPerSet = 1000;
    //
    //
    // std::vector<int> kValues = {2500, 5000, 7500, 9000};
    //
    // std::vector<int> uValues = {2500, 5000, 7500, 9000};
    //
    // const int sqrt_n   = static_cast<int>(std::floor(std::sqrt(n)));
    // const int log2_n   = static_cast<int>(std::floor(std::log2(n)));
    // const int n_div_25 = static_cast<int>(std::lround(n / 2.5));
    //
    // std::vector<int> runsValues = {
    //     13, 100, 200, 1000, 2500,
    //     sqrt_n,
    //     2 * sqrt_n,
    //     log2_n,
    //     n / 4, n / 10, n / 20,
    //     n_div_25
    // };
    //
    //
    // try
    // {
    //     const int S = 100;
    //     ExperimentConfigurator conf(n, S, minValue, maxValue, root);
    //     conf.configure(countPerSet, kValues, uValues, runsValues,
    //                    /*includePermutation*/ true,
    //                    /*includeMultisetPermutation*/ true,
    //                    /*includeRandomK*/     true,
    //                    /*includeRuns*/        true,
    //                    /*overwrite*/          true);
    //     std::cout << "Done: sampleSize = " << S << "\n";
    // }
        //
        //    // S = 200
        //     // /*
    // {
    //     const int S = 200;
    //     ExperimentConfigurator conf(n, S, minValue, maxValue, root);
    //     conf.configure(countPerSet, kValues, uValues, runsValues,
    //                    /*includePermutation*/ true,
    //                    /*includeMultisetPermutation*/ true,
    //                    /*includeRandomK*/     true,
    //                    /*includeRuns*/        true,
    //                    /*overwrite*/          false);
    //     std::cout << "Done: sampleSize = " << S << "\n";
    // }
        //
        //     //S = 300
        //     // /*
    // {
    //     const int S = 300;
    //     ExperimentConfigurator conf(n, S, minValue, maxValue, root);
    //     conf.configure(countPerSet, kValues, uValues, runsValues,
    //                    /*includePermutation*/ true,
    //                    /*includeMultisetPermutation*/ true,
    //                    /*includeRandomK*/     true,
    //                    /*includeRuns*/        true,
    //                    /*overwrite*/          false);
    //     std::cout << "Done: sampleSize = " << S << "\n";
    // }
        //
        //    // S = 400
    // {
    //     const int S = 400;
    //     ExperimentConfigurator conf(n, S, minValue, maxValue, root);
    //     conf.configure(countPerSet, kValues, uValues, runsValues,
    //                    /*includePermutation*/ true,
    //                    /*includeMultisetPermutation*/ true,
    //                    /*includeRandomK*/     true,
    //                    /*includeRuns*/        true,
    //                    /*overwrite*/          false);
    //     std::cout << "Done: sampleSize = " << S << "\n";
    // }
        //
        // S = 500
    // {
    //     const int S = 500;
    //     ExperimentConfigurator conf(n, S, minValue, maxValue, root);
    //     conf.configure(countPerSet, kValues, uValues, runsValues,
    //                    /*includePermutation*/ true,
    //                    /*includeMultisetPermutation*/ true,
    //                    /*includeRandomK*/     true,
    //                    /*includeRuns*/        true,
    //                    /*overwrite*/          false);
    //     std::cout << "Done: sampleSize = " << S << "\n";
    // }
        //
        //    // S = 600
    // {
    //     const int S = 600;
    //     ExperimentConfigurator conf(n, S, minValue, maxValue, root);
    //     conf.configure(countPerSet, kValues, uValues, runsValues,
    //                    /*includePermutation*/ true,
    //                    /*includeMultisetPermutation*/ true,
    //                    /*includeRandomK*/     true,
    //                    /*includeRuns*/        true,
    //                    /*overwrite*/          false);
    //     std::cout << "Done: sampleSize = " << S << "\n";
    // }

        //
        // const std::string inputRoot  = "/src/experiment_data_input";
        // const std::string outputRoot = "/src/experiment_data_output";
        //
        // try {
        //     Evaluator ev;
        //     ev.prepareOutputStructure(inputRoot, outputRoot);
        //     ev.evaluateAll(inputRoot, outputRoot, true);
        //     std::cout << "Analysis completed. Results at: " << outputRoot << "\n";
        // } catch (const std::exception& ex) {
        //     std::cerr << "ERROR: " << ex.what() << "\n";
        //     return 1;
        // }






    return 0;
}











