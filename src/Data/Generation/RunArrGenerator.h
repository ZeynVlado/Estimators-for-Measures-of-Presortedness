
#ifndef RUNARRGENERATOR_H
#define RUNARRGENERATOR_H
#include <random>
#include <vector>

#include "DataGenerator.h"

class RunArrGenerator: public DataGenerator {

    using DataGenerator::DataGenerator;


protected:
    /**
    * @brief Build positive segment lengths whose sum equals @p size (no modulo-bias).
    *
    * Strategy:
    *  - Give 1 to each of @p runsCount segments (ensures positiveness).
    *  - Distribute the remaining @c size - runsCount units uniformly at random
    *    across segments using a uniform integer distribution.
    *
    * Guarantees:
    *  - All lengths > 0.
    *  - Sum(lengths) == size.
    *  - runsCount is clamped to [1..size].
    */
    std::vector<size_t> buildRunLengths(size_t size,
                                        size_t runsCount,
                                        std::mt19937& gen) {

        // Give 1 to each segment so all lengths are positive.
        std::vector<size_t> runLengths(runsCount, 1);

        // Distribute the remaining units uniformly at random across segments.
        size_t remaining = size - runsCount;

        if (remaining > 0) {
            std::uniform_int_distribution<size_t> pick(0, runsCount - 1);
            while (remaining--) {
                size_t idx = pick(gen);
                ++runLengths[idx];
            }
        }

        return runLengths;
    }

     /**
     * @brief Compute start indices for each segment as prefix sums of lengths.
     *        start[i] = sum_{j < i} length[j].
     */
    std::vector<size_t> buildRunStarts(const std::vector<size_t>& runLengths) {
        std::vector<size_t> runStarts(runLengths.size(), 0);
        size_t cursor = 0;
        for (size_t i = 0; i < runLengths.size(); ++i) {
            runStarts[i] = cursor;
            cursor += runLengths[i];
        }
        return runStarts;
    }

    /**
    * @brief Generate a sorted array of unique integers 1..size (1-based).
    * Complexity:
    *  - Time:   O(size)
    *  - Space:  O(1) beyond the output buffer
    *
    * @param size  Number of elements to produce.
    * @return std::vector<int>  The sequence {1, 2, ..., size}.
    */
    std::vector<int> generateSortedUnique(size_t size) {
        std::vector<int> data(size);
        // Fill with a strictly increasing 1-based sequence.
        std::iota(data.begin(), data.end(), 1);

        // Already sorted and unique by construction.
        return data;
    }

};

#endif //RUNARRGENERATOR_H
