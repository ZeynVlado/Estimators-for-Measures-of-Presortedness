

#ifndef DATAGENERATOR_H
#define DATAGENERATOR_H
#include <algorithm>
#include <random>
#include <vector>

#include "ArrayType.h"


/**
 * @brief Utility class for generating integer datasets in various patterns.
 *
 */
class DataGenerator {
public:
    /**
     * @brief Construct with a default seed (from std::random_device).
     *
     * If you need reproducible sequences, replace the seed by a fixed value.
     */
    std::mt19937 rng;
    DataGenerator()
    : rng(std::random_device{}())
    {}

    void printArray(const std::vector<int>& arr) {
        for (size_t i = 0; i < arr.size(); ++i) {
            std::cout << arr[i];
            if (i != arr.size() - 1) {
                std::cout << " ";
            }
        }
        std::cout << std::endl;
    }

};



#endif //DATAGENERATOR_H
