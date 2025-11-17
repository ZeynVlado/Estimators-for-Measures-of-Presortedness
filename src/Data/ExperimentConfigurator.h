

#ifndef EXPSETUPPER_H
#define EXPSETUPPER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <functional>

#include "Generation/DataGenerator.h"
#include "Sampling/Sampler.h"

class ExperimentConfigurator {
public:

    struct Config {
        int n;
        int sampleSize;
        int minValue;
        int maxValue;
        std::string root;

        std::vector<std::string> clusterGroups { "sqrt", "2sqrt", "log2", "2log2", "smlLength" };

        std::function<int(const std::string&, int /*S*/)> clusterSizer =
            [](const std::string& g, int S){
                auto clamp = [&](int v){ return std::max(1, std::min(v, S)); };
                if (g == "sqrt")     return clamp((int)std::floor(std::sqrt((double)S)));
                if (g == "2sqrt")    return clamp(2 * (int)std::floor(std::sqrt((double)S)));
                if (g == "log2")     return clamp((int)std::floor(std::log2((double)S)));
                if (g == "2log2")    return clamp(2 * (int)std::floor(std::log2((double)S)));
                if (g == "smlLength")return S;
                return clamp((int)std::floor(std::sqrt((double)S)));
            };

        std::vector<std::string> stratumGroups{
            "smlLength", "2 smlLength", "n div 4", "n div 10", "n div 20"
        };

        std::function<int(const std::string &, int /*n*/, int /*S*/)> stratumSizer =
            [](const std::string &g, int n, int S) {
                auto clamp = [&](int v){ return std::max(1, std::min(v, n)); };
                if (g == "smlLength")    return clamp(S);
                if (g == "2 smlLength")  return clamp(2*S);
                if (g == "n div 4")      return clamp(n/4);
                if (g == "n div 10")     return clamp(n/10);
                if (g == "n div 20")     return clamp(n/20);
                return clamp((int)std::floor(std::sqrt((double)n)));
            };
    };

    ExperimentConfigurator(int n, int sampleSize, int minValue, int maxValue, const std::string& root)
        : cfg_{n, sampleSize, minValue, maxValue, root} {}

    explicit ExperimentConfigurator(Config cfg) : cfg_(std::move(cfg)) {}

    void setSampleLength(int sampleLength) {
        cfg_.sampleSize = sampleLength;
    }

    void configure(int countPerSet,
                   const std::vector<int>& kValues,
                   const std::vector<int>& uValues,
                   const std::vector<int>& runsValues,
                   bool includePermutation,
                   bool includeMultisetPermutation,
                   bool includeRandomK,
                   bool includeRuns,
                   bool overwrite)

    {
        if (includePermutation) {
            std::cout << "[EXP] permutation n=" << cfg_.n
                      << " count=" << countPerSet << "\n";
            generatePermutationSet(countPerSet, overwrite);
        }

        if (includeMultisetPermutation) {
            for (int u : uValues) {
                generateMultiSetPermutationSet(u, countPerSet, overwrite);
            }

        }

        if (includeRandomK) {
            for (int k : kValues) {
                std::cout << "[EXP] random_array n=" << cfg_.n
                          << " k=" << k
                          << " count=" << countPerSet << "\n";
                generateRandomSet(k, countPerSet, overwrite);
            }
        }

        if (includeRuns) {
            for (int r : runsValues) {
                std::cout << "[EXP] run_array n=" << cfg_.n
                          << " runs=" << r
                          << " (" << runsLabel(cfg_.n, r) << ")"
                          << " count=" << countPerSet << "\n";
                generateRunsSet(r, countPerSet, overwrite);
            }
        }

        std::cout << "[OK] experiment_data_input fully populated under: " << cfg_.root << "\n";
    }

    void generatePermutationSet(int count, bool overwrite = true) const {
        RndArrGenerator gen;
        std::vector<std::vector<int>> arrays;
        arrays.reserve(count);
        for (int i = 0; i < count; ++i) {

            arrays.emplace_back(gen.generatePermutation(static_cast<size_t>(cfg_.n)));
        }
        const std::string base = join(cfg_.root, join("permutation", toStr(cfg_.n)));
        saveArraysAndAllSamples(arrays, base, overwrite);
    }

    void generateRandomSet(int k, int count, bool overwrite = true) const {
        RndArrGenerator gen;
        std::vector<std::vector<int>> arrays; arrays.reserve(count);
        for (int i = 0; i < count; ++i) {
            arrays.emplace_back(gen.generateRandom(
                static_cast<size_t>(cfg_.n),
                cfg_.minValue,
                cfg_.maxValue,
                k
            ));
        }
        const std::string base = join(cfg_.root,
                                      join("random_array", join(toStr(cfg_.n), join("k", toStr(k)))));
        saveArraysAndAllSamples(arrays, base, overwrite);
    }

    void generateMultiSetPermutationSet(int u, int count, bool overwrite = true) const {
        RndArrGenerator gen;
        std::vector<std::vector<int>> arrays;
        arrays.reserve(count);
        for (int i = 0; i < count; ++i) {

            arrays.emplace_back(gen.generateMultisetPermutation(cfg_.n, u, 1));
        }
        const std::string base = join(cfg_.root,
                                      join("multiset_permutation", join(toStr(cfg_.n), join("u", toStr(u)))));
        saveArraysAndAllSamples(arrays, base, overwrite);
    }

void generateRunsSet(int runs, int count, bool overwrite = true) const {

        // runs_perm_sort_segment
    {
        std::vector<std::vector<int>> arrays;
        arrays.reserve(count);

        int fails = 0;
        while ((int)arrays.size() < count) {
            try {
                RunGeneratorSegmentsBased gen;
                // ВАЖНО: именно generateRunsFromPermutationSortedSegments(size, runs)
                auto arr = gen.generateRunsFromPermutationSortedSegments(
                    static_cast<size_t>(cfg_.n),
                    static_cast<size_t>(runs)
                );
                if ((int)arr.size() != cfg_.n) { ++fails; continue; }
                arrays.emplace_back(std::move(arr));
            } catch (...) {
                ++fails;
                continue;
            }
        }

        const std::string base = join(cfg_.root,
            join("run_array",
            join("runs_perm_sort_segment",
            join(toStr(cfg_.n),
            join("r", runsLabel(cfg_.n, runs))))));

        saveArraysAndAllSamples(arrays, base, overwrite);
    }

    // runs_sort_arr_sort_segment
    {
        std::vector<std::vector<int>> arrays;
        arrays.reserve(count);

        int fails = 0;
        while ((int)arrays.size() < count) {
            try {
                RunGeneratorSegmentsBased gen;
                auto arr = gen.generateRunsFromSortedArrayPermutatedSegments(
                    static_cast<size_t>(cfg_.n),
                    static_cast<size_t>(runs),
                    cfg_.minValue,
                    cfg_.maxValue
                );
                if ((int)arr.size() != cfg_.n) { ++fails; continue; }
                arrays.emplace_back(std::move(arr));
            } catch (...) {
                ++fails;
                continue;
            }
        }

        const std::string base = join(cfg_.root,
            join("run_array",
            join("runs_sort_arr_sort_segment",
            join(toStr(cfg_.n),
            join("r", runsLabel(cfg_.n, runs))))));

        saveArraysAndAllSamples(arrays, base, overwrite);
    }

    // runs_noise
    {
        std::vector<std::vector<int>> arrays;
        arrays.reserve(count);

        const size_t noiseCount = (runs > 0) ? static_cast<size_t>(runs - 1) : 0;

        int fails = 0;
        while ((int)arrays.size() < count) {
            try {
                RunGeneratorNoiseBased gen;
                auto arr = gen.generateRunsFromSortedArrayNoise(
                    static_cast<size_t>(cfg_.n),
                    noiseCount,
                    cfg_.minValue,
                    cfg_.maxValue
                );
                if ((int)arr.size() != cfg_.n) { ++fails; continue; }
                arrays.emplace_back(std::move(arr));
            } catch (...) {
                ++fails;
                continue;
            }
        }

        const std::string base = join(cfg_.root,
            join("run_array",
            join("runs_noise",
            join(toStr(cfg_.n),
            join("r", runsLabel(cfg_.n, runs))))));

        saveArraysAndAllSamples(arrays, base, overwrite);
    }
}
private:
    Config cfg_;

    static std::string toStr(int v) { return std::to_string(v); }

    static std::string join(const std::string& a, const std::string& b) {
        if (a.empty()) return b;
        const char sep = '/';
        if (a.back() == '/' || a.back() == '\\') return a + b;
        return a + sep + b;
    }

    static void ensureDir(const std::string& p) {
        std::filesystem::create_directories(p);
    }

    static void writeRowCSV(std::ofstream& ofs, const std::vector<int>& v) {
        for (size_t i = 0; i < v.size(); ++i) {
            ofs << v[i];
            if (i + 1 < v.size()) ofs << ',';
        }
        ofs << '\n';
    }

    static std::string runsLabel(int n, int runs) {
        const int sqrt_n   = std::max(1, (int)std::floor(std::sqrt((double)n)));
        const int log2_n   = std::max(1, (int)std::floor(std::log2((double)n)));
        const int n_div_25 = std::max(1, (int)std::lround(n / 2.5));

        if (runs == sqrt_n)        return "sqrt n";
        if (runs == 2*sqrt_n)      return "2sqrt n";
        if (runs == log2_n)        return "log2 n";
        if (runs == n/4)           return "n div 4";
        if (runs == n/10)          return "n div 10";
        if (runs == n/20)          return "n div 20";
        if (runs == n_div_25)      return "n div 2.5";

        return "r_" + std::to_string(runs);
    }

    static std::string clusterFolderLabel(const std::string& key) {
        if (key == "sqrt")      return "sqrt smlLength";
        if (key == "2sqrt")     return "2sqrt smlLength";
        if (key == "log2")      return "log2 smlLength";
        if (key == "2log2")     return "2log2 smlLength";
        if (key == "smlLength") return "smlLength";
        return key;
    }

    static std::string stratumFolderLabel(const std::string& key) {
        return key;
    }


void saveArraysAndAllSamples(const std::vector<std::vector<int>>& arrays,
                             const std::string& baseDir,
                             bool overwrite) const
{
    // 1) arrays.csv
    const std::string arraysCsv = join(baseDir, "arrays.csv");
    ensureDir(std::filesystem::path(arraysCsv).parent_path().string());

    if (overwrite || !fileExists(arraysCsv)) {
        std::ofstream ofs(arraysCsv, std::ios::trunc);
        if (!ofs) throw std::runtime_error("Cannot open " + arraysCsv);
        for (const auto& v : arrays) writeRowCSV(ofs, v);
    }

    const std::string samplesRoot = join(join(baseDir, "samples"),
                                 sampleLenFolderLabel(cfg_.n, cfg_.sampleSize));

    // Cluster sampling
    for (const auto& cg : cfg_.clusterGroups) {
        const int clSize = cfg_.clusterSizer(cg, cfg_.sampleSize);
        const std::string dir = join(join(samplesRoot, "cluster sampling"),
                                     clusterFolderLabel(cg));
        ensureDir(dir);

        const std::string csv = join(dir, "samples.csv");
        const std::string indices_csv = join(dir, "sample_indices.csv");

        if (overwrite || !fileExists(csv) || !fileExists(indices_csv)) {
            std::ofstream ofs(csv, std::ios::trunc);
            std::ofstream ofs_indices(indices_csv, std::ios::trunc);

            if (!ofs) throw std::runtime_error("Cannot open " + csv);
            if (!ofs_indices) throw std::runtime_error("Cannot open " + indices_csv);

            for (const auto& baseArr : arrays) {

                ClusterSampler sampler(cfg_.sampleSize, clSize);
                sampler.createSample(baseArr);

                writeRowCSV(ofs, sampler.sample());
                writeRowCSV(ofs_indices, sampler.sampleIndices());
            }
        }
    }

    // Stratified sampling
    for (const auto& sg : cfg_.stratumGroups) {
        const int strSize = cfg_.stratumSizer(sg, cfg_.n, cfg_.sampleSize);
        const std::string dir = join(join(samplesRoot, "stratified sampling"),
                                     stratumFolderLabel(sg));
        ensureDir(dir);

        const std::string csv = join(dir, "samples.csv");
        const std::string indices_csv = join(dir, "sample_indices.csv");

        if (overwrite || !fileExists(csv) || !fileExists(indices_csv)) {
            std::ofstream ofs(csv, std::ios::trunc);
            std::ofstream ofs_indices(indices_csv, std::ios::trunc);

            for (const auto& baseArr : arrays) {
                StratifiedSampler sampler(cfg_.sampleSize, strSize);

                sampler.createSample(baseArr);

                writeRowCSV(ofs, sampler.sample());
                writeRowCSV(ofs_indices, sampler.sampleIndices());
            }
        }
    }

    // Combined sampling
    for (const auto& sg : cfg_.stratumGroups) {
        const int strSize = cfg_.stratumSizer(sg, cfg_.n, cfg_.sampleSize);
        const std::string dir = join(join(samplesRoot, "combined sampling"),
                                     stratumFolderLabel(sg));
        ensureDir(dir);

        const std::string csv = join(dir, "samples.csv");
        const std::string indices_csv = join(dir, "sample_indices.csv");

        if (overwrite || !fileExists(csv) || !fileExists(indices_csv)) {
            std::ofstream ofs(csv, std::ios::trunc);
            std::ofstream ofs_indices(indices_csv, std::ios::trunc);

            for (const auto& baseArr : arrays) {

                CombinedSampler sampler(cfg_.sampleSize, strSize);

                sampler.createSample(baseArr);

                writeRowCSV(ofs, sampler.sample());
                writeRowCSV(ofs_indices, sampler.sampleIndices());
            }
        }
    }

    std::cout << "[OK] arrays + samples saved under: " << baseDir << "\n";
}

    static bool fileExists(const std::string& p) {
        return std::filesystem::exists(p);
    }

    static std::string sampleLenFolderLabel(int n, int S) {
        const int sqrt_n = std::max(1, (int)std::floor(std::sqrt((double)n)));
        if (S == sqrt_n)         return "sqrt n";
        if (S == 2*sqrt_n)       return "2sqrt n";
        if (S == 3*sqrt_n)       return "3sqrt n";
        if (S == 4*sqrt_n)       return "4sqrt n";
        if (S == 5*sqrt_n)       return "5sqrt n";
        if (S == 6*sqrt_n)       return "6sqrt n";
        if (S == 7*sqrt_n)       return "7sqrt n";
        if (S == 2*std::max(1, sqrt_n)) return "2sqrt n";

        return "sml " + std::to_string(S);
    }

    static std::vector<std::vector<int>> readArraysCSVAll(const std::string& arraysCsvPath) {
        std::ifstream ifs(arraysCsvPath);
        if (!ifs) throw std::runtime_error("Cannot open arrays csv: " + arraysCsvPath);

        std::vector<std::vector<int>> arrays;
        std::string line;
        while (std::getline(ifs, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string cell;
            std::vector<int> row;
            while (std::getline(ss, cell, ',')) {
                row.push_back(cell.empty() ? 0 : std::stoi(cell));
            }
            if (!row.empty()) arrays.emplace_back(std::move(row));
        }
        return arrays;
    }
};

#endif //EXPSETUPPER_H
