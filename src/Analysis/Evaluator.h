

#ifndef CHARTBUILDER_H
#define CHARTBUILDER_H
#include <filesystem>
#include <fstream>
#include <vector>
#include "PresortednessMetrics.h"

#include "../Data/Sampling/Sampler.h"

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

#include <system_error>

#include "MetricNormalizer.h"


namespace fs = std::filesystem;

class Evaluator {
public:
    Evaluator()  = default;
    ~Evaluator() = default;

    void prepareOutputStructure(const std::string& inputRoot,
                                const std::string& outputRoot) const
    {
        if (!fs::exists(inputRoot)) {
            throw std::runtime_error("Input root does not exist: " + inputRoot);
        }
        ensureDir(outputRoot);

        for (const auto& entry : fs::recursive_directory_iterator(inputRoot)) {
            if (entry.is_directory()) {
                const fs::path rel     = fs::relative(entry.path(), inputRoot);
                const fs::path outDir  = fs::path(outputRoot) / rel;
                ensureDir(outDir.string());

                if (hasFile(entry.path(), "samples.csv")) {
                    const fs::path sampleMetrics = outDir / "sample_metrics.csv";
                    createFileIfNotExists(sampleMetrics.string());
                }
            } else if (entry.is_regular_file() && entry.path().filename() == "arrays.csv") {
                const fs::path outDir        = fs::path(outputRoot) /
                                               fs::relative(entry.path().parent_path(), inputRoot);
                const fs::path arraysMetrics = outDir / "arrays_metrics.csv";
                createFileIfNotExists(arraysMetrics.string());
            }
        }

        std::cout << "[OK] Prepared output structure at: " << outputRoot << "\n";
    }

    void evaluateAll(const std::string& inputRoot,
                 const std::string& outputRoot,
                 bool overwrite) const
{
    if (!fs::exists(inputRoot)) {
        throw std::runtime_error("Input root does not exist: " + inputRoot);
    }
    ensureDir(outputRoot);

    for (const auto& entry : fs::recursive_directory_iterator(inputRoot)) {
        const fs::path rel = fs::relative(entry.path(), inputRoot);
        const fs::path out = fs::path(outputRoot) / rel;

        if (entry.is_regular_file() && entry.path().filename() == "arrays.csv") {
            const fs::path outDir         = out.parent_path();
            const fs::path arraysMetrics  = outDir / "arrays_metrics.csv";
            ensureDir(outDir.string());

            if (!overwrite && fs::exists(arraysMetrics)) {
                std::cout << "[SKIP] " << arraysMetrics.string()
                          << " (exists; overwrite=false)\n";
            } else {
                std::cout << "[WRITE] " << arraysMetrics.string()
                          << (overwrite ? " (overwrite)\n" : " (create)\n");
                evaluateCSVtoNormMetrics(entry.path().string(),
                                         arraysMetrics.string());
            }
            continue;
        }

        if (entry.is_directory() && hasFile(entry.path(), "samples.csv")) {
            const fs::path samplesCsv     = entry.path() / "samples.csv";
            const fs::path outDir         = fs::path(outputRoot) / fs::relative(entry.path(), inputRoot);
            const fs::path sampleMetrics  = outDir / "sample_metrics.csv";
            ensureDir(outDir.string());

            if (!overwrite && fs::exists(sampleMetrics)) {
                std::cout << "[SKIP] " << sampleMetrics.string()
                          << " (exists; overwrite=false)\n";
            } else {
                std::cout << "[WRITE] " << sampleMetrics.string()
                          << (overwrite ? " (overwrite)\n" : " (create)\n");
                evaluateCSVtoNormMetrics(samplesCsv.string(),
                                         sampleMetrics.string());
            }
            continue;
        }
    }

    std::cout << "[OK] Full evaluation finished. Output at: " << outputRoot << "\n";
}

private:

    static void ensureDir(const std::string& path) {
        std::error_code ec;
        fs::create_directories(path, ec);
        if (ec) throw std::runtime_error("Failed to create directories: " + path + " (" + ec.message() + ")");
    }

    static bool hasFile(const fs::path& dir, const std::string& filename) {
        if (!fs::exists(dir) || !fs::is_directory(dir)) return false;
        for (const auto& it : fs::directory_iterator(dir)) {
            if (it.is_regular_file() && it.path().filename() == filename) return true;
        }
        return false;
    }

    static void createFileIfNotExists(const std::string& path) {
        if (!fs::exists(path)) {
            std::ofstream ofs(path);
            if (!ofs) throw std::runtime_error("Cannot create file: " + path);
        }
    }


    static bool readNextRow(std::istream& is, std::vector<int>& out) {
        out.clear();
        std::string line;
        if (!std::getline(is, line)) return false;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            out.push_back(cell.empty() ? 0 : std::stoi(cell));
        }
        return true;
    }

    static void writeHeaderNorm(std::ofstream& ofs) {
        ofs << "n,"
            << "inv_norm,"
            << "runs_norm,"
            << "rem_norm,"
            << "osc_norm,"
            << "dis_norm,"
            << "ham_norm,"
            << "run_entr_norm,"
            << "val_freq_entr_norm,"
            << "distinct_norm,"
            << "distinct_entr_norm"
            << "\n";
    }

    static void evaluateCSVtoNormMetrics(const std::string& inputCsv,
                                         const std::string& outputCsv)
    {
        std::ifstream ifs(inputCsv);
        if (!ifs) {
            throw std::runtime_error("Cannot open input csv: " + inputCsv);
        }

        std::ofstream ofs(outputCsv, std::ios::trunc);
        if (!ofs) {
            throw std::runtime_error("Cannot open output csv: " + outputCsv);
        }

        writeHeaderNorm(ofs);

        std::vector<int> row;
        while (readNextRow(ifs, row)) {
            if (row.empty()) continue;
            writeNormMetricsLine(ofs, row);
        }
    }


  static void writeNormMetricsLine(std::ofstream& ofs, const std::vector<int>& a) {

    const long long n = static_cast<long long>(a.size());

    const long long runs     = PresortednessMetrics::instance().calculateRuns(a);
    const long long inv      = PresortednessMetrics::instance().calculateInversions(a);
    const long long rem      = PresortednessMetrics::instance().calculateRem(a);
    const long long osc      = PresortednessMetrics::instance().calculateOsc(a);
    const long long dis      = PresortednessMetrics::instance().calculateDis(a);
    const long long ham      = PresortednessMetrics::instance().calculateHam(a);

    const double    runsEntr = PresortednessMetrics::instance().calculateRunsEntropy(a);
    const long long distinct = PresortednessMetrics::instance().calculateDistinctElements(a);
    const double    valFreqEntr = PresortednessMetrics::instance().calculateValueFrequencyEntropy(a);
    const double    distinctEntr = PresortednessMetrics::instance().calculateDistinctElementEntropy(a);

    const double runsNorm      = MetricNormalizer::instance().normalizeRuns(runs, n);
    const double invNorm       = MetricNormalizer::instance().normalizeInversions(inv, n);
    const double remNorm       = MetricNormalizer::instance().normalizeRem(rem, n);
    const double oscNorm       = MetricNormalizer::instance().normalizeOsc(osc, n);
    const double disNorm       = MetricNormalizer::instance().normalizeDis(dis, n);
    const double hamNorm       = MetricNormalizer::instance().normalizeHam(ham, n);

    const double runsEntropyNorm      = MetricNormalizer::instance().normalizeRunsEntropy(runsEntr, runs);
    const double valFreqEntropyNorm   = MetricNormalizer::instance().normalizeValueFrequencyEntropy(valFreqEntr, distinct);
    const double distinctNorm       = MetricNormalizer::instance().normalizeDistinct(distinct, n);
    const double distinctEntropyNorm  = MetricNormalizer::instance().normalizeDistinctElementEntropy(distinctEntr, n);

    ofs << n
        << ',' << invNorm
        << ',' << runsNorm
        << ',' << remNorm
        << ',' << oscNorm
        << ',' << disNorm
        << ',' << hamNorm
        << ',' << runsEntropyNorm
        << ',' << valFreqEntropyNorm
        << ',' << distinctNorm
        << ',' << distinctEntropyNorm
        << '\n';
}


};
#endif //CHARTBUILDER_H
