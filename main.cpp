#include <iostream>
#include <fstream>
#include <random>
#include <xsimd/xsimd.hpp>
#include "fasta.hpp"
#include "alignment_result.hpp"
#include "smith_waterman.hpp"

std::string randomSequence(const std::size_t);

int main() {
    Fasta fasta1, fasta2;
    std::string input;
    std::cout << "Select input method:" << std::endl;
    std::cout << "    1. Enter FASTA file path" << std::endl;
    std::cout << "    2. Generate random sequences" << std::endl;
    std::cout << std::endl;
    std::cout << "Enter 1 or 2: ";
    std::getline(std::cin, input);
    
    if (input == "1") {
        std::string filepath;
        std::cout << "The FASTA file path: ";
        std::getline(std::cin, filepath);
        std::ifstream ifs(filepath);
        if (!ifs) {
            throw std::runtime_error("Failed to open file: " + filepath);
        }
        ifs >> fasta1 >> fasta2;
        std::cout << std::endl;
    } else if (input == "2") {
        int len1, len2;
        std::cout << "Enter two sequence lengths (separated by space): ";
        std::cin >> len1 >> len2;
        fasta1 = {"Seq1", randomSequence(len1)};
        fasta2 = {"Seq2", randomSequence(len2)};
    } else {
        throw std::runtime_error("Invalid input. Please enter 1 or 2.");
    }

    std::cout << "Show alignment results? (y/n): ";
    std::cin.ignore();
    std::getline(std::cin, input);
    bool showAlignmentResult = (input == "y" || input == "Y");

    auto report1 = SmithWaterman::run(fasta1, fasta2, SmithWaterman::Mode::Naive);
    std::cout << "[Naive] " << report1.durationMs << " ms, score: " << report1.maxScore << std::endl;
    if (showAlignmentResult) std::cout << report1.alignmentResult << std::endl;
    auto report2 = SmithWaterman::run(fasta1, fasta2, SmithWaterman::Mode::XSIMD);
    std::cout << "[XSIMD] " << report2.durationMs << " ms, score: " << report2.maxScore << std::endl;
    if (showAlignmentResult) std::cout << report2.alignmentResult << std::endl;
    auto report3 = SmithWaterman::run(fasta1, fasta2, SmithWaterman::Mode::CUDA);
    std::cout << "[CUDA] " << report3.durationMs << " ms, score: " << report3.maxScore << std::endl;
    if (showAlignmentResult) std::cout << report3.alignmentResult << std::endl;

    return 0;
}

std::string randomSequence(const std::size_t n) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    static thread_local std::mt19937 rng{std::random_device{}()};
    static thread_local std::uniform_int_distribution<int> dist(0, 3);

    std::string result;
    result.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        result += bases[dist(rng)];
    return result;
}
