#include <iostream>
#include <fstream>
#include <xsimd/xsimd.hpp>
#include "Fasta.hpp"
#include "AlignmentResult.hpp"
#include "SmithWaterman.hpp"

int main() {
    Fasta fasta1, fasta2;
    std::string filepath;
    std::cout << "The FASTA file path: ";
    getline(std::cin, filepath);
    std::ifstream ifs(filepath);
    if (!ifs) {
        throw std::runtime_error("Failed to open file: " + filepath);
    }
    ifs >> fasta1 >> fasta2;
    std::cout << std::endl;
    auto report1 = SmithWaterman::run(fasta1, fasta2, SmithWaterman::Mode::Naive);
    auto report2 = SmithWaterman::run(fasta1, fasta2, SmithWaterman::Mode::SIMD);
    std::cout << "[Naive] " << report1.durationMs << " ms, score: " << report1.maxScore << std::endl << report1.alignmentResult << std::endl;
    std::cout << "[SIMD] " << report2.durationMs << " ms, score: " << report2.maxScore << std::endl << report2.alignmentResult << std::endl;
    std::cout << "Speedup: " << report1.durationMs / report2.durationMs << " x" << std::endl;
    return 0;
}
