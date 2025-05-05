#include "smith_waterman.hpp"

SmithWaterman::AlignmentReport SmithWaterman::run(const Fasta& fasta1, const Fasta& fasta2, Mode mode) {
    switch (mode) {
        case Mode::Naive:
            return naive(fasta1, fasta2);
        case Mode::XSIMD:
            return xsimd(fasta1, fasta2);
        case Mode::CUDA:
            return cuda(fasta1, fasta2);
    }
    throw std::invalid_argument("Invalid SmithWaterman::Mode");
}
