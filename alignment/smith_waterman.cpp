#include <vector>
#include <xsimd/xsimd.hpp>
#include "smith_waterman.hpp"
#include "direction.hpp"
#include "constants.hpp"

SmithWaterman::AlignmentReport SmithWaterman::run(const Fasta& fasta1, const Fasta& fasta2, Mode mode) {
    switch (mode) {
        case Mode::Naive:
            return naive(fasta1, fasta2);
        case Mode::XSIMD:
            return xsimd(fasta1, fasta2);
    }
    throw std::invalid_argument("Invalid SmithWaterman::Mode");
}

AlignmentResult SmithWaterman::traceback(const std::vector<std::vector<Direction>>& trace, const Fasta& fasta1, const Fasta& fasta2, std::pair<std::size_t, std::size_t>& maxPos) {
    AlignmentResult result;
    auto [i, j] = maxPos;
    int endI = i, endJ = j;

    while (i > 0 && j > 0) {
        if (trace[i][j] == Direction::Diag) {
            const char a = fasta1.sequence[--i];
            const char b = fasta2.sequence[--j];
            result.aligned1 = a + result.aligned1;
            result.aligned2 = b + result.aligned2;
            result.alignment = (a == b ? '|' : '*') + result.alignment;
        } else if (trace[i][j] == Direction::Left) {
            result.aligned1 = '-' + result.aligned1;
            result.aligned2 = fasta2.sequence[--j] + result.aligned2;
            result.alignment = ' ' + result.alignment;
        } else if (trace[i][j] == Direction::Up) {
            result.aligned1 = fasta1.sequence[--i] + result.aligned1;
            result.aligned2 = '-' + result.aligned2;
            result.alignment = ' ' + result.alignment;
        } else {
            break;
        }
    }

    result.range1 = {i, endI - 1};
    result.range2 = {j, endJ - 1};

    return result;
}
