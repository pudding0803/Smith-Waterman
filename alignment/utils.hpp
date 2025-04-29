#pragma once

#include <vector>
#include <exception>
#include "AlignmentResult.hpp"
#include "Fasta.hpp"
#include "Direction.hpp"

constexpr uint8_t baseToIndex(const char base) {
    switch (base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
    }
    throw std::invalid_argument("Invalid nucleotide character");
}

AlignmentResult traceback(const std::vector<std::vector<Direction>>& trace, const Fasta& fasta1, const Fasta& fasta2, std::pair<std::size_t, std::size_t>& maxPos) {
    AlignmentResult result;
    auto [i, j] = maxPos;
    int endI = i, endJ = j;

    while (i > 0 && j > 0) {
        if (trace[i][j] == Direction::DIAG) {
            const char a = fasta1.sequence[--i];
            const char b = fasta2.sequence[--j];
            result.aligned1 = a + result.aligned1;
            result.aligned2 = b + result.aligned2;
            result.alignment = (a == b ? '|' : '*') + result.alignment;
        } else if (trace[i][j] == Direction::LEFT) {
            result.aligned1 = '-' + result.aligned1;
            result.aligned2 = fasta2.sequence[--j] + result.aligned2;
            result.alignment = ' ' + result.alignment;
        } else if (trace[i][j] == Direction::UP) {
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