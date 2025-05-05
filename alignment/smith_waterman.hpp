#pragma once

#include <vector>
#include "alignment_result.hpp"
#include "fasta.hpp"
#include "direction.hpp"

class SmithWaterman {
public:
    enum class Mode {
        Naive, XSIMD, CUDA
    };

    struct AlignmentReport {
        AlignmentResult alignmentResult;
        double durationMs;
        int maxScore;
    };

    [[nodiscard]] static AlignmentReport run(const Fasta&, const Fasta&, Mode);

private:
    static inline u_int8_t baseToIndex(const char base) {
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

    template <typename GetTrace>
    static AlignmentResult traceback(const Fasta&, const Fasta&, std::pair<std::size_t, std::size_t>&, GetTrace&&) noexcept;

    static AlignmentReport naive(const Fasta&, const Fasta&);
    static AlignmentReport xsimd(const Fasta&, const Fasta&);
    static AlignmentReport cuda(const Fasta&, const Fasta&);
};

#include "smith_waterman.tpp"
