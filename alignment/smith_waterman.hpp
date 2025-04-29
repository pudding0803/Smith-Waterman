#pragma once

#include <vector>
#include "alignment_result.hpp"
#include "fasta.hpp"
#include "direction.hpp"

class SmithWaterman {
public:
    using value_type = int16_t;

    enum class Mode {
        Naive, XSIMD
    };

    struct AlignmentReport {
        AlignmentResult alignmentResult;
        double durationMs;
        value_type maxScore;
    };

    [[nodiscard]] static AlignmentReport run(const Fasta&, const Fasta&, Mode);

private:
    static AlignmentResult traceback(const std::vector<std::vector<Direction>>&, const Fasta&, const Fasta&, std::pair<std::size_t, std::size_t>&);

    static AlignmentReport naive(const Fasta&, const Fasta&);
    static AlignmentReport xsimd(const Fasta&, const Fasta&);
};
