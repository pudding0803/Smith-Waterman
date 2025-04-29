#pragma once

#include <exception>
#include "AlignmentResult.hpp"
#include "Fasta.hpp"

class SmithWaterman {
public:
    using value_type = int;

    enum class Mode {
        Naive, SIMD
    };

    struct AlignmentReport {
        AlignmentResult alignmentResult;
        double durationMs;
        value_type maxScore;
    };

    [[nodiscard]] static AlignmentReport run(const Fasta&, const Fasta&, Mode);

private:
    static AlignmentReport naive(const Fasta&, const Fasta&);
    static AlignmentReport simd(const Fasta&, const Fasta&);
};
