#include <vector>
#include <algorithm>
#include <chrono>
#include "smith_waterman.hpp"
#include "direction.hpp"
#include "constants.hpp"

SmithWaterman::AlignmentReport SmithWaterman::naive(const Fasta& fasta1, const Fasta& fasta2) {
    auto startTime = std::chrono::high_resolution_clock::now();

    std::vector<int> vPrevH(fasta1.size + 1);
    std::vector<int> vCurrH(fasta1.size + 1);
    std::vector<int> vE(fasta1.size + 1);
    std::vector<std::vector<Direction>> trace(fasta2.size + 1, std::vector<Direction>(fasta1.size + 1, Direction::None));

    int maxScore = 0;
    std::pair<std::size_t, std::size_t> maxPos = {0, 0};

    std::vector<int> dp(fasta1.size * fasta2.size);

    for (std::size_t j = 1; j <= fasta2.size; ++j) {
        int F = 0;

        for (std::size_t i = 1; i <= fasta1.size; ++i) {
            vE[i] = std::max(vPrevH[i] + GAP_OPEN, vE[i] + GAP_EXTEND);
            F = std::max(vCurrH[i - 1] + GAP_OPEN, F + GAP_EXTEND);
    
            int diag = vPrevH[i - 1] + (fasta1.sequence[i - 1] == fasta2.sequence[j - 1] ? MATCH : MISMATCH);
            vCurrH[i] = std::max({static_cast<int>(0), diag, vE[i], F});
            dp[(i - 1) * fasta2.size + (j - 1)] = vCurrH[i];

            if (vCurrH[i] == diag) [[likely]] {
                trace[j][i] = Direction::Diag;
            } else if (vCurrH[i] == vE[i]) {
                trace[j][i] = Direction::Left;
            } else if (vCurrH[i] == F) {
                trace[j][i] = Direction::Up;
            } else [[unlikely]] {
                trace[j][i] = Direction::None;
            }
    
            if (vCurrH[i] > maxScore) {
                maxScore = vCurrH[i];
                maxPos = {i, j};
            }
        }
        std::swap(vPrevH, vCurrH);
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    auto result = traceback(fasta1, fasta2, maxPos, [&](std::size_t i, std::size_t j) {
        return trace[j][i];
    });
    return {result, std::chrono::duration<double, std::milli>(endTime - startTime).count(), maxScore};
}
