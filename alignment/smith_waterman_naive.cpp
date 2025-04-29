#include <vector>
#include <algorithm>
#include <chrono>
#include "smith_waterman.hpp"
#include "direction.hpp"
#include "constants.hpp"

SmithWaterman::AlignmentReport SmithWaterman::naive(const Fasta& fasta1, const Fasta& fasta2) {
    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<value_type>> dp(fasta1.size + 1, std::vector<value_type>(fasta2.size + 1));
    std::vector<std::vector<value_type>> vGap(fasta1.size + 1, std::vector<value_type>(fasta2.size + 1));
    std::vector<std::vector<value_type>> hGap(fasta1.size + 1, std::vector<value_type>(fasta2.size + 1));
    std::vector<std::vector<Direction>> trace(fasta1.size + 1, std::vector<Direction>(fasta2.size + 1));

    value_type maxScore = 0;
    std::pair<std::size_t, std::size_t> maxPos = {0, 0};
    for (std::size_t i = 1; i <= fasta1.size; ++i) {
        for (std::size_t j = 1; j <= fasta2.size; ++j) {
            vGap[i][j] = std::max(dp[i - 1][j] + GAP_OPEN, vGap[i - 1][j] + GAP_EXTEND);
            hGap[i][j] = std::max(dp[i][j - 1] + GAP_OPEN, hGap[i][j - 1] + GAP_EXTEND);
    
            value_type diag = dp[i - 1][j - 1] + (fasta1.sequence[i - 1] == fasta2.sequence[j - 1] ? MATCH : MISMATCH);
            dp[i][j] = std::max({static_cast<value_type>(0), diag, vGap[i][j], hGap[i][j]});
    
            if (dp[i][j] == diag) [[likely]] {
                trace[i][j] = Direction::Diag;
            } else if (dp[i][j] == vGap[i][j]) {
                trace[i][j] = Direction::Up;
            } else if (dp[i][j] == hGap[i][j]) {
                trace[i][j] = Direction::Left;
            }  else [[unlikely]] {
                trace[i][j] = Direction::None;
            }
    
            if (dp[i][j] > maxScore) {
                maxScore = dp[i][j];
                maxPos = {i, j};
            }
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto result = traceback(trace, fasta1, fasta2, maxPos);
    return {result, std::chrono::duration<double, std::milli>(end_time - start_time).count(), maxScore};
}
