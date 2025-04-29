#include <vector>
#include <chrono>
#include <xsimd/xsimd.hpp>
#include "SmithWaterman.hpp"
#include "utils.hpp"
#include "Direction.hpp"
#include "Constants.hpp"

SmithWaterman::AlignmentReport SmithWaterman::run(const Fasta& fasta1, const Fasta& fasta2, Mode mode) {
    switch (mode) {
        case Mode::Naive:
            return naive(fasta1, fasta2);
        case Mode::SIMD:
            return simd(fasta1, fasta2);
    }
    throw std::invalid_argument("Invalid SmithWaterman::Mode");
}

SmithWaterman::AlignmentReport SmithWaterman::naive(const Fasta& fasta1, const Fasta& fasta2) {
    using value_type = int16_t;
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
            dp[i][j] = std::max({static_cast<value_type>(0), diag, hGap[i][j], vGap[i][j]});
    
            if (dp[i][j] == diag) [[likely]] {
                trace[i][j] = Direction::DIAG;
            } else if (dp[i][j] == vGap[i][j]) {
                trace[i][j] = Direction::UP;
            } else if (dp[i][j] == hGap[i][j]) {
                trace[i][j] = Direction::LEFT;
            }  else [[unlikely]] {
                trace[i][j] = Direction::NONE;
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

SmithWaterman::AlignmentReport SmithWaterman::simd(const Fasta& fasta1, const Fasta& fasta2) {
    using batch_t = xsimd::batch<value_type>;

    auto start_time = std::chrono::high_resolution_clock::now();

    std::size_t segLen = (fasta1.size + batch_t::size - 1) / batch_t::size;
    std::array<std::vector<batch_t>, 4> profile;

    batch_t vMatch = batch_t::broadcast(MATCH);
    batch_t vMismatch = batch_t::broadcast(MISMATCH);

    for (const auto& base : {'A', 'C', 'G', 'T'}) {
        std::size_t index = baseToIndex(base); 
        profile[index].resize(segLen);
        const auto vBase = batch_t::broadcast(base);
    
        for (std::size_t i = 0; i < segLen; ++i) {
            std::size_t j = i;
            value_type buffer[batch_t::size];
            for (std::size_t k = 0; k < batch_t::size; ++k, j += segLen) {
                buffer[k] = (j < fasta1.size) ? fasta1.sequence[j] : 0;
            }
            auto sequence = batch_t::load_unaligned(buffer);
            profile[index][i] = xsimd::select(sequence == vBase, vMatch, vMismatch);
        }
    }

    const batch_t zero = batch_t::broadcast(0);
    const batch_t vGapOpen = batch_t::broadcast(GAP_OPEN);
    const batch_t vGapExtend = batch_t::broadcast(GAP_EXTEND);
    const batch_t vDirectionNone = batch_t::broadcast(Direction::NONE);
    const batch_t vDirectionDiag = batch_t::broadcast(Direction::DIAG);
    const batch_t vDirectionUp = batch_t::broadcast(Direction::UP);
    const batch_t vDirectionLeft = batch_t::broadcast(Direction::LEFT);

    std::vector<batch_t> vPrevH(segLen, zero);
    std::vector<batch_t> vCurrH(segLen, zero);
    std::vector<batch_t> vE(segLen, zero);
    batch_t vMax = zero;

    value_type maxScore = 0;
    std::pair<std::size_t, std::size_t> maxPos = {0, 0};
    std::vector<std::vector<Direction>> trace(fasta1.size + 1, std::vector<Direction>(fasta2.size + 1, Direction::NONE));

    for (std::size_t j = 0; j < fasta2.size; ++j) {
        const auto& currProfile = profile[baseToIndex(fasta2.sequence[j])];
        std::vector<batch_t> vF(segLen, zero);
        std::vector<batch_t> vTrace(segLen, zero);
        
        for (std::size_t seg = 0; seg < segLen; ++seg) {
            batch_t vPrevDiag = seg > 0 ? vPrevH[seg - 1] : xsimd::slide_left<sizeof(value_type)>(vPrevH.back());
            batch_t vCurrDiag = vPrevDiag + currProfile[seg];
            vE[seg] = xsimd::max(vPrevH[seg] + vGapOpen, vE[seg] + vGapExtend);
            vCurrH[seg] = xsimd::max(vCurrDiag, vE[seg]);
            vCurrH[seg] = xsimd::max(vCurrH[seg], zero);

            vTrace[seg] = xsimd::select(vCurrH[seg] == vE[seg], vDirectionLeft, vTrace[seg]);
            vTrace[seg] = xsimd::select(vCurrH[seg] == vCurrDiag, vDirectionDiag, vTrace[seg]);
        }

        bool done = false;
        while (!done) {
            done = true;
            for (std::size_t seg = 0; seg < segLen; ++seg) {
                batch_t vCandidateF = seg > 0
                    ? xsimd::max(vCurrH[seg - 1] + vGapOpen, vF[seg - 1] + vGapExtend)
                    : xsimd::max(xsimd::slide_left<sizeof(value_type)>(vCurrH.back()) + vGapOpen, xsimd::slide_left<sizeof(value_type)>(vF.back()) + vGapExtend);
                auto mask = vCandidateF > vCurrH[seg];
                if (xsimd::any(mask)) {
                    done = false;
                }
                vF[seg] = vCandidateF;
                vCurrH[seg] = xsimd::max(vCurrH[seg], vF[seg]);
                vTrace[seg] = xsimd::select(mask || vCandidateF == zero && vTrace[seg] == vDirectionNone, vDirectionUp, vTrace[seg]);
            }
        }

        for (std::size_t seg = 0; seg < segLen; ++seg) {
            for (std::size_t lane = 0; lane < batch_t::size; ++lane) {
                std::size_t i = seg + lane * segLen;
                if (i < fasta1.size) {
                    value_type score = vCurrH[seg].get(lane);
                    if (score > maxScore) {
                        maxScore = score;
                        maxPos = {i + 1, j + 1};
                    }
                    trace[i + 1][j + 1] = static_cast<Direction>(vTrace[seg].get(lane));
                }
            }
        }

        std::swap(vPrevH, vCurrH);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto result = traceback(trace, fasta1, fasta2, maxPos);
    return {result, std::chrono::duration<double, std::milli>(end_time - start_time).count(), maxScore};
}
