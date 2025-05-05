#include <vector>
#include <chrono>
#include "smith_waterman.hpp"
#include "constants.hpp"

namespace {

constexpr std::size_t THREADS_PER_BLOCK = 256;

struct Score {
    int score{};
    int row{};
    int col{};
};

__device__ Score dMaxScore;

__global__ void computeWavefrontDiag(
    int* __restrict__ currH,
    const int* __restrict__ prevH,
    const int* __restrict__ prevPrevH,
    int* __restrict__ currE,
    int* __restrict__ currF,
    const int* __restrict__ prevF,
    Direction* __restrict__ trace,
    const char* __restrict__ seq1,
    const char* __restrict__ seq2,
    std::size_t seq1Size,
    std::size_t seq2Size,
    std::size_t k
) {
    __shared__ Score warpMax[THREADS_PER_BLOCK / 32];

    std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    std::size_t startRow = (k < seq2Size) ? 0 : k - seq2Size + 1;
    std::size_t row = startRow + i;
    std::size_t col = k - row;
    if (row >= seq1Size || col >= seq2Size) return;

    int offset = row + col - seq2Size;
    offset = (offset < -1) ? -1 : (offset > 1 ? 1 : offset);

    int diagI = i + offset;
    int leftI = i + (offset >= 0);
    int upI = leftI - 1;

    int diag = (diagI >= 0 ? prevPrevH[diagI] : 0) + (seq1[row] == seq2[col] ? MATCH : MISMATCH);
    currE[i] = max(prevH[leftI] + GAP_OPEN, currE[leftI] + GAP_EXTEND);
    currF[i] = upI >= 0 ? max(prevH[upI] + GAP_OPEN, prevF[upI] + GAP_EXTEND) : max(GAP_OPEN, GAP_EXTEND);
    currH[i] = max(max(0, diag), max(currE[i], currF[i]));

    std::size_t traceIdx = row * seq2Size + col;
    if (currH[i] == diag) [[likely]] {
        trace[traceIdx] = Direction::Diag;
    } else if (currH[i] == currE[i]) {
        trace[traceIdx] = Direction::Left;
    } else if (currH[i] == currF[i]) {
        trace[traceIdx] = Direction::Up;
    } else [[unlikely]] {
        trace[traceIdx] = Direction::None;
    }

    Score localMax = {currH[i], static_cast<int>(row), static_cast<int>(col)};
    int laneId = threadIdx.x % warpSize;
    unsigned mask = 0xffffffff;
    
    for (int offset = 16; offset > 0; offset /= 2) {
        int otherScore = __shfl_xor_sync(mask, localMax.score, offset);
        int otherRow = __shfl_xor_sync(mask, localMax.row, offset);
        int otherCol = __shfl_xor_sync(mask, localMax.col, offset);
        if (otherScore > localMax.score) {
            localMax.score = otherScore;
            localMax.row = otherRow;
            localMax.col = otherCol;
        }
    }
    
    int warpId = threadIdx.x / 32;
    
    if (laneId == 0) {
        warpMax[warpId] = localMax;
    }
    __syncthreads();
    
    if (threadIdx.x < THREADS_PER_BLOCK / 32) {
        localMax = warpMax[threadIdx.x];
    }
    __syncthreads();
    
    if (threadIdx.x < 32) {
        for (int offset = 16; offset > 0; offset /= 2) {
            int otherScore = __shfl_xor_sync(mask, localMax.score, offset);
            int otherRow = __shfl_xor_sync(mask, localMax.row, offset);
            int otherCol = __shfl_xor_sync(mask, localMax.col, offset);
            if (otherScore > localMax.score) {
                localMax.score = otherScore;
                localMax.row = otherRow;
                localMax.col = otherCol;
            }
        }
    }
    
    if (threadIdx.x == 0) {
        int old = atomicMax(&dMaxScore.score, localMax.score);
        if (localMax.score > old) {
            atomicExch(&dMaxScore.row, localMax.row);
            atomicExch(&dMaxScore.col, localMax.col);
        }
    }
}

}

SmithWaterman::AlignmentReport SmithWaterman::cuda(const Fasta& fasta1, const Fasta& fasta2) {
    cudaFree(0);
    auto startTime = std::chrono::high_resolution_clock::now();

    const std::size_t maxDiagLen = min(fasta1.size, fasta2.size);
    const std::size_t diagNum = fasta1.size + fasta2.size - 1;

    char *dSeq1, *dSeq2;

    cudaMalloc(&dSeq1, fasta1.size * sizeof(char));
    cudaMalloc(&dSeq2, fasta2.size * sizeof(char));

    cudaMemcpy(dSeq1, fasta1.sequence.data(), fasta1.size * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dSeq2, fasta2.sequence.data(), fasta2.size * sizeof(char), cudaMemcpyHostToDevice);

    int *dCurrH, *dPrevH, *dPrevPrevH, *dCurrE, *dCurrF, *dPrevF;
    Direction* dTrace;

    cudaMalloc(&dCurrH, maxDiagLen * sizeof(int));
    cudaMalloc(&dPrevH, maxDiagLen * sizeof(int));
    cudaMalloc(&dPrevPrevH, maxDiagLen * sizeof(int));
    cudaMalloc(&dCurrE, maxDiagLen * sizeof(int));
    cudaMalloc(&dCurrF, maxDiagLen * sizeof(int));
    cudaMalloc(&dPrevF, maxDiagLen * sizeof(int));
    cudaMalloc(&dTrace, fasta1.size * fasta2.size * sizeof(Direction));

    cudaMemset(dCurrH, 0, maxDiagLen * sizeof(int));
    cudaMemset(dPrevH, 0, maxDiagLen * sizeof(int));
    cudaMemset(dPrevPrevH, 0, maxDiagLen * sizeof(int));
    cudaMemset(dCurrE, 0, maxDiagLen * sizeof(int));
    cudaMemset(dCurrF, 0, maxDiagLen * sizeof(int));
    cudaMemset(dPrevF, 0, maxDiagLen * sizeof(int));

    dim3 gridDim((maxDiagLen + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
    dim3 blockDim(THREADS_PER_BLOCK);
    int* temp;

    for (std::size_t k = 0; k < diagNum; ++k) {
        computeWavefrontDiag<<<gridDim, blockDim>>>(
            dCurrH, dPrevH, dPrevPrevH, dCurrE, dCurrF, dPrevF,
            dTrace, dSeq1, dSeq2, fasta1.size, fasta2.size, k
        );
        cudaDeviceSynchronize();
    
        temp = dPrevPrevH;
        dPrevPrevH = dPrevH;
        dPrevH = dCurrH;
        dCurrH = temp;
        temp = dPrevF;
        dPrevF = dCurrF;
        dCurrF = temp;
    }

    std::vector<Direction> hTrace(fasta1.size * fasta2.size);
    cudaMemcpy(hTrace.data(), dTrace, hTrace.size() * sizeof(Direction), cudaMemcpyDeviceToHost);

    Score hMaxScore;
    cudaMemcpyFromSymbol(&hMaxScore, dMaxScore, sizeof(Score), 0, cudaMemcpyDeviceToHost);

    int maxScore = hMaxScore.score;
    std::pair<std::size_t, std::size_t> maxPos{static_cast<std::size_t>(hMaxScore.row + 1), static_cast<std::size_t>(hMaxScore.col + 1)};

    cudaFree(dSeq1);
    cudaFree(dSeq2);
    cudaFree(dCurrH);
    cudaFree(dPrevH);
    cudaFree(dPrevPrevH);
    cudaFree(dCurrE);
    cudaFree(dCurrF);
    cudaFree(dPrevF);
    cudaFree(dTrace);

    auto endTime = std::chrono::high_resolution_clock::now();
    auto result = traceback(fasta1, fasta2, maxPos, [&](std::size_t i, std::size_t j) {
        return hTrace[(i - 1) * fasta2.size + (j - 1)];
    });
    return {result, std::chrono::duration<double, std::milli>(endTime - startTime).count(), maxScore};
}
