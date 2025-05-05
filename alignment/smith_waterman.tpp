template <typename GetTrace>
AlignmentResult SmithWaterman::traceback(
    const Fasta& fasta1,
    const Fasta& fasta2,
    std::pair<std::size_t, std::size_t>& maxPos,
    GetTrace&& getTrace
) noexcept {
    AlignmentResult result;
    auto [i, j] = maxPos;
    int endI = i, endJ = j;

    while (i > 0 && j > 0) {
        auto direction = getTrace(i, j);
        if (direction == Direction::Diag) [[likely]] {
            const char a = fasta1.sequence[--i];
            const char b = fasta2.sequence[--j];
            result.aligned1 = a + result.aligned1;
            result.aligned2 = b + result.aligned2;
            result.alignment = (a == b ? '|' : '*') + result.alignment;
        } else if (direction == Direction::Left) {
            result.aligned1 = '-' + result.aligned1;
            result.aligned2 = fasta2.sequence[--j] + result.aligned2;
            result.alignment = ' ' + result.alignment;
        } else if (direction == Direction::Up) {
            result.aligned1 = fasta1.sequence[--i] + result.aligned1;
            result.aligned2 = '-' + result.aligned2;
            result.alignment = ' ' + result.alignment;
        } else [[unlikely]] {
            break;
        }
    }

    result.range1 = {i, endI - 1};
    result.range2 = {j, endJ - 1};

    return result;
}
