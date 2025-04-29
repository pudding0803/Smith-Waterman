#include <string>
#include <iomanip>
#include "AlignmentResult.hpp"

std::ostream& operator<<(std::ostream& os, const AlignmentResult& result) {
    std::size_t padding = std::to_string(std::max(result.range1.first, result.range2.first) + 1).size();
    os << "Seq1:  " << std::setw(padding) << result.range1.first + 1 << "  " << result.aligned1 << "  " << result.range1.second << std::endl;
    os << "         " << std::string(padding, ' ') << result.alignment << std::endl;
    os << "Seq2:  " << std::setw(padding) << result.range2.first + 1 << "  " << result.aligned2 << "  " << result.range2.second << std::endl;
    return os;
}
