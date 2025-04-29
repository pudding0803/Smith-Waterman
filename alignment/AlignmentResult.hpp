#pragma once

#include <iostream>

struct AlignmentResult {
    std::string aligned1;
    std::string aligned2;
    std::pair<std::size_t, std::size_t> range1;
    std::pair<std::size_t, std::size_t> range2;
    std::string alignment;
};

std::ostream& operator<<(std::ostream&, const AlignmentResult&);
