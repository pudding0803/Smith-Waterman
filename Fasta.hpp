#pragma once

#include <string>
#include <fstream>
#include <stdexcept>

struct Fasta {
    std::string description;
    std::string sequence;
    std::size_t size{};

    friend std::istream& operator>>(std::istream&, Fasta&);
};
