#pragma once

#include <string>
#include <fstream>
#include <stdexcept>

struct Fasta {
    std::string description;
    std::string sequence;
    std::size_t size{};

    Fasta() = default;
    
    Fasta(const std::string& desc, const std::string& seq) : description(desc), sequence(seq), size(seq.size()) {}

    friend std::istream& operator>>(std::istream&, Fasta&);
};
