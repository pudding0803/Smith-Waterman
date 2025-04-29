#include <iostream>
#include "fasta.hpp"

std::istream& operator>>(std::istream& is, Fasta& fasta) {
    std::string line;

    if (is.peek() != '>' || !std::getline(is, line)) {
        return is.setstate(std::ios::failbit), is;
    }

    fasta.description = line.substr(1);
    fasta.sequence.clear();

    while (is.peek() != '>' && std::getline(is, line)) {
        fasta.sequence += line;
    }

    fasta.size = fasta.sequence.size();
    return is;
}