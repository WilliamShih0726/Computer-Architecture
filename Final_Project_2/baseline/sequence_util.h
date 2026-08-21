#include <iostream>
#include <fstream>
#include <string>
#include <vector>

std::vector<std::string> readFastaSequences(const std::string& filename) {
    std::ifstream infile(filename);
    std::vector<std::string> sequences;

    if (!infile) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return sequences;
    }

    std::string line;
    std::string current_sequence;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current_sequence.empty()) {
                sequences.push_back(current_sequence);
                current_sequence.clear();
            }
            // skip storing header, just move to next
        } else {
            current_sequence += line;
        }
    }

    // Push last sequence if file doesn't end with a header
    if (!current_sequence.empty()) {
        sequences.push_back(current_sequence);
    }

    infile.close();
    return sequences;
}

