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

// Writes a vector of strings to a header file as a const std::vector<std::string>
void writeVectorToHeader(const std::vector<std::string> vec, const std::string fileName, const std::string varName) {
    std::ofstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << fileName << " for writing.\n";
        return;
    }

    std::cout << "Writing file " << fileName << std::endl;
    file << "int " << varName << "_count = " << vec.size() << ";\n";
    file << "const char *" << varName << "[] = {\n";

    for (const auto str : vec) {
        file << "    \"" << str << "\",\n";
    }

    file << "};\n";
    file.close();

    std::cout << "✅ Wrote to " << fileName << " successfully.\n";
}

// Wrapper function to export two vectors to two .h files
void exportTwoVectorsToHeaders(const std::vector<std::string> vec1,
                               const std::vector<std::string> vec2,
                               const std::string file1,
                               const std::string file2,
                               const std::string varName1,
                               const std::string varName2) {
    writeVectorToHeader(vec1, file1, varName1);
    writeVectorToHeader(vec2, file2, varName2);
}

