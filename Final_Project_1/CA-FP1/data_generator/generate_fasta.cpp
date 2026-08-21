#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <ctime>

// Helper to clean and concatenate sequences
std::vector<std::string> parseFasta(const std::string& filename) {
    std::ifstream infile(filename);
    std::vector<std::string> sequences;
    std::string line, currentSeq;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!currentSeq.empty()) {
                sequences.push_back(currentSeq);
                currentSeq.clear();
            }
        } else {
            currentSeq += line;
        }
    }
    if (!currentSeq.empty()) {
        sequences.push_back(currentSeq);
    }
    return sequences;
}

std::vector<std::string> extract19mers(const std::vector<std::string>& sequences) {
    std::vector<std::string> kmers;
    const int k = 19;

    for (const auto& seq : sequences) {
        if (seq.length() < k) continue;
        for (size_t i = 0; i <= seq.length() - k; ++i) {
            kmers.push_back(seq.substr(i, k));
        }
    }
    return kmers;
}

std::string buildSequence(const std::vector<std::string>& kmers, int length, std::mt19937& rng) {
    std::uniform_int_distribution<> dist(0, kmers.size() - 1);
    std::string result;
    while (result.size() < static_cast<size_t>(length)) {
        const std::string& kmer = kmers[dist(rng)];
        result += kmer;
    }
    return result.substr(0, length);
}

int main(int argc, char *argv[]) {
    const std::string inputFile = argv[1];
    const std::string outputFile = argv[2];
    const int numSequences = std::stoi(argv[3]);
    const int sequenceLength = std::stoi(argv[4]);

    std::vector<std::string> sequences = parseFasta(inputFile);
    std::vector<std::string> kmers = extract19mers(sequences);

    if (kmers.empty()) {
        std::cerr << "No 19-mers extracted. Check input file.\n";
        return 1;
    }

    std::ofstream outfile(outputFile);
    std::mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));

    for (int i = 1; i <= numSequences; ++i) {
        std::string seq = buildSequence(kmers, sequenceLength, rng);
        outfile << ">seq_" << i << "\n";

        // Wrap lines at 60 characters
        for (size_t j = 0; j < seq.size(); j += 60) {
            outfile << seq.substr(j, 60) << "\n";
        }

        if (i % 100000 == 0) {
            std::cerr << "Generated " << i << " sequences...\n";
        }
    }

    outfile.close();
    std::cout << "Finished writing " << numSequences << " sequences to " << outputFile << "\n";
    return 0;
}

