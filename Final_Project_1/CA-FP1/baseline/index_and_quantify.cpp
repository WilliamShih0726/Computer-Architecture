#include <iostream>
#include <unordered_map>
#include <vector>
#include <thread>
#include <mutex>
#include <string>
//#include <algorithm>
#include <memory>
#include <shared_mutex>
#include "sequence_util.h"


// Vector to store transcript sequence wise quantification result
std::vector<int> result_counts;

// Mutex for updating quantification final results
std::mutex result_counts_mutex;

class ThreadSafeHashMap {
private:
    struct ValueWrapper {
        std::vector<int> data;
        std::mutex mutex;
    };

    std::unordered_map<std::string, std::shared_ptr<ValueWrapper>> map;
    std::mutex map_mutex;

public:
    void insert(const std::string& key, int value) {
        std::shared_ptr<ValueWrapper> wrapper;

        {
            std::lock_guard<std::mutex> lock(map_mutex);
            auto it = map.find(key);
            if (it == map.end()) {
                wrapper = std::make_shared<ValueWrapper>();
                map[key] = wrapper;
            } else {
                wrapper = it->second;
            }
        }

        {
            std::lock_guard<std::mutex> lock(wrapper->mutex);
	    bool exists = false;
	    for(int i=0; i<wrapper->data.size(); ++i){
                if(wrapper->data[i] == value)
		    exists = true;
            }
	    if(!exists)
                wrapper->data.push_back(value);
        }
    }

    std::vector<int> get_values(const std::string& key) {
        std::lock_guard<std::mutex> lock(map_mutex);
        auto it = map.find(key);
        if (it != map.end()) {
            std::lock_guard<std::mutex> lock2(it->second->mutex);
            return it->second->data;
        }
        return {};
    }

    void print_all() {
        std::lock_guard<std::mutex> lock(map_mutex);
        for (const auto& pair : map) {
            std::lock_guard<std::mutex> lock2(pair.second->mutex);
            std::cout << pair.first << ": ";
            for (int val : pair.second->data) {
                std::cout << val << " ";
            }
            std::cout << "\n";
        }
    }
};

void index_kmers(ThreadSafeHashMap& map, const std::vector<std::string> transcript_sequences, int thread_count, int thread_id, int k) {
    int sequences_per_thread = (transcript_sequences.size()+thread_count-1)/thread_count;
    int first_sequence = thread_id*sequences_per_thread;
    int last_sequence = (thread_id+1)*sequences_per_thread;
    last_sequence = last_sequence < transcript_sequences.size() ? last_sequence : transcript_sequences.size();

    for (int i = first_sequence; i < last_sequence; ++i) {
        int kmer_count = transcript_sequences[i].length()+1-k;
	for (int j = 0; j < kmer_count; ++j){
	    std::string kmer = transcript_sequences[i].substr(j, k);
	    map.insert(kmer, i);
	}
    }
}

void quantify(ThreadSafeHashMap& map, const std::vector<std::string> query_sequences, int thread_count, int thread_id, int k) {
    int sequences_per_thread = (query_sequences.size()+thread_count-1)/thread_count;
    int first_sequence = thread_id*sequences_per_thread;
    int last_sequence = (thread_id+1)*sequences_per_thread;
    last_sequence = last_sequence < query_sequences.size() ? last_sequence : query_sequences.size();

    for (int i = first_sequence; i < last_sequence; ++i) {
        int kmer_count = query_sequences[i].length()+1-k;
	std::unordered_map<int, int> transcript_wise_count_map;
	for (int j = 0; j < kmer_count; ++j){
             std::string kmer = query_sequences[i].substr(j, k);
	     std::vector<int> transcripts = map.get_values(kmer);
	     for (int l = 0; l < transcripts.size(); ++l){
	         auto it = transcript_wise_count_map.find(transcripts[l]);
		 if (it != transcript_wise_count_map.end())
                     (it->second)++;
		 else
                     transcript_wise_count_map[transcripts[l]] = 1;
	     }
        }

	int max_count_for_a_transcript = 0;
	for (const auto& pair : transcript_wise_count_map) {
            if (pair.second > max_count_for_a_transcript)
                max_count_for_a_transcript = pair.second;
        }

	std::vector<int> transcript_idxs_with_maximum_count;
	if (max_count_for_a_transcript > 0){
	    for (const auto& pair : transcript_wise_count_map) {
                if (pair.second == max_count_for_a_transcript)
                    transcript_idxs_with_maximum_count.push_back(pair.first);
            }
	}

	for (int j = 0; j < transcript_idxs_with_maximum_count.size(); ++j) {
            // Lock access to the shared final results
            std::lock_guard<std::mutex> lock(result_counts_mutex);
            result_counts[transcript_idxs_with_maximum_count[j]]++;
	}
    }
}

int main(int argc, char *argv[]) {
    // Parse the program inputs
    std::string transcript_filename = argv[1];
    std::string query_filename = argv[2];
    int k = std::stoi(argv[3]);
    int threads_count = std::stoi(argv[4]);

    // Read the reference and query sequences files
    std::vector<std::string> transcript_sequences = readFastaSequences(transcript_filename);
    std::vector<std::string> query_sequences = readFastaSequences(query_filename);

    // Shared map from k-mers to array idx of reference sequences
    ThreadSafeHashMap map;

    // Threads to be launched for generating the k-mer to reference sequence idx
    std::vector<std::thread> indexing_threads;

    // Launch threads for Indexing
    for(int i=0; i<threads_count; ++i){
        indexing_threads.emplace_back(index_kmers, std::ref(map), transcript_sequences, threads_count, i, k);
    }

    // Wait for threads to finish
    for (auto& t : indexing_threads) {
        t.join();
    }

    std::cout << "Indexing done" << std::endl;
    
    // Print the kmer to reference sequence idx map. Uncomment to observe index structure
    //map.print_all();

    // Vector to store transcript sequence wise quantification result
    result_counts.assign(transcript_sequences.size(), 0);

    // Threads to be launched for quantification
    std::vector<std::thread> quantifying_threads;

    // Launch threads for Sequence Quantification
    for(int i=0; i<threads_count; ++i){
        quantifying_threads.emplace_back(quantify, std::ref(map), query_sequences, threads_count, i, k);
    }

    // Wait for threads to finish
    for (auto& t : quantifying_threads) {
        t.join();
    }
    
    // Print final results of quantification
    for (int i = 0; i< result_counts.size(); ++i)
        std::cout << i << "\t" << result_counts[i] << std::endl;
    
    return 0;
}

