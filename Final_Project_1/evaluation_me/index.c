#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>
#include <signal.h>
#include "/home/gem5/include/gem5/m5ops.h"
#include "transcripts_data.h"
#include "query_data.h"
#include "hash-map.h"
#include "gem5_utils.h"

// Sample function for generating hash-map index and correctness check during evaluation
// DO NOT CHANGE
void generate_index_default(const char **index_sequences, int index_sequences_count, int k, HashMap *map){
    for (int i = 0; i < index_sequences_count; i++) {
        const char *seq = index_sequences[i];
        int len = strlen(seq);

        // Skip if sequence is too short for k-mers
        if (len < k) continue;

	
        for (int j = 0; j <= len - k; j++) {
            char *kmer = malloc(k+1);
            strncpy(kmer, &seq[j], k);
            kmer[k] = '\0';

            // Insert into hash map
            insert(map, kmer, i);
        }
    }
}

// Function to be evaluated. Implement this.
// See definition of HashMap and its operations in hash-map.h
// Use your own algorithm to generate the HashMap, functions in hash-map.h are only for your reference and correctness check during evaluations
// DO NOT CHANGE function arguments
void generate_index(const char **index_sequences, int index_sequences_count, int k, HashMap *map){

}

int main(int argc, char *argv[]){
    int k = atoi(argv[1]);
    // Select one from transcript_sequences and query_sequences to generate index. The other one should be treated as query sequences accordingly.
    const char **index_sequences = transcript_sequences;
    int index_sequences_count = transcript_sequences_count;
    // const char **index_sequences = query_sequences;
    // int index_sequences_count = query_sequences_count;

    // Final hash-map to be evaluated
    HashMap *map = create_hashmap();

    // Reset Gem5 simulation statistics
    m5_reset_stats(0,0);

    // Function to implement
    // Call index function to index index_sequences within a single Logic Unit
    // You are free to use Vector Coprocessors (SIMD), which will fetch higher evaluation score if it improves overall runtime
    // Goal of this function is to improve the performance of Indexing by a single Logic Unit
    generate_index(index_sequences, index_sequences_count, k, map);

    // Dump Gem5 simulation final statistics
    m5_dump_stats(0,0);

    // Print Gem5 simulation time for the evaluated indexing function
    print_gem5_simulated_time("m5out/stats.txt");

    // Write the generated index (hash-map) to a file, to be later read for quantification in quantify.c
    write_hashmap_to_file(map, "index_result/hash-map.txt");

    // Function for correctness check. Keep commented to avoid long Gem5 simulation time
    // DO NOT CHANGE
    // HashMap *map_default = create_hashmap();
    // generate_index_default(index_sequences, index_sequences_count, k, map_default);
    // check_correctness(map_default, map);

    free_hashmap(map);
    //free(map_default);

    return 0;
}
