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
// Not implemented yet, will be shared ASAP
// DO NOT CHANGE
void quantify_default(HashMap *map, const char **query_sequences_used, int query_sequences_used_count, int k, int *final_results){

}

void check_correctness(int *final_results_default, int *final_results, int count){
    for (int i = 0; i < count; ++i) {
        if ( final_results[i] != final_results_default[i]) {
            printf("Failed\n");
	    return;
	}
    }
    printf("Passed");
}

// Function to be evaluated. Implement this.
// See definition of HashMap and its operations in hash-map.h
// Use your own algorithm to quantify
// DO NOT CHANGE function arguments
void quantify(HashMap *map, const char **query_sequences_used, int query_sequences_used_count, int k, int *final_results){

}

int main(int argc, char *argv[]){
    int k = atoi(argv[1]);
    // Select one from transcript_sequences and query_sequences to use as query_sequences_used. The other one should be treated as index_sequences accordingly.
    const char **index_sequences = transcript_sequences;
    int index_sequences_count = transcript_sequences_count;
    const char **query_sequences_used = query_sequences;
    int query_sequences_used_count = query_sequences_count;

    // Final quantification results to be evaluated
    int *final_results = (int*)malloc(index_sequences_count*sizeof(int));
    for(int i = 0; i < index_sequences_count; ++i)
        final_results[i] = 0;

    // Read the index (hash-map) saved in index.c
    HashMap *map = create_hashmap();
    read_hashmap_from_file(map, "index_result/hash-map.txt");

    // Reset Gem5 simulation statistics
    m5_reset_stats(0,0);

    // Function to implement
    // Call index function to index index_sequences within a single Logic Unit
    // You are free to use Vector Coprocessors (SIMD), which will fetch higher evaluation score if it improves overall runtime
    // Goal of this function is to improve the performance of Indexing by a single Logic Unit
    quantify(map, query_sequences_used, query_sequences_used_count, k, final_results);

    // Dump Gem5 simulation final statistics
    m5_dump_stats(0,0);

    // Print Gem5 simulation time for the evaluated indexing function
    print_gem5_simulated_time("m5out/stats.txt");

    // Function for correctness check. Keep commented to avoid long Gem5 simulation time
    // DO NOT CHANGE
    // quantify_default(map, query_sequences_used, query_sequences_used_count, k, final_results_default);
    // check_correctness(final_results_default, final_results);

    free_hashmap(map);
    free(final_results);
    //free(final_results_default);

    return 0;
}
