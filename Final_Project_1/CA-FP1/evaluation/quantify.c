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

int index_sequences_count;

// Sample function for generating hash-map index and correctness check during evaluation
// DO NOT CHANGE
void quantify_default(HashMap *map, const char **query_sequences_used, int query_sequences_used_count, int k, int *final_results){
    char *kmer = malloc(k+1);
    kmer[k] = '\0';

    // Keep record of transcript-wise query kmer count
    // Use maximum number of transcripts to minimize memory reallocations
    int *matched_transcripts_counts = (int*)malloc(index_sequences_count*sizeof(int));

    for (int i = 0; i < query_sequences_used_count; ++i) {
        const char *seq = query_sequences_used[i];
        int len = strlen(seq);
        int kmer_count = len+1-k;
	
	// Initialize transcript-wise query kmer counts
	for (int j = 0; j < index_sequences_count; ++j)
	    matched_transcripts_counts[j] = 0;
	
        for (int j = 0; j < kmer_count; ++j) {
	    strncpy(kmer, &seq[j], k);
	    IntArray matched_transcripts_idxs;
	    init_int_array(&matched_transcripts_idxs);
	    get_values(map, kmer, &matched_transcripts_idxs);

	    for (int l = 0; l < matched_transcripts_idxs.size; ++l) {
	        matched_transcripts_counts[matched_transcripts_idxs.data[l]]++;
	    }
	}

	// Find transcript with maximum matches, and record its match count
	int max_count_for_a_transcript = 0;
	for (int j = 0; j < index_sequences_count; ++j) {
	    if (max_count_for_a_transcript < matched_transcripts_counts[j])
		    max_count_for_a_transcript = matched_transcripts_counts[j];
	}

	// Proceed only if any of the transcripts matched atleast once for atleast one kmer
	if (max_count_for_a_transcript == 0)
	    continue;

	// Record all transcripts that matched with the query sequence for maximum number of times
	IntArray transcript_idxs_with_maximum_count;
	init_int_array(&transcript_idxs_with_maximum_count);
	for (int j = 0; j < index_sequences_count; ++j){
	    if(matched_transcripts_counts[j] == max_count_for_a_transcript){
	        add_to_int_array(&transcript_idxs_with_maximum_count, j);
	    }
	}

	// Increase count of final matched transcripts for each query sequence
	for (int j = 0; j < transcript_idxs_with_maximum_count.size; ++j){
	    final_results[transcript_idxs_with_maximum_count.data[j]]++;
	}
    }
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
typedef struct {
    uint64_t   key;      /* 2-bit 編碼 k-mer；0 = Empty          */
    IntArray  *plist;    /* 指向原 HashMap 內的 transcripts 陣列 */
} FastBkt;

void quantify(HashMap *map, const char **query_sequences_used, int query_sequences_used_count, int k, int *final_results){
    /* 1-1 估算 map 的總 bucket 數，設定快速表容量為最近 2 次方      */
    size_t cap = 1;
    while (cap < (size_t)(map->capacity * 2)) cap <<= 1;  /* load < 0.5 */
    const size_t mask = cap - 1;

    FastBkt *ftab = (FastBkt *)calloc(cap, sizeof(FastBkt));

    /* 1-2 把原 HashMap 內所有 Entry 轉成 64-bit key 存進 ftab      */
    static const uint8_t nt2b[256] = { ['A']=0,['C']=1,['G']=2,['T']=3,
                                       ['a']=0,['c']=1,['g']=2,['t']=3 };
    for (int b = 0; b < map->capacity; ++b) {
        Entry *e = map->buckets[b];
        while (e) {
            /* 把字串 key → 64-bit 2-bit 編碼 */
            uint64_t key64 = 0;
            for (int i = 0; i < k; ++i)
                key64 = (key64 << 2) | nt2b[(uint8_t)e->key[i]];

            size_t pos = key64 & mask;
            while (ftab[pos].key && ftab[pos].key != key64)
                pos = (pos + 1) & mask;
            ftab[pos].key   = key64;
            ftab[pos].plist = &e->value;
            e = e->next;
        }
    }

    int *cnt   = (int *)calloc(index_sequences_count, sizeof(int));
    int *touched = (int *)malloc(k * 2 * sizeof(int));
    int  touched_sz;

    const uint64_t k_mask = ((uint64_t)1 << (2 * k)) - 1;

    /* ---------- 逐條 query 進行量化 ---------- */
    for (int qi = 0; qi < query_sequences_used_count; ++qi) {
        const char *s = query_sequences_used[qi];
        int L = (int)strlen(s);
        if (L < k) continue;

        /* 3-1 Rolling key for first k-mer */
        uint64_t key = 0;
        for (int i = 0; i < k; ++i)
            key = (key << 2) | nt2b[(uint8_t)s[i]];

        touched_sz = 0;

        /* lambda-like內聯函式：查表並更新計數 ----------------- */
        #define PROCESS_KEY(K)                                        \
        do {                                                          \
            size_t p = (K) & mask;                                    \
            while (ftab[p].key && ftab[p].key != (K))                 \
                p = (p + 1) & mask;                                   \
            if (ftab[p].key) {                                        \
                IntArray *arr = ftab[p].plist;                        \
                for (int z = 0; z < arr->size; ++z) {                 \
                    int tid = arr->data[z];                           \
                    if (cnt[tid] == 0) touched[touched_sz++] = tid;   \
                    cnt[tid]++;                                       \
                }                                                     \
            }                                                         \
        } while(0)

        PROCESS_KEY(key);

        for (int j = k; j < L; ++j) { 
            key = ((key << 2) & k_mask) | nt2b[(uint8_t)s[j]];
            PROCESS_KEY(key);
        }

        int maxc = 0;
        for (int idx = 0; idx < touched_sz; ++idx) {
            int tid = touched[idx];
            if (cnt[tid] > maxc) maxc = cnt[tid];
        }
        if (maxc == 0) {
            for (int idx = 0; idx < touched_sz; ++idx)
                cnt[touched[idx]] = 0;
            continue;
        }

        for (int idx = 0; idx < touched_sz; ++idx) {
            int tid = touched[idx];
            if (cnt[tid] == maxc) final_results[tid]++;
            cnt[tid] = 0;
        }
    }

    free(ftab);
    free(cnt);
    free(touched);
}

int main(int argc, char *argv[]){
    int k = atoi(argv[1]);
    // Select one from transcript_sequences and query_sequences to use as query_sequences_used. The other one should be treated as index_sequences accordingly.
    const char **index_sequences = transcript_sequences;
    index_sequences_count = transcript_sequences_count;
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

    // Uncomment to print final quantification results
    //for(int i =0; i < index_sequences_count; ++i)
    //    printf("Transcript %d has %d matched queries\n", i, final_results[i]);

    free_hashmap(map);
    free(final_results);
    //free(final_results_default);

    return 0;
}
