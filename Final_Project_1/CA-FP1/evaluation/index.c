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
    char *kmer = malloc(k+1);
    kmer[k] = '\0';
    for (int i = 0; i < index_sequences_count; i++) {
        const char *seq = index_sequences[i];
        int len = strlen(seq);

        // Skip if sequence is too short for k-mers
        if (len < k) continue;
	
        for (int j = 0; j <= len - k; j++) {
            strncpy(kmer, &seq[j], k);

            // Insert into hash map
            insert(map, kmer, i);
        }
    }
    free(kmer);
}

// Function to be evaluated. Implement this.
// See definition of HashMap and its operations in hash-map.h
// Use your own algorithm to generate the HashMap, functions in hash-map.h are only for your reference and correctness check during evaluations
// DO NOT CHANGE function arguments
typedef struct {
    uint64_t key;     /* 2-bit 編碼後的 k-mer；0 表空 bucket */
    int     *ids;     /* 動態陣列存 transcript IDs */
    int      sz;      /* 已存元素數 */
    int      cap;     /* 陣列容量 */
} Bucket;

static inline void push_id(Bucket *b, int tid)
{
    if (b->cap == 0) {
        b->cap = 4;
        b->ids = (int *)malloc(4 * sizeof(int));
    } else if (b->sz == b->cap) {
        b->cap <<= 1;
        b->ids = (int *)realloc(b->ids, b->cap * sizeof(int));
    }
    b->ids[b->sz++] = tid;
}

void generate_index(const char **index_sequences, int index_sequences_count, int k, HashMap *map){
    /* ---------- 1. 轉碼表 & 掩碼 ---------- */
    static const uint8_t nt2b[256] = {
        ['A']=0,['C']=1,['G']=2,['T']=3,
        ['a']=0,['c']=1,['g']=2,['t']=3
    };
    static const char LUT[4] = {'A','C','G','T'};
    const uint64_t k_mask = ((uint64_t)1 << (2 * k)) - 1;

    /* ---------- 2. 粗估 k-mer 數，決定 bucket 容量 ---------- */
    size_t total = 0;
    for (int i = 0; i < index_sequences_count; ++i) {
        int L = (int)strlen(index_sequences[i]);
        if (L >= k) total += (size_t)(L - k + 1);
    }
    size_t cap = 1;
    while (cap < (size_t)(total / 0.7)) cap <<= 1;   /* 目標載入率 ≤0.7 */
    const size_t mask = cap - 1;

    /* ---------- 3. 配置暫存表 ---------- */
    Bucket *tab = (Bucket *)calloc(cap, sizeof(Bucket));

    /* ---------- 4. 掃描所有 transcripts ---------- */
    for (int tid = 0; tid < index_sequences_count; ++tid) {
        const char *s = index_sequences[tid];
        int L = (int)strlen(s);
        if (L < k) continue;

        /* 4-1. 初始 k-mer key */
        uint64_t key = 0;
        for (int i = 0; i < k; ++i)
            key = (key << 2) | nt2b[(uint8_t)s[i]];

        /* 4-2. 插入第 0 個 k-mer */
        size_t pos = key & mask;
        while (tab[pos].key && tab[pos].key != key) pos = (pos + 1) & mask;
        if (tab[pos].key == 0) tab[pos].key = key;
        push_id(&tab[pos], tid);

        /* 4-3. 滑窗處理其餘 k-mers */
        for (int j = k; j < L; ++j) {
            key = ((key << 2) & k_mask) | nt2b[(uint8_t)s[j]];
            pos = key & mask;
            while (tab[pos].key && tab[pos].key != key) pos = (pos + 1) & mask;
            if (tab[pos].key == 0) tab[pos].key = key;
            push_id(&tab[pos], tid);
        }
    }

    /* ---------- 5. 批量寫回原 HashMap ---------- */
    char kmer[64];
    kmer[k] = '\0';

    for (size_t i = 0; i < cap; ++i) {
        if (tab[i].key == 0) continue;

        /* key → 字串 */
        uint64_t ky = tab[i].key;
        for (int p = k - 1; p >= 0; --p) {
            kmer[p] = LUT[ky & 3];
            ky >>= 2;
        }

        /* 第一個 id 建 bucket，其餘追加 */
        insert(map, kmer, tab[i].ids[0]);
        for (int t = 1; t < tab[i].sz; ++t)
            insert(map, kmer, tab[i].ids[t]);

        free(tab[i].ids);   /* 釋放動態陣列 */
    }
    free(tab);              /* 釋放暫存表 */
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