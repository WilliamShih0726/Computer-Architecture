#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INITIAL_CAPACITY 8
#define LOAD_FACTOR_THRESHOLD 0.7

// -------------------- Structs --------------------

typedef struct {
    int *data;
    int size;
    int capacity;
} IntArray;

typedef struct Entry {
    char *key;
    IntArray value;
    struct Entry *next;
} Entry;

typedef struct {
    Entry **buckets;
    int size;       // number of entries
    int capacity;   // number of buckets
} HashMap;

// -------------------- IntArray Functions --------------------

void init_int_array(IntArray *arr) {
    arr->size = 0;
    arr->capacity = 4;
    arr->data = malloc(arr->capacity * sizeof(int));
}

void add_to_int_array(IntArray *arr, int value) {
    if (arr->size >= arr->capacity) {
        arr->capacity += 4;
        arr->data = realloc(arr->data, arr->capacity * sizeof(int));
    }
    arr->data[arr->size++] = value;
}

// -------------------- Hash Functions --------------------

unsigned int hash(const char *str, int capacity) {
    unsigned long hash = 5381;
    int c;
    while ((c = *str++))
        hash = ((hash << 5) + hash) + c;
    return hash % capacity;
}

// -------------------- HashMap Functions --------------------

HashMap* create_hashmap() {
    HashMap *map = malloc(sizeof(HashMap));
    map->capacity = INITIAL_CAPACITY;
    map->size = 0;
    map->buckets = calloc(map->capacity, sizeof(Entry*));
    return map;
}

void free_entry_list(Entry *entry) {
    while (entry) {
        Entry *next = entry->next;
        free(entry->key);
        free(entry->value.data);
        free(entry);
        entry = next;
    }
}

void free_hashmap(HashMap *map) {
    for (int i = 0; i < map->capacity; i++) {
        free_entry_list(map->buckets[i]);
    }
    free(map->buckets);
    free(map);
}

// Rehashes all entries into a new bucket array
void resize_hashmap(HashMap *map) {
    int oldCapacity = map->capacity;
    Entry **oldBuckets = map->buckets;

    map->capacity += INITIAL_CAPACITY;
    map->buckets = calloc(map->capacity, sizeof(Entry*));
    map->size = 0;

    for (int i = 0; i < oldCapacity; i++) {
        Entry *entry = oldBuckets[i];
        while (entry) {
            Entry *next = entry->next;

            // Reinsert key-value into new table
            unsigned int index = hash(entry->key, map->capacity);
            entry->next = map->buckets[index];
            map->buckets[index] = entry;

            map->size++;
            entry = next;
        }
    }
    free(oldBuckets);
}

// Insert or update a key with a new value
void insert(HashMap *map, const char *key, int value) {
    double loadFactor = (double)(map->size + 1) / map->capacity;
    if (loadFactor > LOAD_FACTOR_THRESHOLD) {
        resize_hashmap(map);
    }

    unsigned int index = hash(key, map->capacity);
    Entry *entry = map->buckets[index];

    while (entry) {
        if (strcmp(entry->key, key) == 0) {
            add_to_int_array(&entry->value, value);
            return;
        }
        entry = entry->next;
    }

    // Key not found, create new entry
    Entry *newEntry = malloc(sizeof(Entry));
    newEntry->key = strdup(key);
    init_int_array(&newEntry->value);
    add_to_int_array(&newEntry->value, value);
    newEntry->next = map->buckets[index];
    map->buckets[index] = newEntry;
    map->size++;
}

void print_values(HashMap *map, const char *key) {
    unsigned int index = hash(key, map->capacity);
    Entry *entry = map->buckets[index];

    while (entry) {
        if (strcmp(entry->key, key) == 0) {
            printf("%s : ", key);
            for (int i = 0; i < entry->value.size; i++) {
                printf("%d ", entry->value.data[i]);
            }
            printf("\n");
            return;
        }
        entry = entry->next;
    }

    printf("Key \"%s\" not found.\n", key);
}

void get_values(HashMap *map, const char *key, IntArray *values) {
    unsigned int index = hash(key, map->capacity);
    Entry *entry = map->buckets[index];

    while (entry) {
        if (strcmp(entry->key, key) == 0) {
            for (int i = 0; i < entry->value.size; i++) {
                add_to_int_array(values, entry->value.data[i]);
            }
        }
        entry = entry->next;
    }
}

void write_hashmap_to_file(HashMap *map, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Could not open file for writing");
        return;
    }

    for (int i = 0; i < map->capacity; i++) {
        Entry *entry = map->buckets[i];
        while (entry) {
            fprintf(file, "%s :", entry->key);
            for (int j = 0; j < entry->value.size; j++) {
                fprintf(file, " %d", entry->value.data[j]);
            }
            fprintf(file, "\n");
            entry = entry->next;
        }
    }

    fclose(file);
}

void read_hashmap_from_file(HashMap *map, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Could not open file for reading");
        return;
    }

    char line[10240];

    int i;
    while (fgets(line, sizeof(line), file)) {
        char *key = strtok(line, " :\t\n");
        if (!key) continue;

        char *token = strtok(NULL, " \t\n");
	i = 0;
        while (token) {
            int value = atoi(token);
	    if (i > 0)
                insert(map, key, value);  // Reuse existing insert logic
            token = strtok(NULL, " \t\n");
	    ++i;
        }
    }

    fclose(file);
}
