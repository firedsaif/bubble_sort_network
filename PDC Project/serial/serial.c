#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define MAX_PATH_LENGTH 10000
#define OUTPUT_BUFFER_SIZE (10*1024*1024) // Reduced to 10MB
#define MAX_VERTICES_TO_OUTPUT 100000000 // Limit output for large N

typedef struct {
    char **keys;
    int *values;
    int size;
    int capacity;
} HashMap;

HashMap* create_hashmap(int capacity);
void hashmap_put(HashMap* map, char* key, int value);
int hashmap_get(HashMap* map, char* key, int default_value);
void free_hashmap(HashMap* map);

int factorial(int n);
void swap_adjacent(char *perm, int i, int n);
void get_permutation(int idx, int n, char *out);
int find_position(char *perm, int n);
int index_of(char *perm, int n, char symbol);
void Swap(char *v, int n, char x, char *out);
void FindPosition(char *v, int t, int n, char *out);
void Parent1(char *v, int t, int n, char *parent, char *root);

int factorial(int n) {
    int res = 1;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}

void swap_adjacent(char *perm, int i, int n) {
    char tmp = perm[i];
    perm[i] = perm[i + 1];
    perm[i + 1] = tmp;
}

void get_permutation(int idx, int n, char *out) {
    int *nums = malloc(n * sizeof(int));
    if (!nums) { fprintf(stderr, "Memory allocation failed in get_permutation\n"); exit(1); }
    for (int i = 0; i < n; ++i) nums[i] = i + 1;
    for (int i = 0; i < n; ++i) {
        int fact = factorial(n - 1 - i);
        int selected = idx / fact;
        out[i] = '0' + nums[selected];
        for (int j = selected; j < n - 1 - i; ++j) nums[j] = nums[j + 1];
        idx = idx % fact;
    }
    out[n] = '\0';
    free(nums);
}

int find_position(char *perm, int n) {
    for (int i = n - 1; i >= 0; --i) {
        if (perm[i] != '1' + i) return i + 1;
    }
    return -1;
}

int index_of(char *perm, int n, char symbol) {
    for (int i = 0; i < n; ++i)
        if (perm[i] == symbol) return i;
    return -1;
}

void Swap(char *v, int n, char x, char *out) {
    strcpy(out, v);
    int i = index_of(out, n, x);
    if (i < n - 1) swap_adjacent(out, i, n);
}

void FindPosition(char *v, int t, int n, char *out) {
    char *root = malloc((n + 1) * sizeof(char));
    if (!root) { fprintf(stderr, "Memory allocation failed in FindPosition\n"); exit(1); }
    for (int i = 0; i < n; ++i) root[i] = '1' + i;
    root[n] = '\0';

    char *tmp = malloc((n + 1) * sizeof(char));
    if (!tmp) { fprintf(stderr, "Memory allocation failed in FindPosition\n"); exit(1); }
    Swap(v, n, '0' + t, tmp);
    if (t == 2 && strcmp(tmp, root) == 0) {
        Swap(v, n, '0' + (t - 1), out);
    } else {
        int vn_1 = v[n - 2] - '0';
        if (vn_1 == t || vn_1 == n - 1) {
            int j = find_position(v, n);
            if (j != -1) {
                Swap(v, n, '0' + j, out);
            } else {
                strcpy(out, v);
            }
        } else {
            Swap(v, n, '0' + t, out);
        }
    }
    free(tmp);
    free(root);
}

void Parent1(char *v, int t, int n, char *parent, char *root) {
    if (strcmp(v, root) == 0) {
        strcpy(parent, "ROOT");
        return;
    }

    int vn = v[n - 1] - '0';
    char *tmp = malloc((n + 1) * sizeof(char));
    if (!tmp) { fprintf(stderr, "Memory allocation failed in Parent1\n"); exit(1); }

    if (vn == n) {
        if (t != n - 1) {
            FindPosition(v, t, n, parent);
        } else {
            int vn_1 = v[n - 2] - '0';
            Swap(v, n, '0' + vn_1, parent);
        }
        free(tmp);
        return;
    }

    if (vn == n - 1) {
        int vn_1 = v[n - 2] - '0';
        Swap(v, n, '0' + n, tmp);
        if (vn_1 != n || strcmp(tmp, root) == 0) {
            if (vn == t) {
                Swap(v, n, '0' + n, parent);
            } else {
                Swap(v, n, '0' + t, parent);
            }
        } else {
            if (t == 1) {
                Swap(v, n, '0' + n, parent);
            } else {
                Swap(v, n, '0' + (t - 1), parent);
            }
        }
        free(tmp);
        return;
    }

    if (vn <= n - 2) {
        if (vn == t) {
            Swap(v, n, '0' + n, parent);
        } else {
            Swap(v, n, '0' + t, parent);
        }
        free(tmp);
        return;
    }
    free(tmp);
}

const char *colors[] = {"red", "blue", "green", "orange", "purple", "brown", "cyan", "magenta", "yellow", "black"};

HashMap* create_hashmap(int capacity) {
    HashMap* map = malloc(sizeof(HashMap));
    if (!map) { fprintf(stderr, "Memory allocation failed in create_hashmap\n"); exit(1); }
    map->capacity = capacity * 2;
    map->size = 0;
    map->keys = malloc(map->capacity * sizeof(char*));
    map->values = malloc(map->capacity * sizeof(int));
    if (!map->keys || !map->values) { fprintf(stderr, "Memory allocation failed in create_hashmap\n"); exit(1); }
    for (int i = 0; i < map->capacity; i++) map->keys[i] = NULL;
    return map;
}

void hashmap_put(HashMap* map, char* key, int value) {
    unsigned long hash = 5381;
    for (int i = 0; key[i]; i++) hash = ((hash << 5) + hash) + key[i];
    int index = hash % map->capacity;
    while (map->keys[index] != NULL) {
        if (strcmp(map->keys[index], key) == 0) {
            map->values[index] = value;
            return;
        }
        index = (index + 1) % map->capacity;
    }
    map->keys[index] = strdup(key);
    if (!map->keys[index]) { fprintf(stderr, "Memory allocation failed in hashmap_put\n"); exit(1); }
    map->values[index] = value;
    map->size++;
}

int hashmap_get(HashMap* map, char* key, int default_value) {
    unsigned long hash = 5381;
    for (int i = 0; key[i]; i++) hash = ((hash << 5) + hash) + key[i];
    int index = hash % map->capacity;
    while (map->keys[index] != NULL) {
        if (strcmp(map->keys[index], key) == 0) return map->values[index];
        index = (index + 1) % map->capacity;
    }
    return default_value;
}

void free_hashmap(HashMap* map) {
    for (int i = 0; i < map->capacity; i++) if (map->keys[i]) free(map->keys[i]);
    free(map->keys);
    free(map->values);
    free(map);
}

int main() {
    int n = 10; // Test with N=10, works up to N=11
    int total = factorial(n);
    if (total > 50000000) {
        printf("N=%d results in %d vertices, which exceeds practical limits. Aborting.\n", n, total);
        return 1;
    }

    printf("Constructing spanning trees for bubble-sort network with N=%d (%d vertices)\n", n, total);

    // Allocate and generate all permutations
    char **vertices = malloc(total * sizeof(char*));
    if (!vertices) { fprintf(stderr, "Memory allocation failed for vertices\n"); return 1; }
    for (int i = 0; i < total; ++i) {
        vertices[i] = malloc((n + 1) * sizeof(char));
        if (!vertices[i]) { fprintf(stderr, "Memory allocation failed for vertex %d\n", i); return 1; }
        get_permutation(i, n, vertices[i]);
    }

    // Create hash map for permutation lookup
    HashMap* vertex_map = create_hashmap(total);
    for (int i = 0; i < total; i++) hashmap_put(vertex_map, vertices[i], i);

    // Allocate graph structures
    int *xadj = malloc((total + 1) * sizeof(int));
    int *adj = malloc(total * (n - 1) * sizeof(int));
    if (!xadj || !adj) { fprintf(stderr, "Memory allocation failed for graph structures\n"); return 1; }

    // Build adjacency list sequentially
    printf("Building graph adjacency list...\n");
    xadj[0] = 0;
    int edge_idx = 0;
    for (int i = 0; i < total; i++) {
        char *perm1 = vertices[i];
        int local_edge_count = 0;

        // Try all possible adjacent swaps
        for (int pos = 0; pos < n - 1; pos++) {
            char *perm2 = malloc((n + 1) * sizeof(char));
            if (!perm2) { fprintf(stderr, "Memory allocation failed for perm2\n"); return 1; }
            strcpy(perm2, perm1);
            swap_adjacent(perm2, pos, n);
            int j = hashmap_get(vertex_map, perm2, -1);
            if (j != -1 && j != i) {
                adj[edge_idx++] = j;
                local_edge_count++;
            }
            free(perm2);
        }
        xadj[i + 1] = xadj[i] + local_edge_count;

        //if (i % 10000 == 0) {
        //    printf("Progress: Processed %d/%d vertices for adjacency list\n", i, total);
        //}
    }

    // Allocate parents array
    int *parents = malloc((n - 1) * total * sizeof(int));
    if (!parents) { fprintf(stderr, "Memory allocation failed for parents\n"); return 1; }
    for (int i = 0; i < (n - 1) * total; ++i) parents[i] = -1;

    // Identify root permutation
    char *root = malloc((n + 1) * sizeof(char));
    if (!root) { fprintf(stderr, "Memory allocation failed for root\n"); return 1; }
    for (int i = 0; i < n; ++i) root[i] = '1' + i;
    root[n] = '\0';
    int root_idx = hashmap_get(vertex_map, root, -1);
    if (root_idx == -1) { fprintf(stderr, "Root permutation not found in hash map\n"); return 1; }

    // Compute parent relationships sequentially
    printf("Computing parent relationships...\n");
    char *local_parent = malloc((n + 1) * sizeof(char));
    if (!local_parent) { fprintf(stderr, "Memory allocation failed for local_parent\n"); return 1; }

    for (int t = 0; t < n - 1; t++) {
        for (int i = 0; i < total; i++) {
            if (i == root_idx) continue;
            Parent1(vertices[i], t + 1, n, local_parent, root);
            int parent_idx = (strcmp(local_parent, "ROOT") == 0) ? root_idx : hashmap_get(vertex_map, local_parent, -1);
            if (parent_idx != -1 && parent_idx >= 0 && parent_idx < total) {
                parents[t * total + i] = parent_idx;
            } else if (parent_idx != -1) {
                fprintf(stderr, "Invalid parent_idx %d for vertex %d in tree T_%d\n", parent_idx, i, t + 1);
            }
        }
        printf("Processed tree T_%d\n", t + 1);
    }
    free(local_parent);

    // Write output to file
    FILE *output_file = fopen("spanning_trees_output.txt", "w");
    if (!output_file) { fprintf(stderr, "Error opening output file\n"); return 1; }

    int output_limit = (total < MAX_VERTICES_TO_OUTPUT) ? total : MAX_VERTICES_TO_OUTPUT;
    char *buffer = malloc(OUTPUT_BUFFER_SIZE);
    if (!buffer) { fprintf(stderr, "Memory allocation failed for output buffer\n"); return 1; }
    size_t buffer_pos = 0;

    for (int t = 0; t < n - 1; t++) {
        char temp[1024];
        snprintf(temp, sizeof(temp), "\nTree T_%d (color: %s):\n", t + 1, colors[t]);
        size_t len = strlen(temp);
        if (buffer_pos + len + 1 < OUTPUT_BUFFER_SIZE) {
            strcpy(buffer + buffer_pos, temp);
            buffer_pos += len;
        } else {
            if (fwrite(buffer, 1, buffer_pos, output_file) != buffer_pos) {
                fprintf(stderr, "Error writing to output file\n"); return 1;
            }
            buffer_pos = 0;
            strcpy(buffer, temp);
            buffer_pos += len;
        }

        snprintf(temp, sizeof(temp), "Parents:\n");
        len = strlen(temp);
        if (buffer_pos + len + 1 < OUTPUT_BUFFER_SIZE) {
            strcpy(buffer + buffer_pos, temp);
            buffer_pos += len;
        } else {
            if (fwrite(buffer, 1, buffer_pos, output_file) != buffer_pos) {
                fprintf(stderr, "Error writing to output file\n"); return 1;
            }
            buffer_pos = 0;
            strcpy(buffer, temp);
            buffer_pos += len;
        }

        int vertex_count = 0;
        for (int i = 0; i < total; i++) {
            if (i != root_idx && parents[t * total + i] != -1) vertex_count++;
        }
        snprintf(temp, sizeof(temp), "Total vertices with parents: %d\n", vertex_count);
        len = strlen(temp);
        if (buffer_pos + len + 1 < OUTPUT_BUFFER_SIZE) {
            strcpy(buffer + buffer_pos, temp);
            buffer_pos += len;
        } else {
            if (fwrite(buffer, 1, buffer_pos, output_file) != buffer_pos) {
                fprintf(stderr, "Error writing to output file\n"); return 1;
            }
            buffer_pos = 0;
            strcpy(buffer, temp);
            buffer_pos += len;
        }

        int output_count = 0;
        for (int i = 0; i < total && output_count < output_limit; i++) {
            if (i == root_idx || parents[t * total + i] == -1) continue;
            snprintf(temp, sizeof(temp), "  Vertex %s -> Parent %s\n", vertices[i],
                     parents[t * total + i] == root_idx ? "ROOT" : vertices[parents[t * total + i]]);
            len = strlen(temp);
            if (buffer_pos + len + 1 < OUTPUT_BUFFER_SIZE) {
                strcpy(buffer + buffer_pos, temp);
                buffer_pos += len;
                output_count++;
            } else {
                if (fwrite(buffer, 1, buffer_pos, output_file) != buffer_pos) {
                    fprintf(stderr, "Error writing to output file\n"); return 1;
                }
                buffer_pos = 0;
                strcpy(buffer, temp);
                buffer_pos += len;
                output_count++;
            }
        }

        if (output_count >= output_limit) {
            snprintf(temp, sizeof(temp), "[Output truncated at %d vertices due to large N]\n", output_limit);
            len = strlen(temp);
            if (buffer_pos + len + 1 < OUTPUT_BUFFER_SIZE) {
                strcpy(buffer + buffer_pos, temp);
                buffer_pos += len;
            } else {
                if (fwrite(buffer, 1, buffer_pos, output_file) != buffer_pos) {
                    fprintf(stderr, "Error writing to output file\n"); return 1;
                }
                buffer_pos = 0;
                strcpy(buffer, temp);
                buffer_pos += len;
            }
        }
    }

    if (buffer_pos > 0) {
        if (fwrite(buffer, 1, buffer_pos, output_file) != buffer_pos) {
            fprintf(stderr, "Error writing to output file\n"); return 1;
        }
    }
    if (fclose(output_file) != 0) {
        fprintf(stderr, "Error closing output file\n"); return 1;
    }
    free(buffer);
    printf("Output written to spanning_trees_output.txt\n");

    // Clean up
    free_hashmap(vertex_map);
    for (int i = 0; i < total; i++) free(vertices[i]);
    free(vertices);
    free(parents);
    free(xadj);
    free(adj);
    free(root);
    return 0;
}
