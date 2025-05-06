#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <stdbool.h>

#define IDXTYPEWIDTH 64
#define REALTYPEWIDTH 64
#include <metis.h>

#define MAX_PATH_LENGTH 10000
#define OUTPUT_BUFFER_SIZE (100*1024 * 1024) // 100MB per thread
#define MAX_VERTICES_TO_OUTPUT 10000000 // Limit output for large N

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
void partition_graph(int total, int n, idx_t *xadj, idx_t *adj, idx_t *part, char **vertices, HashMap* vertex_map);

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
    for (int i = 0; i < n; ++i) root[i] = '1' + i;
    root[n] = '\0';

    char *tmp = malloc((n + 1) * sizeof(char));
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

void partition_graph(int total, int n, idx_t *xadj, idx_t *adj, idx_t *part, char **vertices, HashMap* vertex_map) {
    idx_t num_vertices = total;
    idx_t ncon = 1;
    idx_t nparts = omp_get_max_threads();
    if (nparts > 4) nparts = 4;
    idx_t edgecut;
    idx_t options[METIS_NOPTIONS];

    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

    idx_t **adj_lists = malloc(total * sizeof(idx_t*));
    idx_t *adj_sizes = malloc(total * sizeof(idx_t));

    for (idx_t i = 0; i < total; i++) {
        char *perm1 = vertices[i];
        idx_t *local_adj = malloc((n - 1) * sizeof(idx_t));
        idx_t local_edge_idx = 0;

        for (int pos = 0; pos < n - 1; pos++) {
            char *perm2 = malloc((n + 1) * sizeof(char));
            strcpy(perm2, perm1);
            swap_adjacent(perm2, pos, n);
            int j = hashmap_get(vertex_map, perm2, -1);
            if (j != -1 && j != i) {
                local_adj[local_edge_idx++] = j;
            }
            free(perm2);
        }

        adj_lists[i] = local_adj;
        adj_sizes[i] = local_edge_idx;

        // if (i % 10000 == 0) {
        //     #pragma omp critical
        //     printf("Progress: Processed %lld/%d vertices for adjacency list\n", (long long)i, total);
        // }
    }

    xadj[0] = 0;
    for (idx_t i = 0; i < total; i++) {
        xadj[i + 1] = xadj[i] + adj_sizes[i];
    }

    for (idx_t i = 0; i < total; i++) {
        idx_t start = xadj[i];
        for (idx_t k = 0; k < adj_sizes[i]; k++) {
            adj[start + k] = adj_lists[i][k];
        }
        free(adj_lists[i]);
    }
    free(adj_lists);
    free(adj_sizes);

    int ret = METIS_PartGraphKway(&num_vertices, &ncon, xadj, adj, NULL, NULL, NULL,
                                  &nparts, NULL, NULL, options, &edgecut, part);
    if (ret != METIS_OK) {
        printf("METIS_PartGraphKway failed with error code %d\n", ret);
        exit(1);
    }
    printf("METIS partitioning successful. Edgecut: %d\n", (int)edgecut);
    int *partition_sizes = calloc(nparts, sizeof(int));
    for (int i = 0; i < total; i++) partition_sizes[part[i]]++;
    printf("Partition sizes:\n");
    for (int i = 0; i < nparts; i++) printf("  Partition %d: %d vertices\n", i, partition_sizes[i]);
    free(partition_sizes);
}

HashMap* create_hashmap(int capacity) {
    HashMap* map = malloc(sizeof(HashMap));
    map->capacity = capacity * 2;
    map->size = 0;
    map->keys = malloc(map->capacity * sizeof(char*));
    map->values = malloc(map->capacity * sizeof(int));
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
    int n = 10; // Test with N=9, should work up to N=11
    int total = factorial(n);
    if (total > 50000000) {
        printf("N=%d results in %d vertices, which exceeds practical limits. Aborting.\n", n, total);
        return 1;
    }

    printf("Constructing spanning trees for bubble-sort network with N=%d (%d vertices)\n", n, total);
    
    char **vertices = malloc(total * sizeof(char*));
    if (!vertices) { printf("Memory allocation failed for vertices\n"); return 1; }
    for (int i = 0; i < total; ++i) {
        vertices[i] = malloc((n + 1) * sizeof(char));
        if (!vertices[i]) { printf("Memory allocation failed for vertex %d\n", i); return 1; }
        get_permutation(i, n, vertices[i]);
    }

    HashMap* vertex_map = create_hashmap(total);
    for (int i = 0; i < total; i++) hashmap_put(vertex_map, vertices[i], i);

    idx_t *xadj = malloc((total + 1) * sizeof(idx_t));
    idx_t *adj = malloc(total * (n - 1) * sizeof(idx_t));
    idx_t *part = malloc(total * sizeof(idx_t));
    if (!xadj || !adj || !part) { printf("Memory allocation failed for graph structures\n"); return 1; }

    printf("Partitioning graph...\n");
    partition_graph(total, n, xadj, adj, part, vertices, vertex_map);

    char *root = malloc((n + 1) * sizeof(char));
    for (int i = 0; i < n; ++i) root[i] = '1' + i;
    root[n] = '\0';
    int root_idx = hashmap_get(vertex_map, root, -1);

    int *parents = malloc((n - 1) * total * sizeof(int));
    if (!parents) { printf("Memory allocation failed for parents\n"); return 1; }
    for (int i = 0; i < (n - 1) * total; ++i) parents[i] = -1;

    int num_threads = omp_get_max_threads();
    if (num_threads > 4) num_threads = 4;
    printf("Using %d threads\n", num_threads);
    
    int **partition_vertices = malloc(num_threads * sizeof(int*));
    int *partition_counts = calloc(num_threads, sizeof(int));
    if (!partition_vertices || !partition_counts) { printf("Memory allocation failed for partitions\n"); return 1; }

    for (int i = 0; i < total; ++i) partition_counts[part[i] % num_threads]++;
    for (int p = 0; p < num_threads; ++p) {
        partition_vertices[p] = malloc(partition_counts[p] * sizeof(int));
        if (!partition_vertices[p]) { printf("Memory allocation failed for partition_vertices %d\n", p); return 1; }
        partition_counts[p] = 0;
    }
    for (int i = 0; i < total; ++i) {
        int p = part[i] % num_threads;
        partition_vertices[p][partition_counts[p]++] = i;
    }

    printf("Vertex distribution across threads:\n");
    for (int p = 0; p < num_threads; ++p) printf("  Thread %d: %d vertices\n", p, partition_counts[p]);

    printf("Computing parent relationships...\n");
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        int error_flag = 0;
        char *local_parent = malloc((n + 1) * sizeof(char));
        if (!local_parent) {
            printf("Thread %d: Memory allocation failed for local_parent\n", tid);
            error_flag = 1;
        }

        for (int t = 0; t < n - 1 && !error_flag; ++t) {
            for (int p = tid; p < num_threads && !error_flag; p += omp_get_num_threads()) {
                for (int j = 0; j < partition_counts[p] && !error_flag; ++j) {
                    int i = partition_vertices[p][j];
                    if (i == root_idx) continue;

                    Parent1(vertices[i], t + 1, n, local_parent, root);
                    int parent_idx = (strcmp(local_parent, "ROOT") == 0) ? root_idx : hashmap_get(vertex_map, local_parent, -1);
                    if (parent_idx == -1) continue;

                    parents[t * total + i] = parent_idx;

                    /* if ((j % 10000) == 0) {
                        #pragma omp critical
                        printf("Thread %d: Processed %d/%d vertices for tree T_%d\n", tid, j, partition_counts[p], t + 1);
                    }*/
                }
            }
            #pragma omp barrier
        }
        free(local_parent);
        if (error_flag) {
            #pragma omp critical
            printf("Thread %d: Error occurred, terminating computation.\n", tid);
        }
    }

    FILE *output_file = fopen("spanning_trees_output.txt", "w");
    if (!output_file) { printf("Error opening output file\n"); return 1; }

    int output_limit = (total < MAX_VERTICES_TO_OUTPUT) ? total : MAX_VERTICES_TO_OUTPUT;

    for (int t = 0; t < n - 1; ++t) {
        char buffer[1024];
        snprintf(buffer, sizeof(buffer), "\nTree T_%d (color: %s):\n", t + 1, colors[t]);
        fwrite(buffer, 1, strlen(buffer), output_file);
        snprintf(buffer, sizeof(buffer), "Parents:\n");
        fwrite(buffer, 1, strlen(buffer), output_file);

        int vertex_count = 0;
        for (int i = 0; i < total; ++i) {
            if (i != root_idx && parents[t * total + i] != -1) vertex_count++;
        }
        snprintf(buffer, sizeof(buffer), "Total vertices with parents: %d\n", vertex_count);
        fwrite(buffer, 1, strlen(buffer), output_file);

        int output_count = 0;
        #pragma omp parallel num_threads(num_threads)
        {
            int tid = omp_get_thread_num();
            char *local_buffer = malloc(OUTPUT_BUFFER_SIZE);
            if (!local_buffer) {
                #pragma omp critical
                printf("Thread %d: Memory allocation failed for local_buffer\n", tid);
            } else {
                local_buffer[0] = '\0';
                size_t local_pos = 0;

                for (int p = tid; p < num_threads; p += omp_get_num_threads()) {
                    for (int j = 0; j < partition_counts[p]; ++j) {
                        int i = partition_vertices[p][j];
                        if (i == root_idx || parents[t * total + i] == -1) continue;

                        if (output_count < output_limit) {
                            char line[256];
                            snprintf(line, sizeof(line), "  Vertex %s -> Parent %s\n", vertices[i], 
                                     parents[t * total + i] == root_idx ? "ROOT" : vertices[parents[t * total + i]]);
                            size_t line_len = strlen(line);
                            if (local_pos + line_len + 1 < OUTPUT_BUFFER_SIZE) {
                                strcpy(local_buffer + local_pos, line);
                                local_pos += line_len;
                                #pragma omp atomic
                                output_count++;
                            } else {
                                #pragma omp critical
                                {
                                    fwrite(local_buffer, 1, local_pos, output_file);
                                    fflush(output_file);
                                }
                                local_pos = 0;
                                strcpy(local_buffer, line);
                                local_pos += line_len;
                            }
                        }
                    }
                }

                if (local_pos > 0) {
                    #pragma omp critical
                    {
                        fwrite(local_buffer, 1, local_pos, output_file);
                        fflush(output_file);
                    }
                }
            }
            free(local_buffer);
        }

        if (output_count >= output_limit) {
            snprintf(buffer, sizeof(buffer), "[Output truncated at %d vertices due to large N]\n", output_limit);
            fwrite(buffer, 1, strlen(buffer), output_file);
        }
    }
    
    fclose(output_file);
    printf("Output written to spanning_trees_output.txt\n");

    free_hashmap(vertex_map);
    for (int p = 0; p < num_threads; ++p) free(partition_vertices[p]);
    free(partition_vertices);
    free(partition_counts);
    for (int i = 0; i < total; ++i) free(vertices[i]);
    free(vertices);
    free(parents);
    free(xadj);
    free(adj);
    free(part);
    free(root);
    return 0;
}
