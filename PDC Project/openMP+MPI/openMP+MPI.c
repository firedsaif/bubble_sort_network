#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>

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
void serialize_hashmap(HashMap* map, char** buffer, int* size);
HashMap* deserialize_hashmap(char* buffer, int size);

int factorial(int n);
void swap_adjacent(char *perm, int i, int n);
void get_permutation(int idx, int n, char *out);
int find_position(char *perm, int n);
int index_of(char *perm, int n, char symbol);
void Parent1(char *v, int t, int n, char *parent, char *root, int *parent_idx, int rank, HashMap* vertex_map);
void partition_graph(int total, int n, idx_t *xadj, idx_t *adj, idx_t *part, char **vertices, HashMap* vertex_map, int nparts);

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
        out[i] = (nums[selected] == 10) ? '0' : ('0' + nums[selected]);
        for (int j = selected; j < n - 1 - i; ++j) nums[j] = nums[j + 1];
        idx = idx % fact;
    }
    out[n] = '\0';
    free(nums);
}

int find_position(char *perm, int n) {
    for (int i = 0; i < n - 1; ++i) {
        if (perm[i] > perm[i + 1]) return i;
    }
    return -1;
}

int index_of(char *perm, int n, char symbol) {
    for (int i = 0; i < n; ++i)
        if (perm[i] == symbol) return i;
    return -1;
}

void Parent1(char *v, int t, int n, char *parent, char *root, int *parent_idx, int rank, HashMap* vertex_map) {
    if (!v || !root) {
        printf("Rank %d: Null pointer in Parent1 (v=%p, root=%p)\n", rank, (void*)v, (void*)root);
        *parent_idx = -1;
        return;
    }

    if (strcmp(v, root) == 0) {
        strcpy(parent, "ROOT");
        *parent_idx = 0;
        return;
    }

    strcpy(parent, v);
    int pos = t - 1;
    if (pos >= 0 && pos < n - 1 && parent[pos] > parent[pos + 1]) {
        swap_adjacent(parent, pos, n);
    }

    *parent_idx = hashmap_get(vertex_map, parent, -1);
}

void serialize_hashmap(HashMap* map, char** buffer, int* size) {
    int total_size = 0;
    for (int i = 0; i < map->capacity; i++) {
        if (map->keys[i]) {
            total_size += strlen(map->keys[i]) + 1 + sizeof(int);
        }
    }
    *size = total_size;
    *buffer = malloc(*size);
    if (!*buffer) {
        printf("Error: Failed to allocate buffer for hashmap serialization\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int offset = 0;
    for (int i = 0; i < map->capacity; i++) {
        if (map->keys[i]) {
            strcpy(*buffer + offset, map->keys[i]);
            offset += strlen(map->keys[i]) + 1;
            memcpy(*buffer + offset, &map->values[i], sizeof(int));
            offset += sizeof(int);
        }
    }
}

HashMap* deserialize_hashmap(char* buffer, int size) {
    HashMap* map = create_hashmap(size / (sizeof(char*) + sizeof(int)));
    int offset = 0;
    while (offset < size) {
        char* key = buffer + offset;
        offset += strlen(key) + 1;
        int value;
        memcpy(&value, buffer + offset, sizeof(int));
        offset += sizeof(int);
        hashmap_put(map, key, value);
    }
    return map;
}

const char *colors[] = {"red", "blue", "green", "orange", "purple", "brown", "cyan", "magenta", "yellow", "black"};

void partition_graph(int total, int n, idx_t *xadj, idx_t *adj, idx_t *part, char **vertices, HashMap* vertex_map, int nparts) {
    idx_t num_vertices = total;
    idx_t ncon = 1;
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

        adj_lists[i] = malloc(local_edge_idx * sizeof(idx_t));
        for (idx_t k = 0; k < local_edge_idx; k++) {
            adj_lists[i][k] = local_adj[k];
        }
        adj_sizes[i] = local_edge_idx;
        free(local_adj);
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

    printf("Calling METIS_PartGraphKway with %lld vertices, %lld parts\n", (long long)num_vertices, (long long)nparts);
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

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 10;
    int total = factorial(n);
    if (total > 50000000) {
        if (rank == 0) printf("N=%d results in %d vertices, which exceeds practical limits. Aborting.\n", n, total);
        MPI_Finalize();
        return 1;
    }

    double start_time, end_time;

    // Broadcast total and n to all ranks
    start_time = MPI_Wtime();
    MPI_Bcast(&total, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    if (rank == 0) printf("Time for broadcasting total and n: %.3f seconds\n", end_time - start_time);

    // Allocate vertices on all ranks
    char **vertices = malloc(total * sizeof(char*));
    if (!vertices) { printf("Rank %d: Memory allocation failed for vertices\n", rank); MPI_Abort(MPI_COMM_WORLD, 1); }
    for (int i = 0; i < total; ++i) {
        vertices[i] = malloc((n + 1) * sizeof(char));
        if (!vertices[i]) { printf("Rank %d: Memory allocation failed for vertex %d\n", rank, i); MPI_Abort(MPI_COMM_WORLD, 1); }
        vertices[i][0] = '\0';
    }

    // Allocate vertex_buffer on all ranks
    char *vertex_buffer = malloc(total * (n + 1) * sizeof(char));
    if (!vertex_buffer) { printf("Rank %d: Memory allocation failed for vertex_buffer\n", rank); MPI_Abort(MPI_COMM_WORLD, 1); }

    if (rank == 0) {
        for (int i = 0; i < total; ++i) {
            get_permutation(i, n, vertices[i]);
            memcpy(vertex_buffer + i * (n + 1), vertices[i], (n + 1) * sizeof(char));
        }
    }

    // Broadcast the flattened vertex buffer
    start_time = MPI_Wtime();
    MPI_Bcast(vertex_buffer, total * (n + 1), MPI_CHAR, 0, MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    if (rank == 0) printf("Time for broadcasting vertex_buffer: %.3f seconds\n", end_time - start_time);

    // Reconstruct vertices on all ranks
    for (int i = 0; i < total; ++i) {
        memcpy(vertices[i], vertex_buffer + i * (n + 1), (n + 1) * sizeof(char));
    }

    // Debug: Print first few vertices to verify broadcast
    if (rank == 0) printf("Rank %d: First 3 vertices after broadcast:\n", rank);
    for (int i = 0; i < 3 && i < total; ++i) {
        printf("Rank %d: vertices[%d] = %s\n", rank, i, vertices[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    HashMap* vertex_map = NULL;
    idx_t *xadj = NULL, *adj = NULL, *part = NULL;
    int *parents = NULL;
    int **partition_vertices = NULL;
    int *partition_counts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));
    if (!partition_counts || !displs) {
        printf("Rank %d: Memory allocation failed for partition_counts or displs\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int i = 0; i < size; ++i) {
        partition_counts[i] = 0;
        displs[i] = 0;
    }

    if (rank == 0) {
        printf("Constructing spanning trees for bubble-sort network with N=%d (%d vertices)\n", n, total);

        vertex_map = create_hashmap(total);
        for (int i = 0; i < total; i++) {
            hashmap_put(vertex_map, vertices[i], i);
        }

        xadj = malloc((total + 1) * sizeof(idx_t));
        adj = malloc(total * (n - 1) * sizeof(idx_t));
        part = malloc(total * sizeof(idx_t));
        if (!xadj || !adj || !part) { printf("Memory allocation failed for graph structures\n"); MPI_Abort(MPI_COMM_WORLD, 1); }

        printf("Partitioning graph...\n");
        partition_graph(total, n, xadj, adj, part, vertices, vertex_map, size);

        parents = malloc((n - 1) * total * sizeof(int));
        if (!parents) { printf("Memory allocation failed for parents\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
        for (int i = 0; i < (n - 1) * total; ++i) parents[i] = -1;

        partition_vertices = malloc(size * sizeof(int*));
        if (!partition_vertices) { printf("Memory allocation failed for partitions\n"); MPI_Abort(MPI_COMM_WORLD, 1); }

        for (int i = 0; i < total; ++i) partition_counts[part[i] % size]++;
        displs[0] = 0;
        for (int p = 0; p < size; ++p) {
            partition_vertices[p] = malloc(partition_counts[p] * sizeof(int));
            if (!partition_vertices[p]) { printf("Memory allocation failed for partition_vertices %d\n", p); MPI_Abort(MPI_COMM_WORLD, 1); }
            partition_counts[p] = 0;
            if (p > 0) displs[p] = displs[p-1] + partition_counts[p-1];
        }
        for (int i = 0; i < total; ++i) {
            int p = part[i] % size;
            partition_vertices[p][partition_counts[p]++] = i;
        }

        // Flatten partition_vertices into a contiguous array
        int *flat_partition_vertices = malloc(total * sizeof(int));
        if (!flat_partition_vertices) { printf("Memory allocation failed for flat_partition_vertices\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
        int offset = 0;
        int total_vertices = 0;
        for (int p = 0; p < size; p++) {
            memcpy(flat_partition_vertices + displs[p], partition_vertices[p], partition_counts[p] * sizeof(int));
            offset += partition_counts[p];
            total_vertices += partition_counts[p];
            printf("Debug: Process %d partition size: %d, offset: %d\n", p, partition_counts[p], offset);
        }
        printf("Debug: Total vertices in flat_partition_vertices: %d, expected: %d\n", total_vertices, total);

        printf("Vertex distribution across processes:\n");
        for (int p = 0; p < size; ++p) printf("  Process %d: %d vertices\n", p, partition_counts[p]);
    }

    // Broadcast partition counts and displs
    start_time = MPI_Wtime();
    MPI_Bcast(partition_counts, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, size, MPI_INT, 0, MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    if (rank == 0) printf("Time for broadcasting partition_counts and displs: %.3f seconds\n", end_time - start_time);

    int local_count = partition_counts[rank];
    int *local_partition = malloc(local_count * sizeof(int));
    if (!local_partition) { printf("Rank %d: Memory allocation failed for local_partition\n", rank); MPI_Abort(MPI_COMM_WORLD, 1); }

    // Use MPI_Scatterv to distribute flat_partition_vertices
    int *flat_partition_vertices = NULL;
    if (rank == 0) {
        flat_partition_vertices = malloc(total * sizeof(int));
        if (!flat_partition_vertices) { printf("Rank 0: Memory allocation failed for flat_partition_vertices\n"); MPI_Abort(MPI_COMM_WORLD, 1); }
        int offset = 0;
        for (int p = 0; p < size; p++) {
            memcpy(flat_partition_vertices + displs[p], partition_vertices[p], partition_counts[p] * sizeof(int));
            offset += partition_counts[p];
        }
        for (int i = 0; i < size; i++) free(partition_vertices[i]);
        free(partition_vertices);
    }

    start_time = MPI_Wtime();
    MPI_Scatterv(flat_partition_vertices, partition_counts, displs, MPI_INT,
                 local_partition, local_count, MPI_INT, 0, MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    if (rank == 0) printf("Time for MPI_Scatterv flat_partition_vertices: %.3f seconds\n", end_time - start_time);

    if (rank == 0) free(flat_partition_vertices);

    // Compute parent relationships locally with OpenMP
    int *local_parents = malloc((n - 1) * local_count * sizeof(int));
    if (!local_parents) { printf("Rank %d: Memory allocation failed for local_parents\n", rank); MPI_Abort(MPI_COMM_WORLD, 1); }
    for (int i = 0; i < (n - 1) * local_count; ++i) local_parents[i] = -1;

    char **local_vertices = malloc(local_count * sizeof(char*));
    if (!local_vertices) { printf("Rank %d: Memory allocation failed for local_vertices\n", rank); MPI_Abort(MPI_COMM_WORLD, 1); }
    for (int j = 0; j < local_count; ++j) {
        local_vertices[j] = malloc((n + 1) * sizeof(char));
        if (!local_vertices[j]) { printf("Rank %d: Memory allocation failed for local_vertices[%d]\n", rank, j); MPI_Abort(MPI_COMM_WORLD, 1); }
        strcpy(local_vertices[j], vertices[local_partition[j]]);
    }

    // Debug: Print first few local_vertices
    if (local_count > 0) {
        printf("Rank %d: First 3 local_vertices (total %d):\n", rank, local_count);
        for (int j = 0; j < 3 && j < local_count; ++j) {
            printf("Rank %d: local_vertices[%d] = %s (index %d)\n", rank, j, local_vertices[j], local_partition[j]);
        }
    }

    // Broadcast vertex_map
    HashMap* local_vertex_map = NULL;
    if (rank == 0) {
        char* map_buffer;
        int map_size;
        serialize_hashmap(vertex_map, &map_buffer, &map_size);
        start_time = MPI_Wtime();
        MPI_Bcast(&map_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(map_buffer, map_size, MPI_CHAR, 0, MPI_COMM_WORLD);
        end_time = MPI_Wtime();
        if (rank == 0) printf("Time for broadcasting vertex_map: %.3f seconds\n", end_time - start_time);
        free(map_buffer);
    } else {
        int map_size;
        start_time = MPI_Wtime();
        MPI_Bcast(&map_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        char* map_buffer = malloc(map_size);
        if (!map_buffer) { printf("Rank %d: Memory allocation failed for map_buffer\n", rank); MPI_Abort(MPI_COMM_WORLD, 1); }
        MPI_Bcast(map_buffer, map_size, MPI_CHAR, 0, MPI_COMM_WORLD);
        end_time = MPI_Wtime();
        if (rank == 0) printf("Time for broadcasting vertex_map: %.3f seconds\n", end_time - start_time);
        local_vertex_map = deserialize_hashmap(map_buffer, map_size);
        free(map_buffer);
    }
    if (rank == 0) local_vertex_map = vertex_map;

    int num_threads = omp_get_max_threads();
    if (rank == 0) printf("Using %d threads per process\n", num_threads);

    printf("Rank %d: Computing parent relationships for %d vertices...\n", rank, local_count);
    #pragma omp parallel for num_threads(num_threads) collapse(2) schedule(dynamic)
    for (int t = 0; t < n - 1; ++t) {
        for (int j = 0; j < local_count; ++j) {
            char local_parent[n + 1];
            int parent_idx = -1;
            Parent1(local_vertices[j], t + 1, n, local_parent, vertices[0], &parent_idx, rank, local_vertex_map);
            if (parent_idx != -1) {
                #pragma omp critical
                local_parents[t * local_count + j] = parent_idx;
            }
        }
    }

    // Free local_vertices
    for (int j = 0; j < local_count; ++j) free(local_vertices[j]);
    free(local_vertices);

    // Use MPI_Gatherv to collect local_parents
    int *parent_counts = malloc(size * sizeof(int));
    int *parent_displs = malloc(size * sizeof(int));
    if (!parent_counts || !parent_displs) {
        printf("Rank %d: Memory allocation failed for parent_counts or parent_displs\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int i = 0; i < size; i++) {
        parent_counts[i] = (n - 1) * partition_counts[i];
        parent_displs[i] = (i == 0) ? 0 : parent_displs[i-1] + parent_counts[i-1];
    }

    start_time = MPI_Wtime();
    if (rank == 0) {
        MPI_Gatherv(local_parents, (n - 1) * local_count, MPI_INT,
                    parents, parent_counts, parent_displs, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        MPI_Gatherv(local_parents, (n - 1) * local_count, MPI_INT,
                    NULL, NULL, NULL, MPI_DATATYPE_NULL, 0, MPI_COMM_WORLD);
    }
    end_time = MPI_Wtime();
    if (rank == 0) printf("Time for MPI_Gatherv local_parents: %.3f seconds\n", end_time - start_time);

    if (rank == 0) {
        FILE *output_file = fopen("spanning_trees_output.txt", "w");
        if (!output_file) { printf("Error opening output file\n"); MPI_Abort(MPI_COMM_WORLD, 1); }

        int output_limit = (total < MAX_VERTICES_TO_OUTPUT) ? total : MAX_VERTICES_TO_OUTPUT;

        for (int t = 0; t < n - 1; ++t) {
            char buffer[1024];
            snprintf(buffer, sizeof(buffer), "\nTree T_%d (color: %s):\n", t + 1, colors[t]);
            fwrite(buffer, 1, strlen(buffer), output_file);
            snprintf(buffer, sizeof(buffer), "Parents:\n");
            fwrite(buffer, 1, strlen(buffer), output_file);

            int vertex_count = 0;
            for (int i = 0; i < total; ++i) {
                if (i != 0 && parents[t * total + i] != -1) vertex_count++;
            }
            snprintf(buffer, sizeof(buffer), "Total vertices with parents: %d\n", vertex_count);
            fwrite(buffer, 1, strlen(buffer), output_file);

            int output_count = 0;
            for (int j = 0; j < total; ++j) {
                int i = j;
                if (i == 0 || parents[t * total + i] == -1) continue;

                if (output_count < output_limit) {
                    char line[256];
                    snprintf(line, sizeof(line), "  Vertex %s -> Parent %s\n", vertices[i],
                             parents[t * total + i] == 0 ? "ROOT" : vertices[parents[t * total + i]]);
                    fwrite(line, 1, strlen(line), output_file);
                    output_count++;
                }
            }

            if (output_count >= output_limit) {
                snprintf(buffer, sizeof(buffer), "[Output truncated at %d vertices due to large N]\n", output_limit);
                fwrite(buffer, 1, strlen(buffer), output_file);
            }
        }

        fclose(output_file);
        printf("Output written to spanning_trees_output.txt\n");

        free_hashmap(vertex_map);
        free(parents);
        free(xadj);
        free(adj);
        free(part);
    }

    for (int i = 0; i < total; ++i) free(vertices[i]);
    free(vertices);
    free(local_partition);
    free(local_parents);
    free(vertex_buffer);
    free(partition_counts);
    free(displs);
    free(parent_counts);
    free(parent_displs);
    if (rank != 0) free_hashmap(local_vertex_map);

    MPI_Finalize();
    return 0;
}
