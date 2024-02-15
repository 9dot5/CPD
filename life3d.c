#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h> // For 64-bit integers

#define N_SPECIES 9

unsigned int seed;

void init_r4uni(int input_seed) {
    seed = input_seed + 987654321;
}

float r4_uni() {
    int seed_in = seed;

    seed ^= (seed << 13);
    seed ^= (seed >> 17);
    seed ^= (seed << 5);

    return 0.5 + 0.2328306e-09 * (seed_in + (int) seed);
}

// Function to generate the initial grid
char ***gen_initial_grid(long long N, float density, int input_seed) {
    int x, y, z;

    char ***grid = malloc(N * sizeof(char **));
    if (grid == NULL) {
        fprintf(stderr, "Failed to allocate grid\n");
        exit(EXIT_FAILURE);
    }

    for (x = 0; x < N; x++) {
        grid[x] = malloc(N * sizeof(char *));
        if (grid[x] == NULL) {
            fprintf(stderr, "Failed to allocate grid\n");
            exit(EXIT_FAILURE);
        }

        grid[x][0] = calloc(N * N, sizeof(char));
        if (grid[x][0] == NULL) {
            fprintf(stderr, "Failed to allocate grid\n");
            exit(EXIT_FAILURE);
        }

        for (y = 1; y < N; y++)
            grid[x][y] = grid[x][0] + y * N;
    }

    init_r4uni(input_seed);
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
            for (z = 0; z < N; z++)
                if (r4_uni() < density)
                    grid[x][y][z] = (int)(r4_uni() * N_SPECIES) + 1;

    return grid;
}

// Function to count neighbors
void count_neighbors(char ***grid, long long N, int x, int y, int z, int64_t *neighbor_counts) {
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
                if (dx == 0 && dy == 0 && dz == 0)
                    continue;
                int nx = (x + dx + N) % N;
                int ny = (y + dy + N) % N;
                int nz = (z + dz + N) % N;
                if (grid[nx][ny][nz] == grid[x][y][z]) {
                    neighbor_counts[grid[nx][ny][nz] - 1]++;
                }
            }
        }
    }
}

// Function to update the grid and track max
void update_grid_and_track_max(char ***grid, long long N, int64_t *species_max_counts, int64_t *species_max_gens, int generation) {
    char ***temp_grid = malloc(N * sizeof(char **));
    if (temp_grid == NULL) {
        fprintf(stderr, "Failed to allocate temp_grid\n");
        exit(EXIT_FAILURE);
    }

    for (int x = 0; x < N; x++) {
        temp_grid[x] = malloc(N * sizeof(char *));
        if (temp_grid[x] == NULL) {
            fprintf(stderr, "Failed to allocate temp_grid\n");
            exit(EXIT_FAILURE);
        }

        temp_grid[x][0] = calloc(N * N, sizeof(char));
        if (temp_grid[x][0] == NULL) {
            fprintf(stderr, "Failed to allocate temp_grid\n");
            exit(EXIT_FAILURE);
        }

        for (int y = 1; y < N; y++) {
            temp_grid[x][y] = temp_grid[x][0] + y * N;
        }
    }

    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                int64_t neighbor_counts[N_SPECIES] = {0};
                count_neighbors(grid, N, x, y, z, neighbor_counts);

                // Apply transition rules for live cell
                if (grid[x][y][z] > 0) {
                    int alive_neighbors = 0;
                    for (int i = 0; i < N_SPECIES; i++) {
                        if (neighbor_counts[i] > 0) {
                            alive_neighbors += neighbor_counts[i];
                            if (alive_neighbors > 13) {
                                break; // Cell dies due to overcrowding
                            }
                        }
                    }
                    if (alive_neighbors >= 4 && alive_neighbors <= 13) {
                        temp_grid[x][y][z] = grid[x][y][z]; // Remains alive
                    } else {
                        temp_grid[x][y][z] = 0; // Dies
                    }
                } else {
                    // Apply transition rules for dead cell
                    int total_neighbors = 0;
                    for (int i = 0; i < N_SPECIES; i++) {
                        total_neighbors += neighbor_counts[i];
                    }
                    if (total_neighbors >= 7 && total_neighbors <= 10) {
                        int max_count = 0, max_species = 0;
                        for (int i = 0; i < N_SPECIES; i++) {
                            if (neighbor_counts[i] > max_count) {
                                max_count = neighbor_counts[i];
                                max_species = i + 1;
                            }
                        }
                        temp_grid[x][y][z] = max_species; // Becomes the majority species
                    } else {
                        temp_grid[x][y][z] = 0; // Remains dead
                    }
                }

                // Update maximum population and generation for each species
                int species = temp_grid[x][y][z];
                if (species > 0) {
                    if (neighbor_counts[species - 1] > species_max_counts[species - 1]) {
                        species_max_counts[species - 1] = neighbor_counts[species - 1];
                        species_max_gens[species - 1] = generation;
                    } else if (neighbor_counts[species - 1] == species_max_counts[species - 1]) {
                        if (generation < species_max_gens[species - 1]) {
                            species_max_gens[species - 1] = generation;
                        }
                    }
                }
            }
        }
    }

    // Free temporary grid memory
    for (int x = 0; x < N; x++) {
        free(temp_grid[x][0]);
        free(temp_grid[x]);
    }
    free(temp_grid);
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <grid_size> <density> <generations> <seed>\n", argv[0]);
        return EXIT_FAILURE;
    }

    long long N;
    float density;
    int generations;
    int input_seed;

    // Parse and validate input arguments
    if (sscanf(argv[1], "%lld", &N) != 1 || N <= 0) {
        fprintf(stderr, "Invalid grid size\n");
        return EXIT_FAILURE;
    }

    if (sscanf(argv[2], "%f", &density) != 1 || density < 0.0 || density > 1.0) {
        fprintf(stderr, "Invalid density\n");
        return EXIT_FAILURE;
    }

    if (sscanf(argv[3], "%d", &generations) != 1 || generations <= 0) {
        fprintf(stderr, "Invalid number of generations\n");
        return EXIT_FAILURE;
    }

    if (sscanf(argv[4], "%d", &input_seed) != 1) {
        fprintf(stderr, "Invalid seed\n");
        return EXIT_FAILURE;
    }

    // Allocate memory for species maximum population and generation arrays
    int64_t *species_max_counts = calloc(N_SPECIES, sizeof(int64_t));
    if (species_max_counts == NULL) {
        fprintf(stderr, "Failed to allocate species_max_counts\n");
        return EXIT_FAILURE;
    }

    int64_t *species_max_gens = calloc(N_SPECIES, sizeof(int64_t));
    if (species_max_gens == NULL) {
        free(species_max_counts);
        fprintf(stderr, "Failed to allocate species_max_gens\n");
        return EXIT_FAILURE;
    }

    // Generate initial grid
    char ***grid = gen_initial_grid(N, density, input_seed);

    // Run simulation for specified generations
    for (int gen = 0; gen < generations; gen++) {
        // Update grid and track maximum population/generation for each species
        update_grid_and_track_max(grid, N, species_max_counts, species_max_gens, gen);
    }

    // Print final output for each species with maximum population
    for (int i = 0; i < N_SPECIES; i++) {
        if (species_max_counts[i] > 0) {
            printf("%d %lld %lld\n", i + 1, species_max_counts[i], species_max_gens[i]);
        }
    }

    // Free allocated memory
    for (int x = 0; x < N; x++) {
        free(grid[x][0]);
        free(grid[x]);
    }
    free(grid);
    free(species_max_counts);
    free(species_max_gens);

    return EXIT_SUCCESS;
}
