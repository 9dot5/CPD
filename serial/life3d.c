/*
    File: life3d.c
    Description: An implementation of Game of Life 3D (Serial version)
    Authors: Group 22
    Course: Parallel and Distributed Computing - 2023/2024
*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define N_SPECIES 9

/*
    Generate the initial grid
*/

unsigned int seed;

void init_r4uni(int input_seed)
{
    seed = input_seed + 987654321;
}

float r4_uni()
{
    int seed_in = seed;

    seed ^= (seed << 13);
    seed ^= (seed >> 17);
    seed ^= (seed << 5);

    return 0.5 + 0.2328306e-09 * (seed_in + (int) seed);
}

char ***gen_initial_grid(long long N, float density, int input_seed)
{
    int x, y, z;

    char ***grid = (char ***) malloc(N * sizeof(char **));
    if(grid == NULL) {
        printf("Failed to allocate matrix\n");
        exit(1);
    }
    for(x = 0; x < N; x++) {
        grid[x] = (char **) malloc(N * sizeof(char *));
        if(grid[x] == NULL) {
            printf("Failed to allocate matrix\n");
            exit(1);
        }
        grid[x][0] = (char *) calloc(N * N, sizeof(char));
        if(grid[x][0] == NULL) {
            printf("Failed to allocate matrix\n");
            exit(1);
        }
        for (y = 1; y < N; y++)
            grid[x][y] = grid[x][0] + y * N;
    }

    init_r4uni(input_seed);
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
            for (z = 0; z < N; z++)
                if(r4_uni() < density)
                    grid[x][y][z] = (int)(r4_uni() * N_SPECIES) + 1;

    return grid;
}

/*
    Run the simulation
*/

int debug = 0;

void count_neighbors(char ***grid, int cells, int x, int y, int z, int *species, int *alive_count, int *most_common_value) {
    int count = 0;
    
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int k = -1; k <= 1; k++) {

                if (i == 0 && j == 0 && k == 0) continue;

                int new_x = (x + i + cells) % cells;
                int new_y = (y + j + cells) % cells;
                int new_z = (z + k + cells) % cells;

                int neighbor_value = grid[new_x][new_y][new_z];

                if (neighbor_value == 0) continue;
                
                count++;

                species[neighbor_value - 1]++;
            }
        }
    }

    int max = 0;

    for (int i = 0; i < N_SPECIES; i++) {
        if (species[i] > max) {
            max = species[i];
            *most_common_value = i + 1;
        }
    }

    *alive_count = count;
}

void next_gen(char ***grid, int cells) {
    // Do a copy of the grid
    char ***grid_temp = (char ***) malloc(cells * sizeof(char **));
    if (grid_temp == NULL) {
        fprintf(stderr, "Failed to allocate memory for grid temp\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < cells; i++) {
        grid_temp[i] = (char **) malloc(cells * sizeof(char *));
        if (grid_temp[i] == NULL) {
            fprintf(stderr, "Failed to allocate memory for grid temp\n");
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < cells; j++) {
            grid_temp[i][j] = (char *) malloc(cells * sizeof(char));
            if (grid_temp[i][j] == NULL) {
                fprintf(stderr, "Failed to allocate memory for grid temp\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    for (int i = 0; i < cells; i++) {
        for (int j = 0; j < cells; j++) {
            for (int k = 0; k < cells; k++) {
                grid_temp[i][j][k] = grid[i][j][k];
            }
        }
    }

    // Update the grid
    for (int x = 0; x < cells; x++) {
        for (int y = 0; y < cells; y++) {
            for (int z = 0; z < cells; z++) {

                int neighbors_species[N_SPECIES] = {0};
                int neighbors_count = 0, neighbors_most_common_value = 0;
                
                count_neighbors(grid_temp, cells, x, y, z, neighbors_species, &neighbors_count, &neighbors_most_common_value);

                // Dead cell
                if (grid_temp[x][y][z] == 0) {
                    if (neighbors_count >= 7 && neighbors_count <= 10) {
                        grid[x][y][z] = neighbors_most_common_value;
                    }
                }

                // Alive cell
                else {
                    if (neighbors_count <= 4 || neighbors_count > 13) {
                        grid[x][y][z] = 0;
                    }
                }
            }
        }
    }

    // Free the copy
    for (int i = 0; i < cells; i++) {
        for (int j = 0; j < cells; j++) {
            free(grid_temp[i][j]);
        }
        free(grid_temp[i]);
    }
    free(grid_temp);
}

void print_result_verbose(char ***grid, int cells) {
    for (int x = 0; x < cells; x++) {
        printf("Layer %d:\n", x);

        for (int y = 0; y < cells; y++) {
            for (int z = 0; z < cells; z++) {
                if (grid[x][y][z] == 0) {
                    printf("  ");
                } else {
                    printf("%d ", grid[x][y][z]);
                }
            }
            printf("\n");
        }
        printf("\n");
    }
}

void simulation(char ***grid, int cells, int generations, long long *species_max_count, long long *species_max_count_generation) {
    // Generation 0 being the initial state
    for (int i = 0; i < generations + 1; i++) {

        if (i > 0) next_gen(grid, cells);

        // Update species
        long long species_temp[N_SPECIES] = {0};
        
        for (int x = 0; x < cells; x++) {
            for (int y = 0; y < cells; y++) {
                for (int z = 0; z < cells; z++) {
                    if (grid[x][y][z] == 0)
                        continue;
                    
                    species_temp[grid[x][y][z] - 1]++;
                }
            }
        }

        for (int j = 0; j < N_SPECIES; j++) {
            if (species_temp[j] > species_max_count[j]) {
                species_max_count[j] = species_temp[j];
                species_max_count_generation[j] = i;
            }
        }

        if (debug) {
            printf("Generation %d ------------------------------\n", i);
            print_result_verbose(grid, cells);
        }
    }
}

void print_result(char ***grid, int cells, long long *species_max_count, long long *species_max_count_generation) {
    for (int i = 0; i < N_SPECIES; i++) {
        printf("%d %lld %lld\n", i + 1, species_max_count[i], species_max_count_generation[i]);
    }
}

int main(int argc, char *argv[]) {

    // Check for correct usage
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <1: Positive Int> <2: Positive Int> <3: Float between 0 and 1> <4: Int>\n", argv[0]);
        fprintf(stderr, "1: Number of generations\n");
        fprintf(stderr, "2: Number of cells per side of the cube\n");
        fprintf(stderr, "3: Density of the initial population\n");
        fprintf(stderr, "4: Seed for the random number generator\n");
        return EXIT_FAILURE;
    }

    // Parse the input
    int generations = atoi(argv[1]);
    int cells = atoi(argv[2]);
    float density = atof(argv[3]);
    int seed = atoi(argv[4]);

    // Print the input
    if (debug) printf("Arguments: %d %d %f %d\n", generations, cells, density, seed);

    // Allocate memory for the species
    long long *species_max_count = calloc(N_SPECIES, sizeof(long long));
    long long *species_max_count_generation = calloc(N_SPECIES, sizeof(long long));

    if (species_max_count == NULL || species_max_count_generation == NULL) {
        fprintf(stderr, "Failed to allocate memory for species\n");
        return EXIT_FAILURE;
    }

    // Generate the initial grid
    char ***grid = gen_initial_grid(cells, density, seed);

    // Run the simulation and measure its execution time
    double exec_time = -omp_get_wtime();

    simulation(grid, cells, generations, species_max_count, species_max_count_generation);

    exec_time += omp_get_wtime();

    // Print result
    print_result(grid, cells, species_max_count, species_max_count_generation);

    // Print execution time
    fprintf(stderr, "%.1fs\n", exec_time);

    // Free allocated memory
    for (int x = 0; x < cells; x++) {
        free(grid[x][0]);
        free(grid[x]);
    }
    free(grid);

    free(species_max_count);
    free(species_max_count_generation);

    return EXIT_SUCCESS;
}