#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SIZE 5000


// Define clock default starting value
#define CLOCK_MONOTONIC 1

void test_malloc_memset() {
    struct timespec start, stop;
    int nRows = SIZE;

    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Allocate memory with malloc and initialize it with memset
    float *matrix = (float *)malloc(SIZE * sizeof(float));
    if (matrix == NULL) {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }
    memset(matrix, 0, SIZE * sizeof(float)); // Zero initialize the array

    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);

    // Calculate the elapsed time in seconds
    double time_taken = (stop.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (stop.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("malloc + memset time: %e seconds\n", time_taken);

    free(matrix); // Free the allocated memory
}

void test_calloc() {
    struct timespec start, stop;
    int nRows = SIZE;

    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Allocate memory with calloc (automatically zero-initialized)
    float *matrix = (float *)calloc(SIZE, sizeof(float));
    if (matrix == NULL) {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);

    // Calculate the elapsed time in seconds
    double time_taken = (stop.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (stop.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("calloc time: %e seconds\n", time_taken);

    free(matrix); // Free the allocated memory
}

int main() {
    printf("Testing matrix initialization (5000 elements):\n");

    // Test malloc + memset
    test_malloc_memset();

    // Test calloc
    test_calloc();

    return 0;
}
