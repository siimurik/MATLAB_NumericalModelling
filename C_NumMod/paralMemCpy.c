#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

void parallel_memcpy(void *dest, const void *src, size_t n) {
    size_t chunk_size = n / omp_get_max_threads();
    
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        size_t start = thread_id * chunk_size;
        size_t end = (thread_id == omp_get_max_threads() - 1) ? n : start + chunk_size;

        memcpy((char *)dest + start, (const char *)src + start, end - start);
    }
}

int main() {
    size_t size = 100000000; // 100 million bytes
    char *src = malloc(size);
    char *dest = malloc(size);

    // Initialize source array
    for (size_t i = 0; i < size; i++) {
        src[i] = (char)(i % 256);
    }

    // Perform parallel memcpy
    parallel_memcpy(dest, src, size);

    // Verify the copy
    if (memcmp(src, dest, size) == 0) {
        printf("Memory copied successfully!\n");
    } else {
        printf("Memory copy failed!\n");
    }

    free(src);
    free(dest);
    return 0;
}