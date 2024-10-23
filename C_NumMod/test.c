/* Two ways to compie this code.
Method 1:
    $ gcc -o test test.c splash.c -llapacke -lblas -lm
    
Method 2:
    $ gcc -c splash.c -o splash.o
    $ ar rcs libsplash.a splash.o
    $ gcc test.c -L. -lsplash -llapacke -lblas -lm -o test
*/
#include <stdio.h> 		// for `printf` to work
#include <math.h> 		// mathematical functions (exp, cos, log, ...)
#include <stdlib.h> 	// For EXIT_FAILURE to work
#include <lapacke.h> 	// Linear Algebra PACKage in Fortran
#include <cblas.h> 		// Basic Linear Algebra Subprograms in C	
#include <string.h> 	// For the memcpy() command
#include <time.h>
#include "splash.h"

// Define matrix dimesions
#define DIM 10

// Define clock default starting value
#define CLOCK_MONOTONIC 1

/* All handeled by the "splash.h" file
typedef struct {
	float *data;
	int rows;
	int cols;
} Matrix;


// Function to free the allocated memory
void freeMatrix(Matrix *matrix) {
    free(matrix->data);
    matrix->data = NULL; // Avoid dangling pointer
}

// Generate a random floating-point number within the given range.
float rand_range(float min, float max) {
    return min + (max-min) * ((float)rand()/RAND_MAX);
}

// Function to generate a new Matrix with random values
Matrix matrand(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (float *)malloc(rows * cols * sizeof(float));

    float min = 0.0;
    float max = 1.0;

    for (int i = 0; i < rows * cols; i++) {
        mat.data[i] = min + (max-min) * ((float)rand()/RAND_MAX);
    }

    return mat;
}

// Matrix multiplication function
Matrix matmul(const Matrix *A, const Matrix *B) {
    if (A->cols != B->rows) {
        fprintf(stderr, "Matrix dimensions are incompatible for multiplication.\n");
        exit(EXIT_FAILURE);
    }

    Matrix C;
    C.rows = A->rows;
    C.cols = B->cols;
    C.data = (float *)malloc(C.rows * C.cols * sizeof(float));

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, C.rows, C.cols, A->cols, 
					 1.0, A->data, A->cols, B->data, B->cols, 0.0, C.data, C.rows);

    return C;
}

void printMatrix(const Matrix *matvec) {
    int max_print_size = 6;

    printf("[\n");

    if (matvec->rows <= max_print_size) {
        for (int i = 0; i < matvec->rows; i++) {
            printf("  [");
            for (int j = 0; j < matvec->cols; j++) {
                printf("%9.4f", matvec->data[i * matvec->cols + j]);
                if (j < matvec->cols - 1) printf(", ");
            }
            printf("]");
            if (i < matvec->rows - 1) printf(",\n");
        }
    } else {
        for (int i = 0; i < max_print_size; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size; j++) {
                printf("%9.4f", matvec->data[i * matvec->cols + j]);
                if (j < matvec->cols - 1) printf(", ");
            }
            printf(" ...");
            printf("]");
            if (i < max_print_size - 1) printf(",\n");
        }
        printf(",\n  ...\n");
        printf("  ...\n");
        printf("  ...");
    }

    printf("\n]\n");
}
*/
//=================================================================
int main(){
    struct timespec start, stop;
    int rows = DIM;
    int cols = DIM;

    Matrix A = matrand(rows, cols);
    Matrix B = matrand(rows, cols);

    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);
    Matrix C = matmul(&A, &B);

    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);
    // Calculate the elapsed time in seconds
    double time_taken = (stop.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (stop.tv_nsec - start.tv_nsec)) * 1e-9;

    printf("Matrix C:\n");
    printMatrix(&C);
    printf("\nMatrix multiplication took %.3e seconds.\n", time_taken);


    freeMatrix(&A);
    freeMatrix(&B);
    freeMatrix(&C);
    return 0;
}