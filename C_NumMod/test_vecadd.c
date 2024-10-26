// Compile and execute with:
//  $ gcc -O3 -fopenmp -std=c99 test_vecadd.c -o test -llapacke -lblas
//  $ gcc test_vecadd.c  -L. -lsplash -llapacke -lblas -o test
//  $ gcc test_vecadd.c splash.c -o test -llapacke -lblas
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>
#include <time.h>
#include <omp.h>
#include "splash.h"
//#include <immintrin.h> // For SIMD intrinsics

// Define clock default starting value
#define CLOCK_MONOTONIC 1
/*
typedef struct {
    float *data;
    int rows;
    int cols;
} Matrix;

typedef struct {
	float *data;
	int    size;
} Vector;

void freeMatrix(Matrix *matrix) {
    free(matrix->data);
    matrix->data = NULL;    // Avoid dangling pointer.
} 

// Function to free the allocated memory in a vector
void freeVector(Vector *vector) {
    free(vector->data);
    vector->data = NULL; // Avoid dangling pointer
}


// Function to create a vector from array values
Vector createVector(int size, float *values) {
    Vector vector;
    vector.size = size;
    vector.data = (float *)malloc(size * sizeof(float)); // Allocate memory

    if (vector.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Populate the vector with the provided values
    memcpy(vector.data, values, vector.size * sizeof(float));

    return vector;
}

// Vector addition function
Vector vecadd(const Vector *A, const Vector *B) {
    if (A->size != B->size) {
        fprintf(stderr, "Vector dimensions are incompatible for addition.\n");
        exit(EXIT_FAILURE);
    }

    Vector C;
    C.size = A->size;
    C.data = (float *)malloc(C.size * sizeof(float));
    memcpy(C.data, B->data, C.size * sizeof(float));

    // Constant times a vector plus a vector
    // ALPHA * A + C;   ALPHA = 1.0;
    cblas_saxpy(C.size, 1.0, A->data, 1, C.data, 1);

    return C;
}

// Vector scalar multiplication function
Vector vecscal(const Vector *A, float scalar) {
    Vector C;
    C.size = A->size;
    C.data = (float *)malloc(C.size * sizeof(float));
    memcpy(C.data, A->data, C.size * sizeof(float));

    // Scales a vector by a constant
    cblas_sscal(C.size, scalar, C.data, 1);

    return C;
}

Vector vecrand(int dim) {
    Vector vec;
    vec.size = dim;
    vec.data = (float *)malloc(dim * sizeof(float));

    if (vec.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    float min = 0.0;
    float max = 1.0;
	for (int i = 0; i < vec.size; i++) {
	    vec.data[i] = min + (max - min) * ((float)rand() / (float)RAND_MAX);
	}

    return vec;
}

Vector vecrand(int dim) {
    Vector vec;
    vec.size = dim;
    vec.data = (float *)malloc(dim * sizeof(float));

    if (vec.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    float min = 0.0;
    float max = 1.0;
    
    // llelize the random number generation
    #pragma omp llel
    {
        // Each thread will have its own random seed
        unsigned int seed = omp_get_thread_num(); // Use thread number as seed

        #pragma omp for
        for (int i = 0; i < vec.size; i++) {
            // Generate a random number using the thread-local seed
            vec.data[i] = min + (max - min) * ((float)rand_r(&seed) / (float)RAND_MAX);
        }
    }
    
    return vec;
}

int length(const Vector vec){
	return vec.size;
}

void printVector(const Vector *vector) {
    printf("[");
    for (int i = 0; i < vector->size; i++) {
        printf("%.4f", vector->data[i]);
        if (i < vector->size - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

void printVectorBetter(const Vector *vector) {
    printf("[");
    
    if (vector->size > 10) {
        // Print the first three elements
        for (int i = 0; i < 3; i++) {
            printf("%.4f", vector->data[i]);
            if (i < 2) {
                printf(", ");
            }
        }
        // Print ellipsis
        printf(", ...");
        // Print the last three elements
        for (int i = vector->size - 3; i < vector->size; i++) {
            printf(", %.4f", vector->data[i]);
        }
    } else {
        // Print all elements if size is 10 or less
        for (int i = 0; i < vector->size; i++) {
            printf("%.4f", vector->data[i]);
            if (i < vector->size - 1) {
                printf(", ");
            }
        }
    }
    
    printf("]\n");
}


void printMatrix(const Matrix *mat) {
    int max_print_size = 3; // We want to print 3 rows and 3 columns from corners

    printf("[\n");

    if (mat->rows <= 10 && mat->cols <= 10) {
        // Print the entire matrix if both dimensions are 10 or less
        for (int i = 0; i < mat->rows; i++) {
            printf("  [");
            for (int j = 0; j < mat->cols; j++) {
                printf("%7.4f", mat->data[i * mat->cols + j]);
                if (j < mat->cols - 1) printf(", ");
            }
            printf("]");
            if (i < mat->rows - 1) printf(",\n");
        }
    } else {
        // Print the top 3 rows
        for (int i = 0; i < max_print_size && i < mat->rows; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size && j < mat->cols; j++) {
                printf("%7.4f", mat->data[i * mat->cols + j]);
                if (j < max_print_size - 1 && j < mat->cols - 1) printf(", ");
            }
            // Print ellipsis for omitted columns
            if (mat->cols > max_print_size) {
                printf(", ... ");
            }
            // Last elements of the top 3 rows
            for (int j = mat->cols - 3; j < mat->cols; j++) {
                if (j >= 0) { // Ensure we don't go out of bounds
                    printf("%7.4f", mat->data[i * mat->cols + j]);
                    if (j < mat->cols - 1) printf(", ");
                }
            }
            printf("]");
            if (i < max_print_size - 1 && i < mat->rows - 1) printf(",\n");
        }

        // Print ellipsis for omitted rows
        if (mat->rows > max_print_size) {
            printf(",\n\t\t\t      ...\n");
        }

        // Print the bottom 3 rows
        for (int i = mat->rows - max_print_size; i < mat->rows; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size && j < mat->cols; j++) {
                printf("%7.4f", mat->data[i * mat->cols + j]);
                if (j < max_print_size - 1 && j < mat->cols - 1) printf(", ");
            }
            // Print ellipsis for omitted columns
            if (mat->cols > max_print_size) {
                printf(", ... ");
            }
            // Last elements of the bottom 3 rows
            for (int j = mat->cols - 3; j < mat->cols; j++) {
                if (j >= 0) { // Ensure we don't go out of bounds
                    printf("%7.4f", mat->data[i * mat->cols + j]);
                    if (j < mat->cols - 1) printf(", ");
                }
            }
            printf("]");
            if (i < mat->rows - 1) printf(",\n");
        }
    }

    printf("\n]\n");
}

// Function to generate a new Matrix with random values
Matrix matrand(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (float *)malloc(rows * cols * sizeof(float));

    if (mat.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    float min = 0.0;
    float max = 1.0;
	for (int i = 0; i < rows * cols; i++) {
	    mat.data[i] = min + (max - min) * ((float)rand() / (float)RAND_MAX);
	}

    return mat;
}

// Function to generate a new Matrix with random values
Matrix matrand(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (float *)malloc(rows * cols * sizeof(float));

    if (mat.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    float min = 0.0;
    float max = 1.0;
    // llelize the random number generation
    #pragma omp llel
    {
        // Each thread will have its own random seed
        unsigned int seed = omp_get_thread_num(); // Use thread number as seed

        #pragma omp for
        for (int i = 0; i < rows * cols; i++) {
            // Generate a random number using the thread-local seed
            mat.data[i] = min + (max - min) * ((float)rand_r(&seed) / (float)RAND_MAX);
        }
    }

    return mat;
}

// Function to perform matrix element-wise multiplication
Vector vecelem(const Vector *A, const Vector *B) {
    // Check if both vectors have the same dimensions
    if (A->size != B->size) {
        fprintf(stderr, "Vectors must have the same dimensions for element-wise multiplication.\n");
        exit(EXIT_FAILURE);
    }

    Vector C;
    C.size = A->size; // Result has the same dimensions as A and B
    C.data = (float *)malloc(C.size * sizeof(float));

    if (C.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    int i;
    // Perform element-wise multiplication
	for (i = 0; i < C.size; i++) {
        C.data[i] = A->data[i] * B->data[i];
        //printf("Thread %d processing index %d\n", omp_get_thread_num(), i);
    }

    return C;
}

// Function to perform matrix element-wise multiplication in llel
Vector vecelem(const Vector *A, const Vector *B) {
    // Check if both vectors have the same dimensions
    if (A->size != B->size) {
        fprintf(stderr, "Vectors must have the same dimensions for element-wise multiplication.\n");
        exit(EXIT_FAILURE);
    }

    Vector C;
    C.size = A->size; // Result has the same dimensions as A and B
    C.data = (float *)malloc(C.size * sizeof(float));

    if (C.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    int i;
    // Perform element-wise multiplication using OpenMP for llelization
    #pragma omp llel for default(none) shared(A, B, C) private(i)	
	for (i = 0; i < C.size; i++) {
        C.data[i] = A->data[i] * B->data[i];
        //printf("Thread %d processing index %d\n", omp_get_thread_num(), i);
    }

    return C;
}

// Function to perform matrix element-wise multiplication
Matrix matelem(const Matrix *A, const Matrix *B) {
    // Check if both matrices have the same dimensions
    if (A->rows != B->rows || A->cols != B->cols) {
        fprintf(stderr, "Matrices must have the same dimensions for element-wise multiplication.\n");
        exit(EXIT_FAILURE);
    }

    Matrix C;
    C.rows = A->rows;
    C.cols = A->cols; // Result has the same dimensions as A and B
    C.data = (float *)malloc(C.rows * C.cols * sizeof(float));

    if (C.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    int i, j;
    // Perform element-wise multiplication 
    // Use OpenMP for llelization for HUGE matrices
    #pragma omp llel for default(none) shared(A, B, C) private(i, j) collapse(2)
    for (i = 0; i < C.rows; i++) {
        for (j = 0; j < C.cols; j++) {
            C.data[i * C.cols + j] = A->data[i * A->cols + j] * B->data[i * B->cols + j];
        }
    }

    return C;
}
*/
int main(){
    /*
    Importrant notice on lellization:
    
    Possibly Lost Memory Reasons:
    In case of "possibly lost", memory it is likely allocated by the OpenMP runtime 
    for managing threads, and it is not directly freed when your program exits.
    
    Thread Management:
    The memory allocated for thread management (like thread-local storage) by the OpenMP 
    library is not freed until the program terminates. This is typical behavior for many 
    threading libraries, and it does not usually indicate a problem with your code.
    */
    //omp_set_num_threads(8);

    struct timespec start, stop;

	float vec1Values[] = {9.0, 8.0, 7.0, 6.0};
	float vec2Values[] = {1.0, 2.0, 3.0, 4.0};
	int dim = sizeof(vec1Values) / sizeof(vec1Values[0]);
	//printf("Length of vector: %d\n", dim);

	Vector vec1 = createVector(dim, vec1Values);
	Vector vec2 = createVector(dim, vec2Values);
	Vector vec = vecadd(&vec1, &vec2);
	printf("Vector addition result:\n");
	printVector(&vec);

    //vec = vecelem(&vec1, &vec2);
    //printf("Vector element-wise multiplication result:\n");
    //printVector(&vec);
    long int DIMEN = 1000;
    Vector bigVec1 = vecrand(DIMEN);
    Vector bigVec2 = vecrand(DIMEN);

    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);
    Vector bigVec = vecelem(&bigVec1, &bigVec2);
    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);
    // Calculate the elapsed time in seconds
    double elapsed_time = (stop.tv_sec - start.tv_sec) * 1e9;
    elapsed_time = (elapsed_time + (stop.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("\nVector %ld element-wise multiplication result:\n", DIMEN);
    printVector(&bigVec);
    printf("Vector multiplication took %.3e seconds.\n", elapsed_time);

    long int DIMEN2 = 1000;
    Matrix mat1 = matrandPara(DIMEN2, DIMEN2);
    Matrix mat2 = matrandPara(DIMEN2, DIMEN2);

    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);

    Matrix mat  = matelemPara(&mat1, &mat2);

    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);

    // Calculate the elapsed time in seconds
    elapsed_time = (stop.tv_sec - start.tv_sec) * 1e9;
    elapsed_time = (elapsed_time + (stop.tv_nsec - start.tv_nsec)) * 1e-9;

    printf("\nMatrix %ld element-wise multiplication result:\n", DIMEN2);
    printMatrix(&mat);
    printf("Matrix multiplication took %.3e seconds.\n", elapsed_time);

    Matrix matT = transposeMatrixPara(&mat);
    printf("\nTransposed matrix result:\n");
    printMatrix(&matT);

	freeVector(&vec1);
	freeVector(&vec2);
	freeVector(&vec );
	freeVector(&bigVec1);
	freeVector(&bigVec2);
	freeVector(&bigVec);
    freeMatrix(&mat1);
    freeMatrix(&mat2);
    freeMatrix(&mat);
    freeMatrix(&matT);

	return 0;
}
