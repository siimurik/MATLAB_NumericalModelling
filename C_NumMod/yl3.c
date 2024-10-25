// Compile and execute with:
// gcc -o yl3 yl3.c splash.c -llapacke -lblas -lm
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include <time.h>
#include "splash.h"

// Define clock default starting value
#define CLOCK_MONOTONIC 1

// Function pointer type for the iteration function
typedef float (*iterFunc)(float);

float g (float x)
{
    return cbrt(8.0*x - 2.0); // cbrt - cubreroot
}

// The void C function
void him(iterFunc g, float x_init, float tol, int max_iter, float *x_new, int *count) {

    // Initialize variables
    *x_new = x_init;
    float x_old = x_init - 1.0; // Initial value different from x_new
    *count = 0;

    // Iteration process
    while (fabs(*x_new - x_old) >= tol && *count < max_iter) {
        x_old = *x_new;
        *x_new = g(x_old);
        (*count)++;
    }

    // Check if maximum iteration count was exceeded
    if (*count >= max_iter) {
        printf("Warning: Maximum iteration count reachedd.\n");
    }

}

int main()
{
	// Finding the roots of a polynomial x^3 − 8^x + 2 = 0
	float r[4] = {1.0, 0.0, 8.0, 2.0};
	int dim = sizeof(r) / sizeof(r[0]); 
	Matrix matCompan = compan(r, dim);
	printf("Companion matrix:\n");
	printMatrix(&matCompan);
	
	// The eigenvalues are the roots to this polynomial
    printf("\n");
	eig(&matCompan);
	freeMatrix(&matCompan);

    // Harilik iteratsioonimeetod
    float x_init = -3.0;
    float tol    = 1E-4;
    int max_iter = 1000; 
    int count;
    float x_new;
    him(g, x_init, tol, max_iter, &x_new, &count);
    printf("\nSolution: %f, Iterations: %d\n", x_new, count);
    ////////////////////////////////////////////////////////
    // Additional testing
    //omp_set_num_threads(8);
    struct timespec start, stop;
    float vec1Values[] = {9.0, 8.0, 7.0, 6.0};
	float vec2Values[] = {1.0, 2.0, 3.0, 4.0};
	int length = sizeof(vec1Values) / sizeof(vec1Values[0]);
	//printf("Length of vector: %d\n", length);

    // Making the vectors
	Vector vec1 = createVector(length, vec1Values);
	Vector vec2 = createVector(length, vec2Values);
	Vector vec = vecadd(&vec1, &vec2);
	printf("\nAdditional testing below.\nVector addition result:\n");
	printVector(&vec);

    //vec = vecelem(&vec1, &vec2);
    //printf("Vector element-wise multiplication result:\n");
    //printVector(&vec);
    int DIMEN = 1000;
    Vector bigVec1 = vecrand(DIMEN);
    Vector bigVec2 = vecrand(DIMEN);

    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);
    // Test vector element-wise multiplication
    Vector bigVec = vecelem(&bigVec1, &bigVec2);
    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);
    // Calculate the elapsed time in seconds
    double elapsed_time = (stop.tv_sec - start.tv_sec) * 1e9;
    elapsed_time = (elapsed_time + (stop.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("\nVector %d element-wise multiplication result:\n", DIMEN);
    printVector(&bigVec);
    printf("Matrix multiplication took %.3e seconds.\n", elapsed_time);

    freeVector(&vec1);
	freeVector(&vec2);
	freeVector(&vec );
	freeVector(&bigVec1);
	freeVector(&bigVec2);
	freeVector(&bigVec);

	return 0;
}
