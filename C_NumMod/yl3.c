// Compile and execute with:
// gcc -o yl3 yl3.c splash.c -llapacke -lblas
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <lapacke.h>
#include <cblas.h>
#include "splash.h"
/*
typedef struct{
    float *data;
    int    rows;
    int    cols;
} Matrix;

// Function to create the companion matrix:
Matrix compan(const float *coefficients, int size) {
    // Input validation
    if (size < 2) {
        fprintf(stderr, "The input array must have at least 2 elements (degree 1 polynomial).\n");
        return (Matrix){NULL, 0, 0}; // Return an empty matrix
    }

    if (coefficients[0] == 0) {
        fprintf(stderr, "The first coefficient in the array must not be zero.\n");
        return (Matrix){NULL, 0, 0}; // Return an empty matrix
    }

    // Create the companion matrix
    int n = size - 1; // Degree of the polynomial
    Matrix C;
    C.rows = n;
    C.cols = n;
	// Use calloc to initialize all values to 0
    C.data = (float *)calloc(C.rows * C.cols, sizeof(float));

    if (C.data == NULL) {
        fprintf(stderr, "Memory allocation failed for companion matrix.\n");
        return (Matrix){NULL, 0, 0}; // Return an empty matrix
    }

    // Fill the first row
    float leading_coefficient = coefficients[0];
    for (int j = 0; j < n; j++) {
        C.data[0 * C.cols + j] = -coefficients[j + 1] / leading_coefficient; // First row
    }

    // Fill the sub-diagonal with ones
    for (int i = 1; i < n; i++) {
        C.data[i * C.cols + (i - 1)] = 1.0;
    }

    // Fill the rest with zeros (already initialized to 0)
    //for (int i = 1; i < n; i++) {
    //    for (int j = i + 1; j < n; j++) {
    //        C.data[i * C.cols + j] = 0.0; // Explicitly set to 0 (optional)
    //    }
    //}

    return C;
}

// Eigenvalue Calculation Function
void eig(const Matrix *matrix) {
    if (matrix->rows != matrix->cols) {
        fprintf(stderr, "The matrix must be square to compute eigenvalues.\n");
        return;
    }

    char  JOBVL = 'N';
    char  JOBVR = 'V';
    int   N     = matrix->rows;
    int   LDA   = matrix->cols;
    float *A    = (float *)malloc(LDA * N * sizeof(float));
    float *WR   = (float *)malloc(N * sizeof(float)); // Real parts of eigenvalues
    float *WI   = (float *)malloc(N * sizeof(float)); // Imaginary parts of eigenvalues
    int   LDVL  = matrix->rows;
    int   LDVR  = matrix->rows;
    float *VR   = (float *)malloc(LDVR * N * sizeof(float)); // Right eigenvectors
    int INFO;

    if (A == NULL || WR == NULL || WI == NULL || VR == NULL) {
        fprintf(stderr, "Failed to allocate memory.\n");
        goto cleanup;
    }

    // Copy matrix data to LAPACK format
    memcpy(A, matrix->data, N * N * sizeof(float));

    // Compute eigenvalues and right eigenvectors
    //INFO = LAPACKE_sgeev(LAPACK_COL_MAJOR, JOBVL, JOBVR, N, A, LDA, WR, WI, NULL, N, VR, N);
	INFO = LAPACKE_sgeev( LAPACK_COL_MAJOR, JOBVL, JOBVR, N, A, LDA, WR, WI, NULL, LDVL, VR, LDVR);

    if (INFO > 0) {
        fprintf(stderr, "Error in eigenvalue computation: %d\n", INFO);
        goto cleanup;
    }

    // Print eigenvalues
    printf("Eigenvalues:\n");
    for (int i = 0; i < N; i++) {
        if (WI[i] == 0) {
            // Real eigenvalue
            printf("   %7.4f + %7.4fi\n", WR[i], WI[i]);
        } else {
            // Complex eigenvalue
            printf("   %7.4f + %7.4fi\n", WR[i], WI[i]);
            printf("   %7.4f - %7.4fi\n", WR[i], WI[i]);
            i++; // Skip the next one as it's a pair
        }
    }

cleanup:
    free(A);
    free(WR);
    free(WI);
    free(VR);
}


// Function for printing a matrix
void printMatrix(const Matrix *mat) {
    int max_print_size = 6;

    printf("[\n");

    if (mat->rows <= max_print_size) {
        for (int i = 0; i < mat->rows; i++) {
            printf("  [");
            for (int j = 0; j < mat->cols; j++) {
                printf("%9.4f", mat->data[i * mat->cols + j]);
                if (j < mat->cols - 1) printf(", ");
            }
            printf("]");
            if (i < mat->rows - 1) printf(",\n");
        }
    } else {
        for (int i = 0; i < max_print_size; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size; j++) {
                printf("%9.4f", mat->data[i * mat->cols + j]);
                if (j < mat->cols - 1) printf(", ");
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

void freeMatrix(Matrix *matrix) {
    free(matrix->data);
    matrix->data = NULL;    // Avoid dangling pointer.
} 
*/
int main()
{
	// Finding the roots of a polynomial x^3 − 8^x + 2 = 0
	float r[4] = {1.0, 0.0, 8.0, 2.0};
	int dim = sizeof(r) / sizeof(r[0]); 
	Matrix matCompan = compan(r, dim);
	printf("Companion matrix:\n");
	printMatrix(&matCompan);
	
	// The eigenvalues are the roots to this polynomial
	eig(&matCompan);

	freeMatrix(&matCompan);
	return 0;
}
