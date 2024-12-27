//  gcc gsl_debug.c -o deb -lgsl -lm && ./deb
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

typedef struct {
    gsl_matrix *gsl_matrix_ptr;
    unsigned int rows;
    unsigned int cols;
} Matrix;

typedef struct {
    gsl_vector *gsl_vector_ptr;
    unsigned int size;
} Vector;

void freeMatrix(Matrix *matrix) {
    if (matrix->gsl_matrix_ptr != NULL) {
        gsl_matrix_free(matrix->gsl_matrix_ptr);
        matrix->gsl_matrix_ptr = NULL;
    }
}

void freeVector(Vector *vector) {
    if (vector->gsl_vector_ptr != NULL) {
        gsl_vector_free(vector->gsl_vector_ptr); // Free the GSL vector
        vector->gsl_vector_ptr = NULL;           // Avoid dangling pointer
    }
}

void printMatrix(const Matrix *mat) {
    printf("[\n");
    for (unsigned int i = 0; i < mat->rows; i++) {
        printf("  [");
        for (unsigned int j = 0; j < mat->cols; j++) {
            printf("%7.4f", gsl_matrix_get(mat->gsl_matrix_ptr, i, j));
            if (j < mat->cols - 1) {
                printf(", ");
            }
        }
        printf("]\n");
    }
    printf("]\n");
}

void printVector(const Vector *vector) {
    printf("[");

    if (vector->size > 10) {
        // Print the first three elements
        for (int i = 0; i < 3; i++) {
            printf("%.4f", gsl_vector_get(vector->gsl_vector_ptr, i));
            if (i < 2) {
                printf(", ");
            }
        }
        // Print ellipsis
        printf(", ...");
        // Print the last three elements
        for (int i = vector->size - 3; i < vector->size; i++) {
            printf(", %.4f", gsl_vector_get(vector->gsl_vector_ptr, i));
        }
    } else {
        // Print all elements if size is 10 or less
        for (int i = 0; i < vector->size; i++) {
            printf("%.4f", gsl_vector_get(vector->gsl_vector_ptr, i));
            if (i < vector->size - 1) {
                printf(", ");
            }
        }
    }

    printf("]\n");
}

Vector createVector(int size, double *values) {
    Vector vector;
    vector.size = size;
    vector.gsl_vector_ptr = gsl_vector_alloc(size);

    if (vector.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Populate the vector with the provided values
    for (int i = 0; i < size; i++) {
        gsl_vector_set(vector.gsl_vector_ptr, i, values[i]);
    }

    return vector;
}

void qr_decomposition(Matrix *A, Matrix *Q, Matrix *R) {
    gsl_matrix *temp_R = gsl_matrix_alloc(A->rows, A->cols);
    gsl_vector *tau = gsl_vector_alloc(A->cols);  // Use A->cols instead of A->rows

    gsl_matrix_memcpy(temp_R, A->gsl_matrix_ptr);

    gsl_linalg_QR_decomp(temp_R, tau);

    gsl_matrix *Q_temp = gsl_matrix_alloc(A->rows, A->rows); // Temporary Q matrix

    gsl_linalg_QR_unpack(temp_R, tau, Q_temp, R->gsl_matrix_ptr);

    // Copy the relevant part of Q_temp to Q
    for (unsigned int i = 0; i < A->rows; i++) {
        for (unsigned int j = 0; j < A->cols; j++) {
            gsl_matrix_set(Q->gsl_matrix_ptr, i, j, gsl_matrix_get(Q_temp, i, j));
        }
    }

    gsl_matrix_free(temp_R);
    gsl_vector_free(tau);
    gsl_matrix_free(Q_temp);
}

Matrix trimMatrix(const Matrix *input, unsigned int new_rows, unsigned int new_cols) {
    Matrix output;
    output.rows = new_rows;
    output.cols = new_cols;
    output.gsl_matrix_ptr = gsl_matrix_alloc(output.rows, output.cols);

    for (unsigned int i = 0; i < new_rows; i++) {
        for (unsigned int j = 0; j < new_cols; j++) {
            gsl_matrix_set(output.gsl_matrix_ptr, i, j, gsl_matrix_get(input->gsl_matrix_ptr, i, j));
        }
    }
    return output;
}

Matrix transposeMatrix(const Matrix *mat) {
    Matrix transposed;
    transposed.rows = mat->cols;
    transposed.cols = mat->rows;
    transposed.gsl_matrix_ptr = gsl_matrix_alloc(transposed.rows, transposed.cols);

    if (transposed.gsl_matrix_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    int i, j;
    // Perform transposition
    for (i = 0; i < mat->rows; i++) {
        for (j = 0; j < mat->cols; j++) {
            gsl_matrix_set(transposed.gsl_matrix_ptr, j, i, gsl_matrix_get(mat->gsl_matrix_ptr, i, j));
        }
    }

    return transposed;
}

Vector matvec(const Matrix *A, const Vector *x) {
    // Check if the number of columns in A matches the size of the vector x
    if (A->cols != x->size) {
        fprintf(stderr, "Error: Number of columns in matrix A must match the size of vector x.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the result vector
    Vector result;
    result.size = A->rows;
    result.gsl_vector_ptr = gsl_vector_alloc(result.size);
    if (result.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Perform matrix-vector multiplication using GSL
    gsl_blas_dgemv(CblasNoTrans, 1.0, A->gsl_matrix_ptr, x->gsl_vector_ptr, 0.0, result.gsl_vector_ptr);

    return result;
}

Vector polyfitweighted(const Vector *x, const Vector *y, const Vector *w, unsigned int n) {
    unsigned int len = x->size;
    unsigned int cols = n + 1;

    Matrix V;
    V.rows = len;
    V.cols = cols;
    V.gsl_matrix_ptr = gsl_matrix_alloc(V.rows, V.cols);
    Vector wy;
    wy.size = len;
    wy.gsl_vector_ptr = gsl_vector_alloc(wy.size);

    if (V.gsl_matrix_ptr == NULL || wy.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    // Construct the weighted Vandermonde matrix
    for (unsigned int i = 0; i < len; i++) {
        // Set the last column as the weight
        gsl_matrix_set(V.gsl_matrix_ptr, i, n, gsl_vector_get(w->gsl_vector_ptr, i));

        // Fill the Vandermonde row starting from the highest power to the lowest
        for (int j = n - 1; j >= 0; j--) {
            gsl_matrix_set(V.gsl_matrix_ptr, i, j, gsl_vector_get(x->gsl_vector_ptr, i) * gsl_matrix_get(V.gsl_matrix_ptr, i, j + 1));
        }

        // Calculate weighted y values
        gsl_vector_set(wy.gsl_vector_ptr, i, gsl_vector_get(w->gsl_vector_ptr, i) * gsl_vector_get(y->gsl_vector_ptr, i));
    }

    // Perform SVD
    gsl_matrix *U = gsl_matrix_alloc(V.rows, V.cols);
    gsl_matrix *V_svd = gsl_matrix_alloc(V.cols, V.cols);
    gsl_vector *S = gsl_vector_alloc(V.cols);
    gsl_vector *work = gsl_vector_alloc(V.cols);

    gsl_matrix_memcpy(U, V.gsl_matrix_ptr);
    gsl_linalg_SV_decomp(U, V_svd, S, work);

    // Solve for p using SVD
    gsl_vector *p = gsl_vector_alloc(V.cols);
    gsl_linalg_SV_solve(U, V_svd, S, wy.gsl_vector_ptr, p);

    Vector p_result;
    p_result.size = p->size;
    p_result.gsl_vector_ptr = p;

    // Free the allocated matrices and vectors used for SVD
    gsl_matrix_free(U);
    gsl_matrix_free(V_svd);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_matrix_free(V.gsl_matrix_ptr);
    gsl_vector_free(wy.gsl_vector_ptr);

    return p_result;
}

int main() {
    double xData[] = {1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
    double yData[] = {9.08, 10.43, 11.9, 13.48, 15.19, 17.03, 19.01, 21.13, 23.39};
    double wData[] = {1.0, 1.0, 2.0, 5.0, 1.0, 4.0, 2.0, 2.0, 1.0};

    Vector x = createVector(sizeof(xData) / sizeof(xData[0]), xData);
    Vector y = createVector(sizeof(yData) / sizeof(yData[0]), yData);
    Vector w = createVector(sizeof(wData) / sizeof(wData[0]), wData);
    unsigned int n = 3;  // Degree of polynomial

    Vector p = polyfitweighted(&x, &y, &w, n);
    
    printf("\np = ");
    printVector(&p);

    freeVector(&x);
    freeVector(&y);
    freeVector(&w);
    freeVector(&p);

    return 0;
}
