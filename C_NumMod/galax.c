/******************************************************/
//  GALAX - GSL Advanced Analytics eXtension
/******************************************************/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_permutation.h>

typedef struct {
    gsl_matrix *gsl_matrix_ptr;
    unsigned int rows;
    unsigned int cols;
} Matrix;

typedef struct {
    gsl_vector *gsl_vector_ptr;
    unsigned int size;
} Vector;

/////////////////////////////////////////////////////////////////////////////

void freeMatrix(Matrix *matrix) {
    if (matrix->gsl_matrix_ptr != NULL) {
        gsl_matrix_free(matrix->gsl_matrix_ptr); // Free the GSL matrix
        matrix->gsl_matrix_ptr = NULL;           // Avoid dangling pointer
    }
    // No need to free the `matrix` pointer itself as it's often stack-allocated
}


Matrix zerosMatrix(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.gsl_matrix_ptr = gsl_matrix_calloc(rows, cols); // Allocate and initialize with zeros
    return mat;
}


// Function to create and initialize a matrix with specific values
Matrix createMatrix(int rows, int cols, double *values) {
    Matrix matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.gsl_matrix_ptr = gsl_matrix_alloc(rows, cols); // Allocate memory

    if (matrix.gsl_matrix_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Populate the matrix with the provided values
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            gsl_matrix_set(matrix.gsl_matrix_ptr, i, j, values[i * cols + j]);
        }
    }

    return matrix;
}

// Function for printing a matrix
void printMatrix(const Matrix *mat) {
    int max_print_size = 3; // Print 3 rows and 3 columns from corners

    printf("[\n");

    if (mat->rows <= 10 && mat->cols <= 10) {
        // Print the entire matrix if both dimensions are 10 or less
        for (int i = 0; i < mat->rows; i++) {
            printf("  [");
            for (int j = 0; j < mat->cols; j++) {
                printf("%7.4f", gsl_matrix_get(mat->gsl_matrix_ptr, i, j));
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
                printf("%7.4f", gsl_matrix_get(mat->gsl_matrix_ptr, i, j));
                if (j < max_print_size - 1 && j < mat->cols - 1) printf(", ");
            }
            // Print ellipsis for omitted columns
            if (mat->cols > max_print_size) {
                printf(", ... ");
            }
            // Last elements of the top 3 rows
            for (int j = mat->cols - 3; j < mat->cols; j++) {
                if (j >= 0) { // Ensure we don't go out of bounds
                    printf("%7.4f", gsl_matrix_get(mat->gsl_matrix_ptr, i, j));
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
                printf("%7.4f", gsl_matrix_get(mat->gsl_matrix_ptr, i, j));
                if (j < max_print_size - 1 && j < mat->cols - 1) printf(", ");
            }
            // Print ellipsis for omitted columns
            if (mat->cols > max_print_size) {
                printf(", ... ");
            }
            // Last elements of the bottom 3 rows
            for (int j = mat->cols - 3; j < mat->cols; j++) {
                if (j >= 0) { // Ensure we don't go out of bounds
                    printf("%7.4f", gsl_matrix_get(mat->gsl_matrix_ptr, i, j));
                    if (j < mat->cols - 1) printf(", ");
                }
            }
            printf("]");
            if (i < mat->rows - 1) printf(",\n");
        }
    }

    printf("\n]\n");
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

Matrix inverseMatrix(const Matrix *mat) {
    if (mat->rows != mat->cols) {
        fprintf(stderr, "Matrix must be square to compute the inverse.\n");
        exit(EXIT_FAILURE);
    }

    int n = mat->rows;
    Matrix inverse;
    inverse.rows = n;
    inverse.cols = n;
    inverse.gsl_matrix_ptr = gsl_matrix_alloc(n, n);

    if (inverse.gsl_matrix_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Copy the original matrix to the inverse matrix
    gsl_matrix_memcpy(inverse.gsl_matrix_ptr, mat->gsl_matrix_ptr);

    gsl_permutation *p = gsl_permutation_alloc(n);
    int signum;

    // Perform LU decomposition
    gsl_linalg_LU_decomp(inverse.gsl_matrix_ptr, p, &signum);

    // Compute the inverse
    gsl_matrix *inverse_gsl = gsl_matrix_alloc(n, n);
    gsl_linalg_LU_invert(inverse.gsl_matrix_ptr, p, inverse_gsl);
    gsl_matrix_memcpy(inverse.gsl_matrix_ptr, inverse_gsl);

    gsl_matrix_free(inverse_gsl);
    gsl_permutation_free(p);

    return inverse;
}

Matrix matmul(const Matrix *A, const Matrix *B) {
    if (A->cols != B->rows) {
        fprintf(stderr, "Matrix dimensions are incompatible for multiplication.\n");
        exit(EXIT_FAILURE);
    }

    Matrix C;
    C.rows = A->rows;
    C.cols = B->cols;
    C.gsl_matrix_ptr = gsl_matrix_alloc(C.rows, C.cols);

    // Perform matrix multiplication
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A->gsl_matrix_ptr, B->gsl_matrix_ptr, 0.0, C.gsl_matrix_ptr);

    return C;
}

Matrix matadd(const Matrix *A, const Matrix *B) {
    if (A->rows != B->rows || A->cols != B->cols) {
        fprintf(stderr, "Matrix dimensions are incompatible for addition.\n");
        exit(EXIT_FAILURE);
    }

    Matrix C;
    C.rows = A->rows;
    C.cols = A->cols;
    C.gsl_matrix_ptr = gsl_matrix_alloc(C.rows, C.cols);

    // Copy B into C
    gsl_matrix_memcpy(C.gsl_matrix_ptr, B->gsl_matrix_ptr);

    // Perform matrix addition
    gsl_matrix_add(C.gsl_matrix_ptr, A->gsl_matrix_ptr);

    return C;
}

Matrix matscal(const Matrix *A, double scalar) {
    Matrix C;
    C.rows = A->rows;
    C.cols = A->cols;
    C.gsl_matrix_ptr = gsl_matrix_alloc(C.rows, C.cols);

    // Copy A into C
    gsl_matrix_memcpy(C.gsl_matrix_ptr, A->gsl_matrix_ptr);

    // Scale the matrix by a constant
    gsl_matrix_scale(C.gsl_matrix_ptr, scalar);

    return C;
}

double det(const Matrix *mat) {
    // Check if the matrix is square
    if (mat->rows != mat->cols) {
        fprintf(stderr, "Matrix must be square to compute the determinant.\n");
        exit(EXIT_FAILURE);
    }

    int n = mat->rows;
    gsl_matrix *tmp = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(tmp, mat->gsl_matrix_ptr);

    gsl_permutation *p = gsl_permutation_alloc(n);
    int signum;

    // Perform LU decomposition
    gsl_linalg_LU_decomp(tmp, p, &signum);

    // Compute the determinant
    double determinant = gsl_linalg_LU_det(tmp, signum);

    gsl_matrix_free(tmp);
    gsl_permutation_free(p);

    return determinant;
}

Matrix matelem(const Matrix *A, const Matrix *B) {
    // Check if both matrices have the same dimensions
    if (A->rows != B->rows || A->cols != B->cols) {
        fprintf(stderr, "Matrices must have the same dimensions for element-wise multiplication.\n");
        exit(EXIT_FAILURE);
    }

    Matrix C;
    C.rows = A->rows;
    C.cols = A->cols;
    C.gsl_matrix_ptr = gsl_matrix_alloc(C.rows, C.cols);

    // Perform element-wise multiplication
    for (int i = 0; i < C.rows; i++) {
        for (int j = 0; j < C.cols; j++) {
            double value = gsl_matrix_get(A->gsl_matrix_ptr, i, j) * gsl_matrix_get(B->gsl_matrix_ptr, i, j);
            gsl_matrix_set(C.gsl_matrix_ptr, i, j, value);
        }
    }

    return C;
}

Matrix matrand(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.gsl_matrix_ptr = gsl_matrix_alloc(rows, cols);

    if (mat.gsl_matrix_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_matrix *temp = gsl_matrix_alloc(rows, cols);

    gsl_matrix_set_all(temp, 1.0);  // Fill with 1.0 for scaling
    gsl_matrix_scale(temp, gsl_rng_uniform(rng));  // Scale with random values

    gsl_matrix_memcpy(mat.gsl_matrix_ptr, temp);

    gsl_matrix_free(temp);
    gsl_rng_free(rng);

    return mat;
}

Matrix compan(const Vector *coefficients) {
    int size = coefficients->size;

    // Input validation
    if (size < 2) {
        fprintf(stderr, "The input vector must have at least 2 elements (degree 1 polynomial).\n");
        return (Matrix){NULL, 0, 0}; // Return an empty matrix
    }

    if (gsl_vector_get(coefficients->gsl_vector_ptr, 0) == 0) {
        fprintf(stderr, "The first coefficient in the vector must not be zero.\n");
        return (Matrix){NULL, 0, 0}; // Return an empty matrix
    }

    // Create the companion matrix
    int n = size - 1; // Degree of the polynomial
    Matrix C;
    C.rows = n;
    C.cols = n;
    C.gsl_matrix_ptr = gsl_matrix_calloc(C.rows, C.cols);

    if (C.gsl_matrix_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for companion matrix.\n");
        return (Matrix){NULL, 0, 0}; // Return an empty matrix
    }

    // Fill the first row
    double leading_coefficient = gsl_vector_get(coefficients->gsl_vector_ptr, 0);
    for (int j = 0; j < n; j++) {
        gsl_matrix_set(C.gsl_matrix_ptr, 0, j, -gsl_vector_get(coefficients->gsl_vector_ptr, j + 1) / leading_coefficient);
    }

    // Fill the sub-diagonal with ones
    for (int i = 1; i < n; i++) {
        gsl_matrix_set(C.gsl_matrix_ptr, i, i - 1, 1.0);
    }

    return C;
}


void eig(const Matrix *matrix) {
    if (matrix->rows != matrix->cols) {
        fprintf(stderr, "The matrix must be square to compute eigenvalues.\n");
        exit(EXIT_FAILURE);
    }

    int n = matrix->rows;
    Matrix result;
    result.rows = n;
    result.cols = 2; // To store real and imaginary parts
    result.gsl_matrix_ptr = gsl_matrix_alloc(n, 2);

    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(n);
    gsl_vector_complex *eval = gsl_vector_complex_alloc(n);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(n, n);
    gsl_matrix *temp = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(temp, matrix->gsl_matrix_ptr);

    gsl_eigen_nonsymmv(temp, eval, evec, w);

    gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        gsl_complex z = gsl_vector_complex_get(eval, i);
        double real = GSL_REAL(z);
        double imag = GSL_IMAG(z);
        gsl_matrix_set(result.gsl_matrix_ptr, i, 0, real);
        gsl_matrix_set(result.gsl_matrix_ptr, i, 1, imag);
        printf("%7.4f + %7.4fi\n", real, imag);
    }

    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_matrix_free(temp);
    gsl_eigen_nonsymmv_free(w);
    freeMatrix(&result);
    //return result;
}

Matrix eigB(const Matrix *matrix, bool save_result) {
    if (matrix->rows != matrix->cols) {
        fprintf(stderr, "The matrix must be square to compute eigenvalues.\n");
        exit(EXIT_FAILURE);
    }

    int n = matrix->rows;
    Matrix result;
    result.rows = n;
    result.cols = 2; // To store real and imaginary parts
    result.gsl_matrix_ptr = gsl_matrix_alloc(n, 2);

    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(n);
    gsl_vector_complex *eval = gsl_vector_complex_alloc(n);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(n, n);
    gsl_matrix *temp = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(temp, matrix->gsl_matrix_ptr);

    gsl_eigen_nonsymmv(temp, eval, evec, w);

    gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    if (save_result == true){
        printf("Eigenvalues:\n");
        for (int i = 0; i < n; i++) {
            gsl_complex z = gsl_vector_complex_get(eval, i);
            double real = GSL_REAL(z);
            double imag = GSL_IMAG(z);
            gsl_matrix_set(result.gsl_matrix_ptr, i, 0, real);
            gsl_matrix_set(result.gsl_matrix_ptr, i, 1, imag);
            printf("%7.4f + %7.4fi\n", real, imag);
        }
        gsl_vector_complex_free(eval);
        gsl_matrix_complex_free(evec);
        gsl_matrix_free(temp);
        gsl_eigen_nonsymmv_free(w);
        freeMatrix(&result);
    } else {
        return result;
        gsl_vector_complex_free(eval);
        gsl_matrix_complex_free(evec);
        gsl_matrix_free(temp);
        gsl_eigen_nonsymmv_free(w);
    }
}


double normMatrix(const Matrix *mat, double p) {
    double norm_value = 0.0;

    if (p == 1) {
        // 1-norm (maximum column sum)
        for (int j = 0; j < mat->cols; j++) {
            double col_sum = 0.0;
            for (int i = 0; i < mat->rows; i++) {
                col_sum += fabs(gsl_matrix_get(mat->gsl_matrix_ptr, i, j));
            }
            if (col_sum > norm_value) {
                norm_value = col_sum;
            }
        }
    } else if (p == 2) {
        // Frobenius norm (sum of squares of all elements, then square root)
        double sum = 0.0;
        for (int i = 0; i < mat->rows; i++) {
            for (int j = 0; j < mat->cols; j++) {
                double value = gsl_matrix_get(mat->gsl_matrix_ptr, i, j);
                sum += value * value;
            }
        }
        norm_value = sqrt(sum);
    } else if (p == INFINITY) {
        // Infinity-norm (maximum row sum)
        for (int i = 0; i < mat->rows; i++) {
            double row_sum = 0.0;
            for (int j = 0; j < mat->cols; j++) {
                row_sum += fabs(gsl_matrix_get(mat->gsl_matrix_ptr, i, j));
            }
            if (row_sum > norm_value) {
                norm_value = row_sum;
            }
        }
    } else {
        fprintf(stderr, "Unsupported norm type for matrices.\n");
        exit(EXIT_FAILURE);
    }

    return norm_value;
}




Matrix linsolve_overdet(const Matrix *A, const Matrix *F) {
    // Ensure the dimensions are compatible
    if (A->rows != F->rows) {
        fprintf(stderr, "Matrix dimensions are incompatible for solving.\n");
        exit(EXIT_FAILURE);
    }

    int n = A->rows;  // Number of equations
    int m = A->cols;  // Number of variables
    int nrhs = F->cols;  // Number of right-hand sides

    Matrix x;
    x.rows = m;
    x.cols = nrhs;
    x.gsl_matrix_ptr = gsl_matrix_alloc(m, nrhs);

    gsl_matrix *gslA = gsl_matrix_alloc(n, m);
    gsl_vector *gslB = gsl_vector_alloc(n); // Vector instead of matrix for y (right-hand side)
    gsl_matrix *cov = gsl_matrix_alloc(m, m);
    gsl_vector *c = gsl_vector_alloc(m);
    double chisq;

    gsl_matrix_memcpy(gslA, A->gsl_matrix_ptr);

    // We need to solve for each column of F separately
    for (int j = 0; j < nrhs; j++) {
        for (int i = 0; i < n; i++) {
            gsl_vector_set(gslB, i, gsl_matrix_get(F->gsl_matrix_ptr, i, j));
        }

        gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, m);

        // Solve the system using GSL's least squares solver
        gsl_multifit_linear(gslA, gslB, c, cov, &chisq, work);

        // Copy the results to the output matrix
        for (int i = 0; i < m; i++) {
            gsl_matrix_set(x.gsl_matrix_ptr, i, j, gsl_vector_get(c, i));
        }

        gsl_multifit_linear_free(work);
    }

    gsl_matrix_free(gslA);
    gsl_vector_free(gslB);
    gsl_matrix_free(cov);
    gsl_vector_free(c);

    return x;
}


////////////////////////////////////////////////////////////////////////////

void freeVector(Vector *vector) {
    if (vector->gsl_vector_ptr != NULL) {
        gsl_vector_free(vector->gsl_vector_ptr); // Free the GSL vector
        vector->gsl_vector_ptr = NULL;           // Avoid dangling pointer
    }
    // No need to free the `vector` pointer itself as it's often stack-allocated
}


Vector zerosVector(int size) {
    Vector vec;
    vec.size = size;
    vec.gsl_vector_ptr = gsl_vector_calloc(size); // Allocate and initialize with zeros
    return vec;
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

Vector linsolve(const Matrix *A, const Vector *b) {
    // Check for dimension compatibility
    if (A == NULL || b == NULL || A->gsl_matrix_ptr == NULL || b->gsl_vector_ptr == NULL) {
        fprintf(stderr, "Error: Null pointer input.\n");
        exit(EXIT_FAILURE);
    }
    if (A->cols != b->size) {
        fprintf(stderr, "Error: Incompatible dimensions for system.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the result vector
    Vector x;
    x.size = b->size;
    x.gsl_vector_ptr = gsl_vector_alloc(x.size);
    if (x.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for result vector.\n");
        exit(EXIT_FAILURE);
    }

    gsl_matrix *tempA = gsl_matrix_alloc(A->rows, A->cols);
    gsl_vector *tempb = gsl_vector_alloc(b->size);
    gsl_matrix_memcpy(tempA, A->gsl_matrix_ptr);
    gsl_vector_memcpy(tempb, b->gsl_vector_ptr);

    gsl_permutation *p = gsl_permutation_alloc(A->rows);
    int signum;

    gsl_linalg_LU_decomp(tempA, p, &signum);
    gsl_linalg_LU_solve(tempA, p, tempb, x.gsl_vector_ptr);

    gsl_matrix_free(tempA);
    gsl_vector_free(tempb);
    gsl_permutation_free(p);

    return x; // Return the result vector
}

Vector vecadd(const Vector *A, const Vector *B) {
    if (A->size != B->size) {
        fprintf(stderr, "Vector dimensions are incompatible for addition.\n");
        exit(EXIT_FAILURE);
    }

    Vector C;
    C.size = A->size;
    C.gsl_vector_ptr = gsl_vector_alloc(C.size);
    if (C.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for result vector.\n");
        exit(EXIT_FAILURE);
    }

    // Copy B into C
    gsl_vector_memcpy(C.gsl_vector_ptr, B->gsl_vector_ptr);

    // Perform vector addition
    gsl_vector_add(C.gsl_vector_ptr, A->gsl_vector_ptr);

    return C;
}

Vector vecscal(const Vector *A, double scalar) {
    Vector C;
    C.size = A->size;
    C.gsl_vector_ptr = gsl_vector_alloc(C.size);

    // Copy A into C
    gsl_vector_memcpy(C.gsl_vector_ptr, A->gsl_vector_ptr);

    // Scale the vector by a constant
    gsl_vector_scale(C.gsl_vector_ptr, scalar);

    return C;
}

Vector vecrand(int dim) {
    Vector vec;
    vec.size = dim;
    vec.gsl_vector_ptr = gsl_vector_alloc(dim);

    if (vec.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    for (int i = 0; i < dim; i++) {
        gsl_vector_set(vec.gsl_vector_ptr, i, gsl_rng_uniform(rng));
    }

    gsl_rng_free(rng);

    return vec;
}

Vector vecelem(const Vector *A, const Vector *B) {
    // Check if both vectors have the same dimensions
    if (A->size != B->size) {
        fprintf(stderr, "Vectors must have the same dimensions for element-wise multiplication.\n");
        exit(EXIT_FAILURE);
    }

    Vector C;
    C.size = A->size;
    C.gsl_vector_ptr = gsl_vector_alloc(C.size);

    if (C.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Perform element-wise multiplication
    for (int i = 0; i < C.size; i++) {
        double value = gsl_vector_get(A->gsl_vector_ptr, i) * gsl_vector_get(B->gsl_vector_ptr, i);
        gsl_vector_set(C.gsl_vector_ptr, i, value);
    }

    return C;
}

double norm(const Vector *vec, double p) {
    double norm_value = 0.0;

    if (p == 1) {
        // 1-norm using GSL
        norm_value = gsl_blas_dasum(vec->gsl_vector_ptr);
    } else if (p == 2) {
        // 2-norm using GSL
        norm_value = gsl_blas_dnrm2(vec->gsl_vector_ptr);
    } else if (p == INFINITY) {
        // Infinity-norm (Maximum norm)
        for (int i = 0; i < vec->size; i++) {
            double abs_val = fabs(gsl_vector_get(vec->gsl_vector_ptr, i));
            if (abs_val > norm_value) {
                norm_value = abs_val;
            }
        }
    } else {
        // p-norm
        for (int i = 0; i < vec->size; i++) {
            norm_value += pow(fabs(gsl_vector_get(vec->gsl_vector_ptr, i)), p);
        }
        norm_value = pow(norm_value, 1.0 / p);
    }

    return norm_value;
}

Vector createArray(double start, double end, double step) {
    // Check if step size is valid
    if (step == 0.0) {
        fprintf(stderr, "Error: Step size cannot be zero.\n");
        exit(EXIT_FAILURE);
    }
    if ((end - start) / step < 0) {
        fprintf(stderr, "Error: Step size is inconsistent with start and end points.\n");
        exit(EXIT_FAILURE);
    }

    // Calculate the number of elements
    int size = (int)ceil((end - start) / step) + 1;

    // Allocate memory for the array
    Vector result;
    result.size = size;
    result.gsl_vector_ptr = gsl_vector_alloc(size);
    if (result.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for result array.\n");
        exit(EXIT_FAILURE);
    }

    // Fill the array with values
    double current = start;
    for (int i = 0; i < size; i++) {
        gsl_vector_set(result.gsl_vector_ptr, i, current);
        current += step;
        if ((step > 0 && current > end) || (step < 0 && current < end)) {
            current = end; // Ensure last value is exactly the endpoint
        }
    }

    return result;
}

Vector polyfitweighted(const Vector *x, const Vector *y, const Vector *w, int n) {
    int rows = x->size;
    int cols = n + 1;

    Matrix V;
    V.rows = rows;
    V.cols = cols;
    V.gsl_matrix_ptr = gsl_matrix_alloc(V.rows, V.cols);
    Vector wy;
    wy.size = rows;
    wy.gsl_vector_ptr = gsl_vector_alloc(wy.size);

    if (V.gsl_matrix_ptr == NULL || wy.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    // Construct the weighted Vandermonde matrix
    for (int i = 0; i < rows; i++) {
        gsl_matrix_set(V.gsl_matrix_ptr, i, n, gsl_vector_get(w->gsl_vector_ptr, i));

        for (int j = n - 1; j >= 0; j--) {
            gsl_matrix_set(V.gsl_matrix_ptr, i, j, gsl_vector_get(x->gsl_vector_ptr, i) * gsl_matrix_get(V.gsl_matrix_ptr, i, j + 1));
        }

        gsl_vector_set(wy.gsl_vector_ptr, i, gsl_vector_get(w->gsl_vector_ptr, i) * gsl_vector_get(y->gsl_vector_ptr, i));
    }

    gsl_matrix *Q = gsl_matrix_alloc(V.rows, V.cols);
    gsl_matrix *R = gsl_matrix_alloc(V.cols, V.cols);
    gsl_vector *tau = gsl_vector_alloc(V.cols);

    // Perform QR decomposition
    gsl_linalg_QR_decomp(V.gsl_matrix_ptr, tau);
    gsl_linalg_QR_unpack(V.gsl_matrix_ptr, tau, Q, R);

    gsl_vector *qwy = gsl_vector_alloc(V.cols);
    gsl_blas_dgemv(CblasTrans, 1.0, Q, wy.gsl_vector_ptr, 0.0, qwy);

    gsl_vector *p = gsl_vector_alloc(V.cols);
    gsl_linalg_R_solve(R, qwy, p);

    Vector result;
    result.size = V.cols;
    result.gsl_vector_ptr = p;

    gsl_matrix_free(Q);
    gsl_matrix_free(R);
    gsl_vector_free(tau);
    gsl_vector_free(qwy);
    gsl_matrix_free(V.gsl_matrix_ptr);
    gsl_vector_free(wy.gsl_vector_ptr);

    return result;
}

Vector linspace(double start, double end, int num) {
    Vector result;
    result.size = num;
    result.gsl_vector_ptr = gsl_vector_alloc(num);

    if (result.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    double step = (end - start) / (num - 1);

    for (int i = 0; i < num; i++) {
        gsl_vector_set(result.gsl_vector_ptr, i, start + i * step);
    }

    return result;
}

// Function to find the roots of the polynomial using iterative refinement
Matrix roots(const Vector *coeff)
{
    int n = coeff->size - 1; // Degree of the polynomial
    Matrix result;
    result.rows = n;
    result.cols = 2; // Two columns for real and imaginary parts

    // Allocate memory for the gsl_matrix (n rows, 2 columns)
    result.gsl_matrix_ptr = gsl_matrix_alloc(result.rows, result.cols);
    if (result.gsl_matrix_ptr == NULL) 
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Initial guess for roots (uniformly spaced on the unit circle)
    double pi = 3.14159265358979323846;
    double angle;
    for (int i = 0; i < n; i++) 
    {
        angle = 2 * pi * i / n; // Angle for unit circle
        gsl_matrix_set(result.gsl_matrix_ptr, i, 0, cos(angle)); // Real part
        gsl_matrix_set(result.gsl_matrix_ptr, i, 1, sin(angle)); // Imaginary part
    }

    // Iterative method to refine the root estimates
    double real, imag, p, dp, denominator;
    for (int iter = 0; iter < 100; iter++) // Number of iterations
    { 
        for (int i = 0; i < n; i++) 
        {
            real = gsl_matrix_get(result.gsl_matrix_ptr, i, 0);  // Real part
            imag = gsl_matrix_get(result.gsl_matrix_ptr, i, 1);  // Imaginary part

            // Evaluate the polynomial and its derivative at the current estimate
            p  = 0.0; // Polynomial value
            dp = 0.0; // Derivative value

            for (int j = 0; j < coeff->size; j++) 
            {
                p += gsl_vector_get(coeff->gsl_vector_ptr, j) * pow(real, n - j);
                if (j < n) {
                    dp += (n - j) * gsl_vector_get(coeff->gsl_vector_ptr, j) * pow(real, n - j - 1);
                }
            }

            // Update the root estimate
            denominator = p / (dp + 1e-10); // Avoid division by zero
            gsl_matrix_set(result.gsl_matrix_ptr, i, 0, real - denominator); // Update real part
            gsl_matrix_set(result.gsl_matrix_ptr, i, 1, imag - denominator); // Update imaginary part
        }
    }

    return result;
}

Vector polycoefs(const Vector *roots) {
    int n = roots->size; // Number of roots
    Vector coeff;
    coeff.size = n + 1; // Coefficients array size is degree + 1
    coeff.gsl_vector_ptr = gsl_vector_alloc(coeff.size);

    if (coeff.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    gsl_poly_complex_workspace *workspace = gsl_poly_complex_workspace_alloc(n + 1);
    gsl_vector_set(coeff.gsl_vector_ptr, 0, 1.0); // Leading coefficient (x^n)

    for (int i = 0; i < n; i++) {
        double root = gsl_vector_get(roots->gsl_vector_ptr, i);

        for (int j = i; j >= 0; j--) {
            gsl_vector_set(coeff.gsl_vector_ptr, j + 1,
                           gsl_vector_get(coeff.gsl_vector_ptr, j + 1) - root * gsl_vector_get(coeff.gsl_vector_ptr, j));
        }
    }

    gsl_poly_complex_workspace_free(workspace);

    return coeff;
}

Vector conv(const Vector *a, const Vector *b) {
    int n = a->size + b->size - 1; // Size of the result vector
    Vector result;
    result.size = n;
    result.gsl_vector_ptr = gsl_vector_alloc(result.size);

    if (result.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Perform convolution
    for (int i = 0; i < result.size; i++) {
        gsl_vector_set(result.gsl_vector_ptr, i, 0.0);
    }

    for (int i = 0; i < a->size; i++) {
        for (int j = 0; j < b->size; j++) {
            gsl_vector_set(result.gsl_vector_ptr, i + j,
                           gsl_vector_get(result.gsl_vector_ptr, i + j) +
                           gsl_vector_get(a->gsl_vector_ptr, i) * gsl_vector_get(b->gsl_vector_ptr, j));
        }
    }

    return result;
}
