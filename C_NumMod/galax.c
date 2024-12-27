/******************************************************/
//  GALAX - GSL Advanced Analytics eXtension
/******************************************************/

//#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>

typedef struct {
    gsl_matrix *gsl_matrix_ptr;
    unsigned int rows;
    unsigned int cols;
} Matrix;

typedef struct {
    gsl_vector *gsl_vector_ptr;
    unsigned int size;
} Vector;

// Define a struct to hold ODE parameters
typedef struct {
    gsl_odeiv2_system system;
    gsl_odeiv2_driver *driver;
} ode_solver;

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

    gsl_matrix_memcpy(inverse.gsl_matrix_ptr, mat->gsl_matrix_ptr);

    gsl_permutation *p = gsl_permutation_alloc(n);
    int signum;

    // Perform LU decomposition
    gsl_linalg_LU_decomp(inverse.gsl_matrix_ptr, p, &signum);

    // Compute the inverse directly
    gsl_linalg_LU_invert(inverse.gsl_matrix_ptr, p, inverse.gsl_matrix_ptr);

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

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            gsl_matrix_set(mat.gsl_matrix_ptr, i, j, gsl_rng_uniform(rng));
        }
    }

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

// Function to compute the roots of a polynomial using the Durand-Kerner method
Matrix roots(const Vector *coeff) {
    unsigned int n = coeff->size - 1; // Degree of the polynomial
    if (n < 1) {
        fprintf(stderr, "Polynomial degree must be at least 1.\n");
        exit(EXIT_FAILURE);
    }

    // Initialize the output matrix
    Matrix result;
    result.rows = n;
    result.cols = 2; // Real and imaginary parts
    result.gsl_matrix_ptr = gsl_matrix_alloc(n, 2);
    if (result.gsl_matrix_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for result matrix.\n");
        exit(EXIT_FAILURE);
    }

    // Initial guess: roots uniformly spaced on the unit circle
    double pi = 3.14159265358979323846;
    for (unsigned int i = 0; i < n; i++) {
        double angle = 2 * pi * i / n; // Angle for the unit circle
        gsl_matrix_set(result.gsl_matrix_ptr, i, 0, cos(angle)); // Real part
        gsl_matrix_set(result.gsl_matrix_ptr, i, 1, sin(angle)); // Imaginary part
    }

    // Iterative refinement using the Durand-Kerner method
    int max_iter = 1000;
    double tolerance = 1e-14;
    double real_i, imag_i, p_real, p_imag, denom_real, denom_imag;
    double real_j, imag_j, diff_real, diff_imag, mag_sq, new_real, new_imag;
    double coeff_real, temp_real, temp_imag, delta_real, delta_imag;
    for (unsigned int iter = 0; iter < max_iter; iter++) {
        int converged = 1; // Flag to check for convergence

        for (unsigned int i = 0; i < n; i++) {
            real_i = gsl_matrix_get(result.gsl_matrix_ptr, i, 0);
            imag_i = gsl_matrix_get(result.gsl_matrix_ptr, i, 1);

            // Compute the value of the polynomial and its correction term
            p_real = gsl_vector_get(coeff->gsl_vector_ptr, 0); // Leading coefficient
            p_imag = 0.0; // No imaginary part initially
            denom_real = 1.0, denom_imag = 0.0; // Product term denominator

            for (unsigned int j = 0; j < n; j++) {
                if (j == i) continue; // Skip the current root
                real_j = gsl_matrix_get(result.gsl_matrix_ptr, j, 0);
                imag_j = gsl_matrix_get(result.gsl_matrix_ptr, j, 1);

                // Compute (z_i - z_j) and its magnitude squared
                diff_real = real_i - real_j;
                diff_imag = imag_i - imag_j;
                mag_sq = diff_real * diff_real + diff_imag * diff_imag;

                // Update the denominator as a complex product
                new_real = denom_real * diff_real - denom_imag * diff_imag;
                new_imag = denom_real * diff_imag + denom_imag * diff_real;
                denom_real = new_real;
                denom_imag = new_imag;
            }

            // Evaluate the polynomial at z_i
            for (unsigned int k = 1; k <= n; k++) {
                coeff_real = gsl_vector_get(coeff->gsl_vector_ptr, k);
                temp_real = p_real * real_i - p_imag * imag_i + coeff_real;
                temp_imag = p_real * imag_i + p_imag * real_i;
                p_real = temp_real;
                p_imag = temp_imag;
            }

            // Compute the correction term
            mag_sq = denom_real * denom_real + denom_imag * denom_imag;
            delta_real = (p_real * denom_real + p_imag * denom_imag) / mag_sq;
            delta_imag = (p_imag * denom_real - p_real * denom_imag) / mag_sq;

            // Update the root estimate
            new_real = real_i - delta_real;
            new_imag = imag_i - delta_imag;

            gsl_matrix_set(result.gsl_matrix_ptr, i, 0, new_real);
            gsl_matrix_set(result.gsl_matrix_ptr, i, 1, new_imag);

            // Check for convergence
            if (fabs(delta_real) > tolerance || fabs(delta_imag) > tolerance) {
                converged = 0;
            }
        }

        if (converged) break; // Stop if all roots have converged
    }
    return result;
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

// Function to convert roots to polynomial coefficients
Vector polycoefs(const Vector *roots) {
    int n = roots->size; // Number of roots
    Vector coeff;
    coeff.size = n + 1; // Coefficients array size is degree + 1
    coeff.gsl_vector_ptr = gsl_vector_alloc(coeff.size);

    if (coeff.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Initialize all coefficients to zero
    for (int i = 0; i < coeff.size; i++) {
        gsl_vector_set(coeff.gsl_vector_ptr, i, 0.0);
    }

    gsl_poly_complex_workspace *workspace = gsl_poly_complex_workspace_alloc(n + 1);
    gsl_vector_set(coeff.gsl_vector_ptr, 0, 1.0); // Leading coefficient (x^n)

    for (int i = 0; i < n; i++) {
        double root = gsl_vector_get(roots->gsl_vector_ptr, i);

        for (int j = i; j >= 0; j--) {
            double new_val = gsl_vector_get(coeff.gsl_vector_ptr, j + 1) - root * gsl_vector_get(coeff.gsl_vector_ptr, j);
            gsl_vector_set(coeff.gsl_vector_ptr, j + 1, new_val);

            // Debug print to verify each step
            //printf("Debug: coeff[%d] = %.4f after setting with root %.4f\n", j + 1, new_val, root);
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

// Function to interpolate using specified method ("linear" or "spline")
Vector interp1(const Vector *s, const Vector *t, const Vector *ss, const char *method)
{
    // Check for dimension compatibility
    if (s == NULL || t == NULL || ss == NULL || s->size != t->size) {
        fprintf(stderr, "Error: Input vectors are incompatible.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the result vector
    Vector result;
    result.size = ss->size;
    result.gsl_vector_ptr = gsl_vector_alloc(result.size);
    if (result.gsl_vector_ptr == NULL) {
        fprintf(stderr, "Memory allocation failed for result vector.\n");
        exit(EXIT_FAILURE);
    }

    // Determine the interpolation method
    const gsl_interp_type *interp_type;
    if (strcmp(method, "linear") == 0) {
        interp_type = gsl_interp_linear;
    } else if (strcmp(method, "spline") == 0) {
        interp_type = gsl_interp_cspline;
    } else {
        fprintf(stderr, "Error: Unsupported interpolation method '%s'.\n", method);
        gsl_vector_free(result.gsl_vector_ptr);
        exit(EXIT_FAILURE);
    }

    // Prepare the data for GSL (convert double to double for GSL compatibility)
    gsl_vector *s_double = gsl_vector_alloc(s->size);
    gsl_vector *t_double = gsl_vector_alloc(t->size);
    gsl_vector *ss_double = gsl_vector_alloc(ss->size);
    if (!s_double || !t_double || !ss_double) {
        fprintf(stderr, "Error: Memory allocation failed for temporary arrays.\n");
        gsl_vector_free(result.gsl_vector_ptr);
        gsl_vector_free(s_double);
        gsl_vector_free(t_double);
        gsl_vector_free(ss_double);
        exit(EXIT_FAILURE);
    }

    // Copy data into the gsl_vectors
    for (unsigned int i = 0; i < s->size; i++) {
        gsl_vector_set(s_double, i, s->gsl_vector_ptr->data[i]);
        gsl_vector_set(t_double, i, t->gsl_vector_ptr->data[i]);
    }
    for (unsigned int i = 0; i < ss->size; i++) {
        gsl_vector_set(ss_double, i, ss->gsl_vector_ptr->data[i]);
    }

    // Create a GSL interpolation object and accelerator
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(interp_type, s->size);

    // Initialize the interpolation object
    if (gsl_spline_init(spline, s_double->data, t_double->data, s->size) != GSL_SUCCESS) {
        fprintf(stderr, "Error: Failed to initialize interpolation.\n");
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
        gsl_vector_free(result.gsl_vector_ptr);
        gsl_vector_free(s_double);
        gsl_vector_free(t_double);
        gsl_vector_free(ss_double);
        exit(EXIT_FAILURE);
    }

    // Perform interpolation for each value in ss
    for (unsigned int i = 0; i < ss->size; i++) {
        if (gsl_vector_get(ss_double, i) < gsl_vector_get(s_double, 0) || 
            gsl_vector_get(ss_double, i) > gsl_vector_get(s_double, s->size - 1)) {
            fprintf(stderr, "Error: Interpolation point ss[%d] = %f is out of bounds.\n", i, gsl_vector_get(ss_double, i));
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            gsl_vector_free(result.gsl_vector_ptr);
            gsl_vector_free(s_double);
            gsl_vector_free(t_double);
            gsl_vector_free(ss_double);
            exit(EXIT_FAILURE);
        }
        gsl_vector_set(result.gsl_vector_ptr, i, gsl_spline_eval(spline, gsl_vector_get(ss_double, i), acc));
    }

    // Free GSL resources
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    gsl_vector_free(s_double);
    gsl_vector_free(t_double);
    gsl_vector_free(ss_double);

    return result; // Return the interpolated vector
}


////////////////////////////////////////////////////////////////////////////////////////////////////////

double integral(double (*func)(double, void *), double xmin, double xmax) {
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

    double result, error;
    gsl_function F;
    F.function = func;
    F.params = NULL;
    double epsabs = 1.49e-08;
    double epsrel = 1.49e-08;

    // Perform the integration
    gsl_integration_qags(&F, xmin, xmax, epsabs, epsrel, 1000, workspace, &result, &error);

    gsl_integration_workspace_free(workspace);
    return result;
}

// Function to create an ODE solver
ode_solver *create_ode_solver(int (*func)(double, const double[], double[], void *), void *params, size_t dim, double hstart, double epsabs, double epsrel) {
    ode_solver *solver = (ode_solver *)malloc(sizeof(ode_solver));
    solver->system.function = func;
    solver->system.jacobian = NULL;
    solver->system.dimension = dim;
    solver->system.params = params;
    solver->driver = gsl_odeiv2_driver_alloc_y_new(&solver->system, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);
    return solver;
}

// Function to solve ODEs with a given time step
void solve_ode(ode_solver *solver, double *t, double t1, double y[], double dt) {
    while (*t < t1) {
        int status = gsl_odeiv2_driver_apply(solver->driver, t, *t + dt, y);
        if (status != GSL_SUCCESS) {
            printf("error: return value=%d\n", status);
            break;
        }

        printf("%.5e", *t);  // Print time
        for (size_t i = 0; i < solver->system.dimension; i++) {
            printf(" %.5e", y[i]);  // Print each element of y
        }
        printf("\n");
    }
}

// Function to free the ODE solver
void free_ode_solver(ode_solver *solver) {
    gsl_odeiv2_driver_free(solver->driver);
    free(solver);
}

// Function to solve ODE with default parameters and specified dimension
void ode45(int (*func)(double, const double[], double[], void *), void *params, size_t dimension, double *t, double t1, double y[], double dt) {
    // Default parameters for the solver
    double initial_step_size = 1e-6;
    double absolute_error = 1e-6;
    double relative_error = 0.0;

    ode_solver *solver = create_ode_solver(func, params, dimension, initial_step_size, absolute_error, relative_error);
    solve_ode(solver, t, t1, y, dt);
    free_ode_solver(solver);
}






