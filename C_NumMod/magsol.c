/*******************************************/
/* MAGSOL - Matrix Algebra and GSL SOLvers */
/*******************************************/

//      gcc magsol.c -o magsol -lgsl -lm

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
    gsl_odeiv2_system system;
    gsl_odeiv2_driver *driver;
} ode_solver;

//////////////////////////////////////////////////////////////////////////////////////

void          freeMatrix(gsl_matrix *matrix);
gsl_matrix*   zerosMatrix(int rows, int cols);
gsl_matrix*   createMatrix(int rows, int cols, double *values);
void          printMatrix(const gsl_matrix *mat);
gsl_matrix*   transposeMatrix(const gsl_matrix *mat);
gsl_matrix*   inverseMatrix(const gsl_matrix *mat);
gsl_matrix*   matmul(const gsl_matrix *A, const gsl_matrix *B);
gsl_matrix*   matadd(const gsl_matrix *A, const gsl_matrix *B);
gsl_matrix*   matscal(const gsl_matrix *A, double scalar);
double        det(const gsl_matrix *mat);
gsl_matrix*   matelem(const gsl_matrix *A, const gsl_matrix *B);
gsl_matrix*   matrand(int rows, int cols);
gsl_matrix*   compan(const gsl_vector *coefficients);
void          eig(const gsl_matrix *matrix);
double        normMatrix(const gsl_matrix *mat, double p);
gsl_matrix*   linsolve_overdet(const gsl_matrix *A, const gsl_matrix *F);
gsl_matrix*   roots(const gsl_vector *coeff);
void          freeVector(gsl_vector *vector);
gsl_vector*   zerosVector(int size);
gsl_vector*   createVector(int size, double *values);
void          printVector(const gsl_vector *vector);
gsl_vector*   matvec(const gsl_matrix *A, const gsl_vector *x);
gsl_vector*   linsolve(const gsl_matrix *A, const gsl_vector *b);
gsl_vector*   vecadd(const gsl_vector *A, const gsl_vector *B);
gsl_vector*   vecscal(const gsl_vector *A, double scalar);
gsl_vector*   vecrand(int dim);
gsl_vector*   vecelem(const gsl_vector *A, const gsl_vector *B);
double        norm(const gsl_vector *vec, double p);
gsl_vector*   createArray(double start, double end, double step);
gsl_vector*   polyfitweighted(const gsl_vector *x, const gsl_vector *y, const gsl_vector *w, unsigned int n);
gsl_vector*   linspace(double start, double end, int num);
gsl_vector*   polycoefs(const gsl_vector *roots);
gsl_vector*   conv(const gsl_vector *a, const gsl_vector *b);
gsl_vector*   interp1(const gsl_vector *s, const gsl_vector *t, const gsl_vector *ss, const char *method);
double        integral(double (*func)(double, void *), double xmin, double xmax);
ode_solver*   create_ode_solver(int (*func)(double, const double[], double[], void *), void *params, size_t dim, double hstart, double epsabs, double epsrel);
void          solve_ode(ode_solver *solver, double *t, double t1, double y[], double dt);
void          free_ode_solver(ode_solver *solver);
void          ode45(int (*func)(double, const double[], double[], void *), void *params, size_t dimension, double *t, double t1, double y[], double dt);

// integral() test function
double fun(double x, void *params){
    (void)params;
    return (x*x)/(1.0 + sin(x) + cos(x));
}

// ODE Solver test function
int ode_func(double t, const double y[], double f[], void *params) {
    (void)(t);  // avoid unused parameter warning
    double mu = *(double *)params;
    f[0] = y[1];
    f[1] = -mu * (1 - y[0] * y[0]) * y[1] - y[0];
    return GSL_SUCCESS;
}

int main() {
    // Matrix and vector creation for testing
    double matrix_data_A[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double matrix_data_B[] = {9, 8, 7, 6, 5, 4, 3, 2, 1};
    double vector_data[] = {1, 2, 3};
    
    gsl_matrix *A = createMatrix(3, 3, matrix_data_A);
    gsl_matrix *B = createMatrix(3, 3, matrix_data_B);
    gsl_vector *v = createVector(3, vector_data);

    // Test printing matrix A
    printf("Matrix A:\n");
    printMatrix(A);

    // Test printing matrix B
    printf("\nMatrix B:\n");
    printMatrix(B);

    // Test matrix transpose
    gsl_matrix *At = transposeMatrix(A);
    printf("\nTranspose of A:\n");
    printMatrix(At);

    // Test matrix inverse (for a simple example, use a 2x2 invertible matrix)
    double matrix_data_C[] = {4, 7, 2, 6};
    gsl_matrix *C = createMatrix(2, 2, matrix_data_C);
    printf("\nMatrix C:\n");
    printMatrix(C);

    gsl_matrix *C_inv = inverseMatrix(C);
    printf("\nInverse of C:\n");
    printMatrix(C_inv);
    gsl_matrix_free(C);
    gsl_matrix_free(C_inv);

    // Test matrix multiplication
    gsl_matrix *AB = matmul(A, B);
    printf("\nA * B:\n");
    printMatrix(AB);

    // Test matrix addition
    gsl_matrix *A_plus_B = matadd(A, B);
    printf("\nA + B:\n");
    printMatrix(A_plus_B);

    // Test matrix scalar multiplication
    gsl_matrix *A_scaled = matscal(A, 2.0);
    printf("\n2 * A:\n");
    printMatrix(A_scaled);

    // Test matrix determinant (for the same 2x2 matrix)
    C = createMatrix(2, 2, matrix_data_C);
    double det_C = det(C);
    printf("\nDeterminant of C: %.4f\n", det_C);
    gsl_matrix_free(C);

    // Test element-wise matrix multiplication
    gsl_matrix *A_elem_B = matelem(A, B);
    printf("\nElement-wise multiplication of A and B:\n");
    printMatrix(A_elem_B);

    // Test random matrix generation
    gsl_matrix *R = matrand(3, 3);
    printf("\nRandom matrix R:\n");
    printMatrix(R);

    // Free allocated memory
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(At);
    gsl_matrix_free(AB);
    gsl_matrix_free(A_plus_B);
    gsl_matrix_free(A_scaled);
    gsl_matrix_free(A_elem_B);
    gsl_matrix_free(R);
    gsl_vector_free(v);

    // Test the companion matrix function
    double coeff_data[] = {1, -6, 11, -6}; // Polynomial: x^3 - 6x^2 + 11x - 6
    gsl_vector *coefficients = createVector(4, coeff_data);
    gsl_matrix *compan_matrix = compan(coefficients);
    printf("\nCompanion matrix:\n");
    printMatrix(compan_matrix);
    gsl_matrix_free(compan_matrix);
    gsl_vector_free(coefficients);

    // Test the eigenvalues function
    double eig_data[] = {4, -1, -1, 4};
    gsl_matrix *eig_matrix = createMatrix(2, 2, eig_data);
    printf("\nEigenvalues of matrix:\n");
    eig(eig_matrix);
    gsl_matrix_free(eig_matrix);

    // Test the norm function
    double norm_data[] = {  1, 2, 3, 
                            4, 5, 6};
    gsl_matrix *norm_matrix = createMatrix(2, 3, norm_data);
    double norm_val = normMatrix(norm_matrix, 2);
    printf("\nFrobenius norm of matrix: %.4f\n", norm_val);
    gsl_matrix_free(norm_matrix);

    // Test the overdetermined system solver
    double A_data[] = {1, 1, 
                       1, 1, 
                       2, 3};
    double F_data[] = {1, 
                       2, 
                       3};
    gsl_matrix *A_matrix = createMatrix(3, 2, A_data);
    gsl_matrix *F_matrix = createMatrix(3, 1, F_data);
    gsl_matrix *sol_matrix = linsolve_overdet(A_matrix, F_matrix);
    printf("\nSolution to overdetermined system:\n");
    printMatrix(sol_matrix);
    gsl_matrix_free(A_matrix);
    gsl_matrix_free(F_matrix);
    gsl_matrix_free(sol_matrix);

    // Test the roots function using Durand-Kerner method
    double roots_data[] = {1, -6, 11, -6}; // Polynomial: x^3 - 6x^2 + 11x - 6
    gsl_vector *roots_coeff = createVector(4, roots_data);
    gsl_matrix *roots_matrix = roots(roots_coeff);
    printf("\nRoots of polynomial:\n");
    printMatrix(roots_matrix);
    gsl_matrix_free(roots_matrix);
    gsl_vector_free(roots_coeff);


    // Test zerosVector function
    gsl_vector *zero_vec = zerosVector(5);
    printf("\nZero Vector:\n");
    printVector(zero_vec);

    // Test createVector function
    double vec_data[] = {1.1, 2.2, 3.3, 4.4, 5.5};
    gsl_vector *created_vec = createVector(5, vec_data);
    printf("\nCreated Vector:\n");
    printVector(created_vec);

    // Test vector addition
    gsl_vector *added_vec = vecadd(created_vec, created_vec);
    printf("\nVector Addition (v + v):\n");
    printVector(added_vec);

    // Test scalar multiplication of vector
    gsl_vector *scaled_vec = vecscal(created_vec, 2.0);
    printf("\nScalar Multiplication of Vector (2 * v):\n");
    printVector(scaled_vec);

    // Test random vector generation
    gsl_vector *random_vec = vecrand(5);
    printf("\nRandom Vector:\n");
    printVector(random_vec);

    // Test element-wise vector multiplication
    gsl_vector *elem_mult_vec = vecelem(created_vec, created_vec);
    printf("\nElement-wise Multiplication of Vector (v .* v):\n");
    printVector(elem_mult_vec);

    // Test vector norm calculation
    double vec_norm = norm(created_vec, 2);
    printf("\n2-Norm of the Vector: %.4f\n", vec_norm);

    // Test matrix-vector multiplication
    double Adata[] = {1, 2, 3, 
                      4, 5, 6};
    gsl_matrix *A_mat = createMatrix(2, 3, Adata);
    double cdata[] = {5, 3, 1};
    gsl_vector *c_vec = createVector(3, cdata);
    gsl_vector *mat_vec_mult = matvec(A_mat, c_vec);
    printf("\nMatrix-Vector Multiplication (A * v):\n");
    printVector(mat_vec_mult);
    gsl_matrix_free(A_mat);
    gsl_vector_free(c_vec);

    // Test solving linear system
    double A_Data[] = {3, 2, -1, 2, -2, 4, -1, 0.5, -1};
    double b_Data[] = {1, -2, 0};
    gsl_matrix *A_sys = createMatrix(3, 3, A_Data);
    gsl_vector *b_sys = createVector(3, b_Data);
    gsl_vector *x = linsolve(A_sys, b_sys);
    printf("\nSolution to Linear System (Ax = b):\n");
    printVector(x);
    gsl_matrix_free(A_sys);
    gsl_vector_free(b_sys);
    gsl_vector_free(x);

    // Test creating array with start, end, and step
    gsl_vector *array = createArray(0, 10, 2);
    printf("\nCreated Array (start=0, end=10, step=2):\n");
    printVector(array);
    gsl_vector_free(array);

    // Test linspace function
    gsl_vector *linspace_vec = linspace(0, 1, 5);
    printf("\nLinspace (start=0, end=1, num=5):\n");
    printVector(linspace_vec);
    gsl_vector_free(linspace_vec);

    // Test converting roots to polynomial coefficients
    double root_data[] = {-1, 1, 2}; // Roots of the polynomial
    gsl_vector *roots = createVector(3, root_data);
    gsl_vector *coefs = polycoefs(roots);
    printf("\nPolynomial Coefficients from Roots:\n");
    printVector(coefs);
    gsl_vector_free(roots);
    gsl_vector_free(coefs);

    // Test polynomial fit with weights
    double x_data[] = {1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
    double y_data[] = {9.08, 10.43, 11.9, 13.48, 15.19, 17.03, 19.01, 21.13, 23.39};
    double w_data[] = {1.0, 1.0, 2.0, 5.0, 1.0, 4.0, 2.0, 2.0, 1.0};
    gsl_vector *x_vec = createVector(9, x_data);
    gsl_vector *y_vec = createVector(9, y_data);
    gsl_vector *w_vec = createVector(9, w_data);
    gsl_vector *poly_fit = polyfitweighted(x_vec, y_vec, w_vec, 3);
    printf("\nPolynomial Fit with Weights:\n");
    printVector(poly_fit);
    gsl_vector_free(x_vec);
    gsl_vector_free(y_vec);
    gsl_vector_free(w_vec);
    gsl_vector_free(poly_fit);

    // Test convolution
    double a_data[] = {1, 2, 3};
    double b_data[] = {0, 1, 0.5};
    gsl_vector *a_vec = createVector(3, a_data);
    gsl_vector *b_vec = createVector(3, b_data);
    gsl_vector *conv_result = conv(a_vec, b_vec);
    printf("\nConvolution of Vectors a and b:\n");
    printVector(conv_result);
    gsl_vector_free(a_vec);
    gsl_vector_free(b_vec);
    gsl_vector_free(conv_result);

    // Free allocated memory
    gsl_vector_free(zero_vec);
    gsl_vector_free(created_vec);
    gsl_vector_free(added_vec);
    gsl_vector_free(scaled_vec);
    gsl_vector_free(random_vec);
    gsl_vector_free(elem_mult_vec);
    gsl_vector_free(mat_vec_mult);

    // Test interpolation function
    double s_data[] = {0, 1, 2, 3, 4};
    double t_data[] = {0, 1, 4, 9, 16};
    double ss_data[] = {0.5, 1.5, 2.5, 3.5};
    gsl_vector *s = createVector(5, s_data);
    gsl_vector *t = createVector(5, t_data);
    gsl_vector *ss = createArray(0, 4, 0.4);
    printf("\nVector ss length: %zu\n", ss->size);
    printf("ss = "); printVector(ss);

    gsl_vector *interp_linear = interp1(s, t, ss, "linear");
    printf("\nLinear Interpolation:\n");
    printVector(interp_linear);

    gsl_vector *interp_spline = interp1(s, t, ss, "spline");
    printf("\nSpline Interpolation:\n");
    printVector(interp_spline);

    gsl_vector_free(s);
    gsl_vector_free(t);
    gsl_vector_free(ss);
    gsl_vector_free(interp_linear);
    gsl_vector_free(interp_spline);

    // Test integral function
    double integral_result = integral(fun, 0, M_PI/2.0);
    printf("\nIntegral of fun(x) from 0 to pi/2: %.4f\n", integral_result);

    // Test ODE solver
    double params = 10.0;
    size_t dimension = 2;
    double tt = 0.0, t1 = 1.0;
    double y[2] = {1.0, 0.0};
    double dt = 1e-1;

    printf("\nODE Solver Results:\n");
    ode45(ode_func, &params, dimension, &tt, t1, y, dt);

    // Free allocated memory for remaining variables
    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////

void freeMatrix(gsl_matrix *matrix) {
    if (matrix != NULL) {
        gsl_matrix_free(matrix);
    }
}

gsl_matrix* zerosMatrix(int rows, int cols) {
    return gsl_matrix_calloc(rows, cols);
}

gsl_matrix* createMatrix(int rows, int cols, double *values) {
    gsl_matrix *matrix = gsl_matrix_alloc(rows, cols);

    if (matrix == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Populate the matrix with the provided values
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            gsl_matrix_set(matrix, i, j, values[i * cols + j]);
        }
    }

    return matrix;
}

void printMatrix(const gsl_matrix *mat) {
    int rows = mat->size1;
    int cols = mat->size2;
    int max_print_size = 3; // Print 3 rows and 3 columns from corners

    printf("[\n");

    if (rows <= 10 && cols <= 10) {
        // Print the entire matrix if both dimensions are 10 or less
        for (int i = 0; i < rows; i++) {
            printf("  [");
            for (int j = 0; j < cols; j++) {
                printf("%7.4f", gsl_matrix_get(mat, i, j));
                if (j < cols - 1) printf(", ");
            }
            printf("]");
            if (i < rows - 1) printf(",\n");
        }
    } else {
        // Print the top 3 rows
        for (int i = 0; i < max_print_size && i < rows; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size && j < cols; j++) {
                printf("%7.4f", gsl_matrix_get(mat, i, j));
                if (j < max_print_size - 1 && j < cols - 1) printf(", ");
            }
            // Print ellipsis for omitted columns
            if (cols > max_print_size) {
                printf(", ... ");
            }
            // Last elements of the top 3 rows
            for (int j = cols - 3; j < cols; j++) {
                if (j >= 0) { // Ensure we don't go out of bounds
                    printf("%7.4f", gsl_matrix_get(mat, i, j));
                    if (j < cols - 1) printf(", ");
                }
            }
            printf("]");
            if (i < max_print_size - 1 && i < rows - 1) printf(",\n");
        }

        // Print ellipsis for omitted rows
        if (rows > max_print_size) {
            printf(",\n\t\t\t      ...\n");
        }

        // Print the bottom 3 rows
        for (int i = rows - max_print_size; i < rows; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size && j < cols; j++) {
                printf("%7.4f", gsl_matrix_get(mat, i, j));
                if (j < max_print_size - 1 && j < cols - 1) printf(", ");
            }
            // Print ellipsis for omitted columns
            if (cols > max_print_size) {
                printf(", ... ");
            }
            // Last elements of the bottom 3 rows
            for (int j = cols - 3; j < cols; j++) {
                if (j >= 0) { // Ensure we don't go out of bounds
                    printf("%7.4f", gsl_matrix_get(mat, i, j));
                    if (j < cols - 1) printf(", ");
                }
            }
            printf("]");
            if (i < rows - 1) printf(",\n");
        }
    }

    printf("\n]\n");
}

gsl_matrix* transposeMatrix(const gsl_matrix *mat) {
    size_t rows = mat->size2;
    size_t cols = mat->size1;

    gsl_matrix *transposed = gsl_matrix_alloc(rows, cols);

    if (transposed == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < mat->size1; i++) {
        for (size_t j = 0; j < mat->size2; j++) {
            gsl_matrix_set(transposed, j, i, gsl_matrix_get(mat, i, j));
        }
    }

    return transposed;
}

gsl_matrix* inverseMatrix(const gsl_matrix *mat) {
    if (mat->size1 != mat->size2) {
        fprintf(stderr, "Matrix must be square to compute the inverse.\n");
        exit(EXIT_FAILURE);
    }

    int n = mat->size1;
    gsl_matrix *inverse = gsl_matrix_alloc(n, n);

    if (inverse == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    gsl_matrix_memcpy(inverse, mat);

    gsl_permutation *p = gsl_permutation_alloc(n);
    int signum;

    // Perform LU decomposition
    gsl_linalg_LU_decomp(inverse, p, &signum);

    // Compute the inverse directly
    gsl_linalg_LU_invert(inverse, p, inverse);

    gsl_permutation_free(p);

    return inverse;
}

gsl_matrix* matmul(const gsl_matrix *A, const gsl_matrix *B) {
    if (A->size2 != B->size1) {
        fprintf(stderr, "Matrix dimensions are incompatible for multiplication.\n");
        exit(EXIT_FAILURE);
    }

    gsl_matrix *C = gsl_matrix_alloc(A->size1, B->size2);

    // Perform matrix multiplication
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, C);

    return C;
}

gsl_matrix* matadd(const gsl_matrix *A, const gsl_matrix *B) {
    if (A->size1 != B->size1 || A->size2 != B->size2) {
        fprintf(stderr, "Matrix dimensions are incompatible for addition.\n");
        exit(EXIT_FAILURE);
    }

    gsl_matrix *C = gsl_matrix_alloc(A->size1, A->size2);

    // Copy B into C
    gsl_matrix_memcpy(C, B);

    // Perform matrix addition
    gsl_matrix_add(C, A);

    return C;
}

gsl_matrix* matscal(const gsl_matrix *A, double scalar) {
    gsl_matrix *C = gsl_matrix_alloc(A->size1, A->size2);

    // Copy A into C
    gsl_matrix_memcpy(C, A);

    // Scale the matrix by a constant
    gsl_matrix_scale(C, scalar);

    return C;
}

double det(const gsl_matrix *mat) {
    // Check if the matrix is square
    if (mat->size1 != mat->size2) {
        fprintf(stderr, "Matrix must be square to compute the determinant.\n");
        exit(EXIT_FAILURE);
    }

    int n = mat->size1;
    gsl_matrix *tmp = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(tmp, mat);

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

gsl_matrix* matelem(const gsl_matrix *A, const gsl_matrix *B) {
    // Check if both matrices have the same dimensions
    if (A->size1 != B->size1 || A->size2 != B->size2) {
        fprintf(stderr, "Matrices must have the same dimensions for element-wise multiplication.\n");
        exit(EXIT_FAILURE);
    }

    gsl_matrix *C = gsl_matrix_alloc(A->size1, A->size2);

    // Perform element-wise multiplication
    for (size_t i = 0; i < C->size1; i++) {
        for (size_t j = 0; j < C->size2; j++) {
            double value = gsl_matrix_get(A, i, j) * gsl_matrix_get(B, i, j);
            gsl_matrix_set(C, i, j, value);
        }
    }

    return C;
}

gsl_matrix* matrand(int rows, int cols) {
    gsl_matrix *mat = gsl_matrix_alloc(rows, cols);

    if (mat == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            gsl_matrix_set(mat, i, j, gsl_rng_uniform(rng));
        }
    }

    gsl_rng_free(rng);

    return mat;
}


gsl_matrix* compan(const gsl_vector *coefficients) {
    int size = coefficients->size;

    // Input validation
    if (size < 2) {
        fprintf(stderr, "The input vector must have at least 2 elements (degree 1 polynomial).\n");
        return NULL; // Return NULL for error case
    }

    if (gsl_vector_get(coefficients, 0) == 0) {
        fprintf(stderr, "The first coefficient in the vector must not be zero.\n");
        return NULL; // Return NULL for error case
    }

    // Create the companion matrix
    int n = size - 1; // Degree of the polynomial
    gsl_matrix *C = gsl_matrix_calloc(n, n);

    if (C == NULL) {
        fprintf(stderr, "Memory allocation failed for companion matrix.\n");
        return NULL; // Return NULL for error case
    }

    // Fill the first row
    double leading_coefficient = gsl_vector_get(coefficients, 0);
    for (int j = 0; j < n; j++) {
        gsl_matrix_set(C, 0, j, -gsl_vector_get(coefficients, j + 1) / leading_coefficient);
    }

    // Fill the sub-diagonal with ones
    for (int i = 1; i < n; i++) {
        gsl_matrix_set(C, i, i - 1, 1.0);
    }

    return C;
}

void eig(const gsl_matrix *matrix) {
    if (matrix->size1 != matrix->size2) {
        fprintf(stderr, "The matrix must be square to compute eigenvalues.\n");
        exit(EXIT_FAILURE);
    }

    int n = matrix->size1;
    gsl_matrix *result = gsl_matrix_alloc(n, 2); // To store real and imaginary parts

    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(n);
    gsl_vector_complex *eval = gsl_vector_complex_alloc(n);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(n, n);
    gsl_matrix *temp = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(temp, matrix);

    gsl_eigen_nonsymmv(temp, eval, evec, w);

    gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        gsl_complex z = gsl_vector_complex_get(eval, i);
        double real = GSL_REAL(z);
        double imag = GSL_IMAG(z);
        gsl_matrix_set(result, i, 0, real);
        gsl_matrix_set(result, i, 1, imag);
        printf("%7.4f + %7.4fi\n", real, imag);
    }

    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_matrix_free(temp);
    gsl_eigen_nonsymmv_free(w);
    gsl_matrix_free(result);
}

double normMatrix(const gsl_matrix *mat, double p) {
    double norm_value = 0.0;
    size_t rows = mat->size1;
    size_t cols = mat->size2;

    if (p == 1) {
        // 1-norm (maximum column sum)
        for (size_t j = 0; j < cols; j++) {
            double col_sum = 0.0;
            for (size_t i = 0; i < rows; i++) {
                col_sum += fabs(gsl_matrix_get(mat, i, j));
            }
            if (col_sum > norm_value) {
                norm_value = col_sum;
            }
        }
    } else if (p == 2) {
        // Frobenius norm (sum of squares of all elements, then square root)
        double sum = 0.0;
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                double value = gsl_matrix_get(mat, i, j);
                sum += value * value;
            }
        }
        norm_value = sqrt(sum);
    } else if (p == INFINITY) {
        // Infinity-norm (maximum row sum)
        for (size_t i = 0; i < rows; i++) {
            double row_sum = 0.0;
            for (size_t j = 0; j < cols; j++) {
                row_sum += fabs(gsl_matrix_get(mat, i, j));
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

gsl_matrix* linsolve_overdet(const gsl_matrix *A, const gsl_matrix *F) {
    // Ensure the dimensions are compatible
    if (A->size1 != F->size1) {
        fprintf(stderr, "Matrix dimensions are incompatible for solving.\n");
        exit(EXIT_FAILURE);
    }

    size_t n = A->size1;  // Number of equations
    size_t m = A->size2;  // Number of variables
    size_t nrhs = F->size2;  // Number of right-hand sides

    gsl_matrix *x = gsl_matrix_alloc(m, nrhs);

    gsl_matrix *gslA = gsl_matrix_alloc(n, m);
    gsl_vector *gslB = gsl_vector_alloc(n); // Vector instead of matrix for y (right-hand side)
    gsl_matrix *cov = gsl_matrix_alloc(m, m);
    gsl_vector *c = gsl_vector_alloc(m);
    double chisq;

    gsl_matrix_memcpy(gslA, A);

    // We need to solve for each column of F separately
    for (size_t j = 0; j < nrhs; j++) {
        for (size_t i = 0; i < n; i++) {
            gsl_vector_set(gslB, i, gsl_matrix_get(F, i, j));
        }

        gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, m);

        // Solve the system using GSL's least squares solver
        gsl_multifit_linear(gslA, gslB, c, cov, &chisq, work);

        // Copy the results to the output matrix
        for (size_t i = 0; i < m; i++) {
            gsl_matrix_set(x, i, j, gsl_vector_get(c, i));
        }

        gsl_multifit_linear_free(work);
    }

    gsl_matrix_free(gslA);
    gsl_vector_free(gslB);
    gsl_matrix_free(cov);
    gsl_vector_free(c);

    return x;
}

gsl_matrix* roots(const gsl_vector *coeff) {
    unsigned int n = coeff->size - 1; // Degree of the polynomial
    if (n < 1) {
        fprintf(stderr, "Polynomial degree must be at least 1.\n");
        exit(EXIT_FAILURE);
    }

    // Initialize the output matrix
    gsl_matrix *result = gsl_matrix_alloc(n, 2); // Real and imaginary parts
    if (result == NULL) {
        fprintf(stderr, "Memory allocation failed for result matrix.\n");
        exit(EXIT_FAILURE);
    }

    // Initial guess: roots uniformly spaced on the unit circle
    double pi = 3.14159265358979323846;
    for (unsigned int i = 0; i < n; i++) {
        double angle = 2 * pi * i / n; // Angle for the unit circle
        gsl_matrix_set(result, i, 0, cos(angle)); // Real part
        gsl_matrix_set(result, i, 1, sin(angle)); // Imaginary part
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
            real_i = gsl_matrix_get(result, i, 0);
            imag_i = gsl_matrix_get(result, i, 1);

            // Compute the value of the polynomial and its correction term
            p_real = gsl_vector_get(coeff, 0); // Leading coefficient
            p_imag = 0.0; // No imaginary part initially
            denom_real = 1.0, denom_imag = 0.0; // Product term denominator

            for (unsigned int j = 0; j < n; j++) {
                if (j == i) continue; // Skip the current root
                real_j = gsl_matrix_get(result, j, 0);
                imag_j = gsl_matrix_get(result, j, 1);

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
                coeff_real = gsl_vector_get(coeff, k);
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

            gsl_matrix_set(result, i, 0, new_real);
            gsl_matrix_set(result, i, 1, new_imag);

            // Check for convergence
            if (fabs(delta_real) > tolerance || fabs(delta_imag) > tolerance) {
                converged = 0;
            }
        }

        if (converged) break; // Stop if all roots have converged
    }
    return result;
}

//////////////////////////////////////////////////////////////////////////////////////

void freeVector(gsl_vector *vector) {
    if (vector != NULL) {
        gsl_vector_free(vector); // Free the GSL vector
    }
}

gsl_vector* zerosVector(int size) {
    return gsl_vector_calloc(size); // Allocate and initialize with zeros
}

gsl_vector* createVector(int size, double *values) {
    gsl_vector *vector = gsl_vector_alloc(size);

    if (vector == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Populate the vector with the provided values
    for (int i = 0; i < size; i++) {
        gsl_vector_set(vector, i, values[i]);
    }

    return vector;
}

void printVector(const gsl_vector *vector) {
    size_t size = vector->size;

    printf("[");

    if (size > 10) {
        // Print the first three elements
        for (size_t i = 0; i < 3; i++) {
            printf("%.4f", gsl_vector_get(vector, i));
            if (i < 2) {
                printf(", ");
            }
        }
        // Print ellipsis
        printf(", ...");
        // Print the last three elements
        for (size_t i = size - 3; i < size; i++) {
            printf(", %.4f", gsl_vector_get(vector, i));
        }
    } else {
        // Print all elements if size is 10 or less
        for (size_t i = 0; i < size; i++) {
            printf("%.4f", gsl_vector_get(vector, i));
            if (i < size - 1) {
                printf(", ");
            }
        }
    }

    printf("]\n");
}

gsl_vector* matvec(const gsl_matrix *A, const gsl_vector *x) {
    // Check if the number of columns in A matches the size of the vector x
    if (A->size2 != x->size) {
        fprintf(stderr, "Error: Number of columns in matrix A must match the size of vector x.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the result vector
    gsl_vector *result = gsl_vector_alloc(A->size1);
    if (result == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Perform matrix-vector multiplication using GSL
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, result);

    return result;
}

gsl_vector* linsolve(const gsl_matrix *A, const gsl_vector *b) {
    // Check for dimension compatibility
    if (A == NULL || b == NULL) {
        fprintf(stderr, "Error: Null pointer input.\n");
        exit(EXIT_FAILURE);
    }
    if (A->size2 != b->size) {
        fprintf(stderr, "Error: Incompatible dimensions for system.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the result vector
    gsl_vector *x = gsl_vector_alloc(b->size);
    if (x == NULL) {
        fprintf(stderr, "Memory allocation failed for result vector.\n");
        exit(EXIT_FAILURE);
    }

    gsl_matrix *tempA = gsl_matrix_alloc(A->size1, A->size2);
    gsl_vector *tempb = gsl_vector_alloc(b->size);
    gsl_matrix_memcpy(tempA, A);
    gsl_vector_memcpy(tempb, b);

    gsl_permutation *p = gsl_permutation_alloc(A->size1);
    int signum;

    gsl_linalg_LU_decomp(tempA, p, &signum);
    gsl_linalg_LU_solve(tempA, p, tempb, x);

    gsl_matrix_free(tempA);
    gsl_vector_free(tempb);
    gsl_permutation_free(p);

    return x; // Return the result vector
}

gsl_vector* vecadd(const gsl_vector *A, const gsl_vector *B) {
    if (A->size != B->size) {
        fprintf(stderr, "Vector dimensions are incompatible for addition.\n");
        exit(EXIT_FAILURE);
    }

    gsl_vector *C = gsl_vector_alloc(A->size);
    if (C == NULL) {
        fprintf(stderr, "Memory allocation failed for result vector.\n");
        exit(EXIT_FAILURE);
    }

    // Copy B into C
    gsl_vector_memcpy(C, B);

    // Perform vector addition
    gsl_vector_add(C, A);

    return C;
}

gsl_vector* vecscal(const gsl_vector *A, double scalar) {
    gsl_vector *C = gsl_vector_alloc(A->size);

    // Copy A into C
    gsl_vector_memcpy(C, A);

    // Scale the vector by a constant
    gsl_vector_scale(C, scalar);

    return C;
}

gsl_vector* vecrand(int dim) {
    gsl_vector *vec = gsl_vector_alloc(dim);

    if (vec == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    for (int i = 0; i < dim; i++) {
        gsl_vector_set(vec, i, gsl_rng_uniform(rng));
    }

    gsl_rng_free(rng);

    return vec;
}

gsl_vector* vecelem(const gsl_vector *A, const gsl_vector *B) {
    // Check if both vectors have the same dimensions
    if (A->size != B->size) {
        fprintf(stderr, "Vectors must have the same dimensions for element-wise multiplication.\n");
        exit(EXIT_FAILURE);
    }

    gsl_vector *C = gsl_vector_alloc(A->size);

    if (C == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Perform element-wise multiplication
    for (size_t i = 0; i < C->size; i++) {
        double value = gsl_vector_get(A, i) * gsl_vector_get(B, i);
        gsl_vector_set(C, i, value);
    }

    return C;
}

double norm(const gsl_vector *vec, double p) {
    double norm_value = 0.0;

    if (p == 1) {
        // 1-norm using GSL
        norm_value = gsl_blas_dasum(vec);
    } else if (p == 2) {
        // 2-norm using GSL
        norm_value = gsl_blas_dnrm2(vec);
    } else if (p == INFINITY) {
        // Infinity-norm (Maximum norm)
        for (size_t i = 0; i < vec->size; i++) {
            double abs_val = fabs(gsl_vector_get(vec, i));
            if (abs_val > norm_value) {
                norm_value = abs_val;
            }
        }
    } else {
        // p-norm
        for (size_t i = 0; i < vec->size; i++) {
            norm_value += pow(fabs(gsl_vector_get(vec, i)), p);
        }
        norm_value = pow(norm_value, 1.0 / p);
    }

    return norm_value;
}

gsl_vector* createArray(double start, double end, double step) {
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
    gsl_vector *result = gsl_vector_alloc(size);
    if (result == NULL) {
        fprintf(stderr, "Memory allocation failed for result array.\n");
        exit(EXIT_FAILURE);
    }

    // Fill the array with values
    double current = start;
    for (int i = 0; i < size; i++) {
        gsl_vector_set(result, i, current);
        current += step;
        if ((step > 0 && current > end) || (step < 0 && current < end)) {
            current = end; // Ensure last value is exactly the endpoint
        }
    }

    return result;
}

gsl_vector* polyfitweighted(const gsl_vector *x, const gsl_vector *y, const gsl_vector *w, unsigned int n) {
    unsigned int len = x->size;
    unsigned int cols = n + 1;

    gsl_matrix *V = gsl_matrix_alloc(len, cols);
    gsl_vector *wy = gsl_vector_alloc(len);

    if (V == NULL || wy == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    // Construct the weighted Vandermonde matrix
    for (unsigned int i = 0; i < len; i++) {
        // Set the last column as the weight
        gsl_matrix_set(V, i, n, gsl_vector_get(w, i));

        // Fill the Vandermonde row starting from the highest power to the lowest
        for (int j = n - 1; j >= 0; j--) {
            gsl_matrix_set(V, i, j, gsl_vector_get(x, i) * gsl_matrix_get(V, i, j + 1));
        }

        // Calculate weighted y values
        gsl_vector_set(wy, i, gsl_vector_get(w, i) * gsl_vector_get(y, i));
    }

    // Perform SVD
    gsl_matrix *U = gsl_matrix_alloc(len, cols);
    gsl_matrix *V_svd = gsl_matrix_alloc(cols, cols);
    gsl_vector *S = gsl_vector_alloc(cols);
    gsl_vector *work = gsl_vector_alloc(cols);

    gsl_matrix_memcpy(U, V);
    gsl_linalg_SV_decomp(U, V_svd, S, work);

    // Solve for p using SVD
    gsl_vector *p = gsl_vector_alloc(cols);
    gsl_linalg_SV_solve(U, V_svd, S, wy, p);

    // Free the allocated matrices and vectors used for SVD
    gsl_matrix_free(U);
    gsl_matrix_free(V_svd);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_matrix_free(V);
    gsl_vector_free(wy);

    return p;
}

gsl_vector* linspace(double start, double end, int num) {
    gsl_vector *result = gsl_vector_alloc(num);

    if (result == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    double step = (end - start) / (num - 1);

    for (int i = 0; i < num; i++) {
        gsl_vector_set(result, i, start + i * step);
    }

    return result;
}

gsl_vector* polycoefs(const gsl_vector *roots) {
    int n = roots->size; // Number of roots
    gsl_vector *coeff = gsl_vector_alloc(n + 1); // Coefficients array size is degree + 1

    if (coeff == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Initialize all coefficients to zero
    gsl_vector_set_zero(coeff);

    gsl_poly_complex_workspace *workspace = gsl_poly_complex_workspace_alloc(n + 1);
    gsl_vector_set(coeff, 0, 1.0); // Leading coefficient (x^n)

    for (int i = 0; i < n; i++) {
        double root = gsl_vector_get(roots, i);

        for (int j = i; j >= 0; j--) {
            double new_val = gsl_vector_get(coeff, j + 1) - root * gsl_vector_get(coeff, j);
            gsl_vector_set(coeff, j + 1, new_val);
        }
    }

    gsl_poly_complex_workspace_free(workspace);

    return coeff;
}

gsl_vector* conv(const gsl_vector *a, const gsl_vector *b) {
    int n = a->size + b->size - 1; // Size of the result vector
    gsl_vector *result = gsl_vector_alloc(n);

    if (result == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Perform convolution
    gsl_vector_set_zero(result);

    for (size_t i = 0; i < a->size; i++) {
        for (size_t j = 0; j < b->size; j++) {
            gsl_vector_set(result, i + j,
                           gsl_vector_get(result, i + j) +
                           gsl_vector_get(a, i) * gsl_vector_get(b, j));
        }
    }

    return result;
}

gsl_vector* interp1(const gsl_vector *s, const gsl_vector *t, const gsl_vector *ss, const char *method) {
    // Check for dimension compatibility
    if (s == NULL || t == NULL || ss == NULL || s->size != t->size) {
        fprintf(stderr, "Error: Input vectors are incompatible.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the result vector
    gsl_vector *result = gsl_vector_alloc(ss->size);
    if (result == NULL) {
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
        gsl_vector_free(result);
        exit(EXIT_FAILURE);
    }

    // Create a GSL interpolation object and accelerator
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(interp_type, s->size);

    // Initialize the interpolation object
    if (gsl_spline_init(spline, gsl_vector_const_ptr(s, 0), gsl_vector_const_ptr(t, 0), s->size) != GSL_SUCCESS) {
        fprintf(stderr, "Error: Failed to initialize interpolation.\n");
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
        gsl_vector_free(result);
        exit(EXIT_FAILURE);
    }

    // Perform interpolation for each value in ss
    for (size_t i = 0; i < ss->size; i++) {
        double ss_val = gsl_vector_get(ss, i);
        if (ss_val < gsl_vector_get(s, 0) || ss_val > gsl_vector_get(s, s->size - 1)) {
            fprintf(stderr, "Error: Interpolation point ss[%zu] = %f is out of bounds.\n", i, ss_val);
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            gsl_vector_free(result);
            exit(EXIT_FAILURE);
        }
        gsl_vector_set(result, i, gsl_spline_eval(spline, ss_val, acc));
    }

    // Free GSL resources
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return result; // Return the interpolated vector
}

//////////////////////////////////////////////////////////////////////////////////////

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

ode_solver* create_ode_solver(int (*func)(double, const double[], double[], void *), void *params, size_t dim, double hstart, double epsabs, double epsrel) {
    ode_solver *solver = (ode_solver *)malloc(sizeof(ode_solver));
    solver->system.function = func;
    solver->system.jacobian = NULL;
    solver->system.dimension = dim;
    solver->system.params = params;
    solver->driver = gsl_odeiv2_driver_alloc_y_new(&solver->system, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);
    return solver;
}

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

void free_ode_solver(ode_solver *solver) {
    gsl_odeiv2_driver_free(solver->driver);
    free(solver);
}

void ode45(int (*func)(double, const double[], double[], void *), void *params, size_t dimension, double *t, double t1, double y[], double dt) {
    // Default parameters for the solver
    double initial_step_size = 1e-6;
    double absolute_error = 1e-6;
    double relative_error = 0.0;

    ode_solver *solver = create_ode_solver(func, params, dimension, initial_step_size, absolute_error, relative_error);
    solve_ode(solver, t, t1, y, dt);
    free_ode_solver(solver);
}

