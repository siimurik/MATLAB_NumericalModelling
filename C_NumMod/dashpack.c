/***********************************************************************************/
/*    DASHPACK - Double-precision Algebraic Solutions for High-performance PACKage */
/***********************************************************************************/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>
#include <cblas.h>

typedef struct{
    double *data;
    int    rows;
    int    cols;
} Matrix;


typedef struct{
    double *data;
    int    size;
} Vector;

/* Function to free the allocated memory
--------------------------------------------------------
 You can also use:
 ```
    void matvecFree(MatrixVector matvec) {
        free(matvec.data);
    }
```
 But it can lead to potential memory leaks if you 
 need to manage the struct after freeing its data.

 This one below is more robust and aligns with 
 best practices in C programming for memory management.*/
void freeMatrix(Matrix *matrix) {
    free(matrix->data);
    matrix->data = NULL;    // Avoid dangling pointer.
}                           // Good practice to avoid accidental
                            // access to freed memory.

// Function to create and initialize a matrix with specific values
Matrix createMatrix(int rows, int cols, double *values) {
    Matrix matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.data = (double *)malloc(rows * cols * sizeof(double)); // Allocate memory

    if (matrix.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Populate the matrix with the provided values
    memcpy(matrix.data, values, matrix.rows * matrix.cols * sizeof(double));
    
    // Perform in parallel if necessary
    //#pragma omp parallel for    
    //for (int i = 0; i < rows * cols; i++) {
    //    matrix.data[i] = values[i];
    //}
    

    return matrix;
}

// Function for printing a matrix
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

// Function to find the transposed matrix
Matrix transposeMatrix(const Matrix *mat) {
    Matrix transposed;
    transposed.rows = mat->cols;
    transposed.cols = mat->rows;
    transposed.data = (double *)malloc(transposed.rows * transposed.cols * sizeof(double));

    if (transposed.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Perform transposition in parallel if necessary
    //#pragma omp parallel for
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            transposed.data[j * transposed.rows + i] = mat->data[i * mat->cols + j];
        }
    }

    return transposed;
}

// Function to compute the inverse matrix using LAPACK
Matrix inverseMatrix(const Matrix *mat) {
    if (mat->rows != mat->cols) {
        fprintf(stderr, "Matrix must be square to compute the inverse.\n");
        exit(EXIT_FAILURE);
    }

    int n = mat->rows;
    int *ipiv = (int *)malloc(n * sizeof(int)); // Pivot indices
    Matrix inverse;
    inverse.data = (double *)malloc(n * n * sizeof(double));
    inverse.rows = n;
    inverse.cols = n;

    if (inverse.data == NULL || ipiv == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Copy the original matrix to the inverse matrix
    memcpy(inverse.data, mat->data, inverse.rows * inverse.cols * sizeof(double));
    
    // Extra option for parallelization:
    // If not done this way, it will run, but dump the core and say
    // "free(): double free detected in tcache 2
    // Aborted (core dumped)"
    // #pragma omp parallel for
    //for (int i = 0; i < n; i++) {
    //    for (int j = 0; j < n; j++) {
    //        inverse.data[i * n + j] = mat->data[i * n + j];
    //    }
    //}


    // Perform LU decomposition
    int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, inverse.data, n, ipiv);
    if (info != 0) {
        fprintf(stderr, "Matrix is singular and cannot be inverted.\n");
        free(inverse.data);
        free(ipiv);
        exit(EXIT_FAILURE);
    }

    // Compute the inverse
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, inverse.data, n, ipiv);
    if (info != 0) {
        fprintf(stderr, "Error occurred during inversion.\n");
        free(inverse.data);
        free(ipiv);
        exit(EXIT_FAILURE);
    }

    free(ipiv);
    return inverse;
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
    C.data = (double *)malloc(C.rows * C.cols * sizeof(double));

    // Perform matrix multiplication
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, C.rows, C.cols, A->cols, 
                    1.0, A->data, A->cols, B->data, B->cols, 0.0, C.data, C.rows);

    return C;
}

// Matrix addition function
Matrix matadd(const Matrix *A, const Matrix *B) {
    if (A->rows != B->rows || A->cols != B->cols) {
        fprintf(stderr, "Matrix dimensions are incompatible for addition.\n");
        exit(EXIT_FAILURE);
    }

    Matrix C;
    C.rows = A->rows;
    C.cols = A->cols;
    C.data = (double *)malloc(C.rows * C.cols * sizeof(double));
    memcpy(C.data, B->data, C.rows * C.cols * sizeof(double));

    // Constant times a vector plus a vector
    // ALPHA * A + C;   ALPHA = 1.0;
    cblas_daxpy(C.rows * C.cols, 1.0, A->data, 1, C.data, 1);

    return C;
}

// Matrix scalar multiplication function
Matrix matscal(const Matrix *A, double scalar) {
    Matrix C;
    C.rows = A->rows;
    C.cols = A->cols;
    C.data = (double *)malloc(C.rows * C.cols * sizeof(double));
    memcpy(C.data, A->data, C.rows * C.cols * sizeof(double));

    // Scales a vector by a constant
    cblas_dscal(C.rows * C.cols, scalar, C.data, 1);

    return C;
}

// Function to compute the determinant of a matrix using LU-decomposition
double det(const Matrix *mat) {
    // Check if the matrix is square
    if (mat->rows != mat->cols) {
        fprintf(stderr, "Matrix must be square to compute the determinant.\n");
        exit(EXIT_FAILURE);
    }

    int n = mat->rows;
    int *ipiv = (int *)malloc(n * sizeof(int)); // Pivot indices
    double determinant = 1.0;

    if (ipiv == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Need to make a copy of 'mat', bc sgetrf() modifies the input matrix
    Matrix matCopy;
    matCopy.rows = mat->rows;
    matCopy.cols = mat->cols;
    matCopy.data = (double *)malloc(matCopy.rows * matCopy.cols * sizeof(double));

    if (matCopy.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(ipiv);
        exit(EXIT_FAILURE);
    }    

    // Copy the original matrix values to the temporary one to avoid overwriting
    memcpy(matCopy.data, mat->data, matCopy.rows * matCopy.cols * sizeof(double));

    // Perform in-place LU decomposition
    int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, matCopy.data, n, ipiv);
    if (info != 0) {
        fprintf(stderr, "Matrix is singular and cannot compute the determinant.\n");
        free(ipiv);
        exit(EXIT_FAILURE);
    }

    // Calculate the determinant as the product of the diagonal elements
    for (int i = 0; i < n; i++) {
        determinant *= matCopy.data[i * n + i]; // Diagonal element
        if (ipiv[i] != i + 1) { // Check for row swaps
            determinant = -determinant; // Adjust sign for row swaps
        }
    }
    
    // Free temporary allocated memory, which is not necessary
    // to keep outside the function.
    free(ipiv);
    freeMatrix(&matCopy);
    return determinant;
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
    C.data = (double *)malloc(C.rows * C.cols * sizeof(double));

    if (C.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Perform element-wise multiplication 
    // Use OpenMP for parallelization in necessary
    //#pragma omp parallel for
    for (int i = 0; i < C.rows; i++) {
        for (int j = 0; j < C.cols; j++) {
            C.data[i * C.cols + j] = A->data[i * A->cols + j] * B->data[i * B->cols + j];
        }
    }

    return C;
}

// Function to solve an overdetermined system of linear equations
Matrix linsolve_overdet(const Matrix *A, const Matrix *F) {
    // Prepare the LAPACK parameters
    int n = A->rows;  // Size of the matrix
    int nrhs = F->cols;  // Number of right-hand sides
    int lda = n;  // Leading dimension of A
    int ldb = n;  // Leading dimension of B
    int ipiv[n];  // Array for pivot indices
    int info;  // Variable to store the info from LAPACK

    // Ensure x has the correct dimensions
    Matrix x;
    x.rows = n;
    x.cols = nrhs;
    x.data = (double *)malloc(n * nrhs * sizeof(double));
    
    if (x.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Write F matrix values into x, since they get overwritten
    memcpy(x.data, F->data, x.rows * x.cols * sizeof(double));
    //for (int i = 0; i < n; i++) {
    //    for (int j = 0; j < nrhs; j++) {
    //        x.data[i * nrhs + j] = F->data[i * nrhs + j];  // Copy each element
    //    }
    //}

    // Solve the system using LAPACK's sgesv function
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A->data, lda, ipiv, x.data, ldb);

    if (info != 0) {
        fprintf(stderr, "Error: LAPACK sgesv failed with info = %d\n", info);
        free(x.data);
        exit(EXIT_FAILURE);
    }

    return x;
}

// Function to generate a new Matrix with random values
Matrix matrand(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (double *)malloc(rows * cols * sizeof(double));

    double min = 0.0;
    double max = 1.0;
	for (int i = 0; i < rows * cols; i++) {
	    mat.data[i] = min + (max - min) * ((double)rand() / (double)RAND_MAX);
	}

    return mat;
}

// Function to create the companion matrix:
Matrix compan(const double *coefficients, int size) {
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
    C.data = (double *)calloc(C.rows * C.cols, sizeof(double));

    if (C.data == NULL) {
        fprintf(stderr, "Memory allocation failed for companion matrix.\n");
        return (Matrix){NULL, 0, 0}; // Return an empty matrix
    }

    // Fill the first row
    double leading_coefficient = coefficients[0];
    for (int j = 0; j < n; j++) {
        C.data[0 * C.cols + j] = -coefficients[j + 1] / leading_coefficient; // First row
    }

    // Fill the sub-diagonal with ones
    for (int i = 1; i < n; i++) {
        C.data[i * C.cols + (i - 1)] = 1.0;
    }

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
    double *A    = (double *)malloc(LDA * N * sizeof(double));
    double *WR   = (double *)malloc(N * sizeof(double)); // Real parts of eigenvalues
    double *WI   = (double *)malloc(N * sizeof(double)); // Imaginary parts of eigenvalues
    int   LDVL  = matrix->rows;
    int   LDVR  = matrix->rows;
    double *VR   = (double *)malloc(LDVR * N * sizeof(double)); // Right eigenvectors
    int INFO;

    if (A == NULL || WR == NULL || WI == NULL || VR == NULL) {
        fprintf(stderr, "Failed to allocate memory.\n");
        goto cleanup;
    }

    // Copy matrix data to LAPACK format
    memcpy(A, matrix->data, N * N * sizeof(double));

    // Compute eigenvalues and right eigenvectors
    //INFO = LAPACKE_dgeev(LAPACK_COL_MAJOR, JOBVL, JOBVR, N, A, LDA, WR, WI, NULL, N, VR, N);
	INFO = LAPACKE_dgeev( LAPACK_COL_MAJOR, JOBVL, JOBVR, N, A, LDA, WR, WI, NULL, LDVL, VR, LDVR);

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

//=====================================================================================
//=====================================================================================
//=====================================================================================

// Functions for vectors

// Function to free the allocated memory in a vector
void freeVector(Vector *vector) {
    free(vector->data);
    vector->data = NULL; // Avoid dangling pointer
}

// Function to create a vector from array values
Vector createVector(int size, double *values) {
    Vector vector;
    vector.size = size;
    vector.data = (double *)malloc(size * sizeof(double)); // Allocate memory

    if (vector.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Populate the vector with the provided values
    memcpy(vector.data, values, vector.size * sizeof(double));
    
    //for (int i = 0; i < size; i++) {
    //    vector.data[i] = values[i];
    //}


    return vector;
}

void printVector(const Vector *vector) {
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

// Function to perform matrix-vector multiplication
Vector matvec(const Matrix *A, const Vector *x) {
    // Check if the number of columns in A matches the size of the vector x
    if (A->cols != x->size) {
        fprintf(stderr, "Error: Number of columns in matrix A must match the size of vector x.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the result vector
    Vector result;
    result.size = A->rows;
    result.data = (double *)malloc(result.size * sizeof(double));
    if (result.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Perform matrix-vector multiplication using BLAS
    // result = alpha * A * x + beta * result
    double alpha = 1.0f; // Scalar multiplier for A * x
    double beta = 0.0f;  // Scalar multiplier for the initial value of result
    cblas_dgemv(CblasRowMajor, CblasNoTrans, A->rows, A->cols, alpha, A->data, 
                                    A->cols, x->data, 1, beta, result.data, 1);

    return result;
}

// Function to solve the linear system Ax = b
Vector linsolve(const Matrix *A, const Vector *b) {
    // Check for dimension compatibility
    if (A->cols != b->size) {
        fprintf(stderr, "Error: Incompatible dimensions for system.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the result vector
    Vector x;
    x.size = A->cols;
    x.data = (double *)malloc(x.size * sizeof(double));

    if (x.data == NULL) {
        fprintf(stderr, "Memory allocation failed for result vector.\n");
        exit(EXIT_FAILURE);
    }

    // Copy the input vector b to the result vector x
    memcpy(x.data, b->data, x.size * sizeof(double));
    //for (int i = 0; i < x.size; i++) {
    //    x.data[i] = b->data[i];
    //}

    // Perform LU factorization and solve the system using LAPACKE
    int *ipiv = (int *)malloc(A->rows * sizeof(int));
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, // Matrix storage order
                             A->rows, // Number of rows in A
                             1, // Number of right-hand sides
                             A->data, // Pointer to matrix A
                             A->cols, // Leading dimension of A
                             ipiv, // Pivot indices
                             x.data, // Pointer to result vector x
                             A->cols); // Leading dimension of x

    if (info != 0) {
        fprintf(stderr, "Error: LAPACKE_dgesv failed with info = %d.\n", info);
        exit(EXIT_FAILURE);
    }

    // Free allocated memory for pivot indices
    free(ipiv);

    return x; // Return the result vector
}

// Vector addition function
Vector vecadd(const Vector *A, const Vector *B) {
    if (A->size != B->size) {
        fprintf(stderr, "Vector dimensions are incompatible for addition.\n");
        exit(EXIT_FAILURE);
    }

    Vector C;
    C.size = A->size;
    C.data = (double *)malloc(C.size * sizeof(double));
    memcpy(C.data, B->data, C.size * sizeof(double));

    // Constant times a vector plus a vector
    // ALPHA * A + C;   ALPHA = 1.0;
    cblas_saxpy(C.size, 1.0, A->data, 1, C.data, 1);

    return C;
}

// Vector scalar multiplication function
Vector vecscal(const Vector *A, double scalar) {
    Vector C;
    C.size = A->size;
    C.data = (double *)malloc(C.size * sizeof(double));
    memcpy(C.data, A->data, C.size * sizeof(double));

    // Scales a vector by a constant
    cblas_sscal(C.size, scalar, C.data, 1);

    return C;
}

// Generate a vectro full of random values ranging from 0...1
Vector vecrand(int dim){
	Vector vec;
	vec.size = dim;
	vec.data = (double *)malloc(dim * sizeof(double));

    double min = 0.0;
    double max = 1.0;
	for (int i = 0; i < vec.size; i++) {
	    vec.data[i] = min + (max - min) * ((double)rand() / (double)RAND_MAX);
	}

    return vec;
}