/****************************************************************************/
/*  Officaly known as                                                       */
/*    SPLASH - Single Precision for Linear Algebra and Scientific Handling  */
/*--------------------------------------------------------------------------*/
/*  Or unoffically known as                                                 */
/*    SPLASH - Siim’s Package for Linear Algebra and Scientific Handling    */
/****************************************************************************/
/*  
    $ gcc -c splash.c -o splash.o
    $ ar rcs libsplash.a splash.o
*/
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include <omp.h> 

typedef struct{
    float *data;
    unsigned int rows;
    unsigned int cols;
} Matrix;


typedef struct{
    float *data;
    unsigned int size;
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
Matrix createMatrix(int rows, int cols, float *values) {
    Matrix matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.data = (float *)malloc(rows * cols * sizeof(float)); // Allocate memory

    if (matrix.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Populate the matrix with the provided values
    memcpy(matrix.data, values, matrix.rows * matrix.cols * sizeof(float));
    
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
    transposed.data = (float *)malloc(transposed.rows * transposed.cols * sizeof(float));

    if (transposed.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    int i, j;
    // Perform transposition in parallel if necessary
    //#pragma omp parallel for
    for (i = 0; i < mat->rows; i++) {
        for (j = 0; j < mat->cols; j++) {
            transposed.data[j * transposed.cols + i] = mat->data[i * mat->cols + j];
        }
    }

    return transposed;
}

// Function to find the transposed matrix
Matrix transposeMatrixPara(const Matrix *mat) {
    Matrix transposed;
    transposed.rows = mat->cols;
    transposed.cols = mat->rows;
    transposed.data = (float *)malloc(transposed.rows * transposed.cols * sizeof(float));

    if (transposed.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    int i, j;
    // Perform transposition in parallel on the rows
    #pragma omp parallel for private(j) 
    for (i = 0; i < mat->rows; i++) {
        for (j = 0; j < mat->cols; j++) {
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
    inverse.data = (float *)malloc(n * n * sizeof(float));
    inverse.rows = n;
    inverse.cols = n;

    if (inverse.data == NULL || ipiv == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Copy the original matrix to the inverse matrix
    memcpy(inverse.data, mat->data, inverse.rows * inverse.cols * sizeof(float));
    
    // Extra option for parallelization:
    // If not done this way, it will run, but dump the core and say
    // "free(): float free detected in tcache 2
    // Aborted (core dumped)"
    // #pragma omp parallel for
    //for (int i = 0; i < n; i++) {
    //    for (int j = 0; j < n; j++) {
    //        inverse.data[i * n + j] = mat->data[i * n + j];
    //    }
    //}

    // Perform LU decomposition
    int info = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, n, n, inverse.data, n, ipiv);
    if (info != 0) {
        fprintf(stderr, "Matrix is singular and cannot be inverted.\n");
        free(inverse.data);
        free(ipiv);
        exit(EXIT_FAILURE);
    }

    // Compute the inverse
    info = LAPACKE_sgetri(LAPACK_ROW_MAJOR, n, inverse.data, n, ipiv);
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
    C.data = (float *)malloc(C.rows * C.cols * sizeof(float));

    // Perform matrix multiplication
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, C.rows, C.cols, A->cols, 
                    1.0f, A->data, A->cols, B->data, B->cols, 0.0f, C.data, C.rows);

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
    C.data = (float *)malloc(C.rows * C.cols * sizeof(float));
    memcpy(C.data, B->data, C.rows * C.cols * sizeof(float));

    // Constant times a vector plus a vector
    // ALPHA * A + C;   ALPHA = 1.0;
    cblas_saxpy(C.rows * C.cols, 1.0f, A->data, 1, C.data, 1);

    return C;
}

// Matrix scalar multiplication function
Matrix matscal(const Matrix *A, float scalar) {
    Matrix C;
    C.rows = A->rows;
    C.cols = A->cols;
    C.data = (float *)malloc(C.rows * C.cols * sizeof(float));
    memcpy(C.data, A->data, C.rows * C.cols * sizeof(float));

    // Scales a vector by a constant
    cblas_sscal(C.rows * C.cols, scalar, C.data, 1);

    return C;
}

// Function to compute the determinant of a matrix using LU-decomposition
float det(const Matrix *mat) {
    // Check if the matrix is square
    if (mat->rows != mat->cols) {
        fprintf(stderr, "Matrix must be square to compute the determinant.\n");
        exit(EXIT_FAILURE);
    }

    int n = mat->rows;
    int *ipiv = (int *)malloc(n * sizeof(int)); // Pivot indices
    float determinant = 1.0f;

    if (ipiv == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Need to make a copy of 'mat', bc sgetrf() modifies the input matrix
    Matrix matCopy;
    matCopy.rows = mat->rows;
    matCopy.cols = mat->cols;
    matCopy.data = (float *)malloc(matCopy.rows * matCopy.cols * sizeof(float));

    if (matCopy.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(ipiv);
        exit(EXIT_FAILURE);
    }    

    // Copy the original matrix values to the temporary one to avoid overwriting
    memcpy(matCopy.data, mat->data, matCopy.rows * matCopy.cols * sizeof(float));

    // Perform in-place LU decomposition
    int info = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, n, n, matCopy.data, n, ipiv);
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
    C.data = (float *)malloc(C.rows * C.cols * sizeof(float));

    if (C.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Perform element-wise multiplication 
    for (int i = 0; i < C.rows; i++) {
        for (int j = 0; j < C.cols; j++) {
            C.data[i * C.cols + j] = A->data[i * A->cols + j] * B->data[i * B->cols + j];
        }
    }

    return C;
}
///*
// Function to perform matrix element-wise multiplication
Matrix matelemPara(const Matrix *A, const Matrix *B) {
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
    // Use OpenMP for parallelization for HUGE matrices
    #pragma omp parallel for default(none) shared(A, B, C) private(i, j) collapse(2)
    for (i = 0; i < C.rows; i++) {
        for (j = 0; j < C.cols; j++) {
            C.data[i * C.cols + j] = A->data[i * A->cols + j] * B->data[i * B->cols + j];
        }
    }

    return C;
}
//*/
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
    x.data = (float *)malloc(n * nrhs * sizeof(float));
    
    if (x.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Write F matrix values into x, since they get overwritten
    memcpy(x.data, F->data, x.rows * x.cols * sizeof(float));
    //for (int i = 0; i < n; i++) {
    //    for (int j = 0; j < nrhs; j++) {
    //        x.data[i * nrhs + j] = F->data[i * nrhs + j];  // Copy each element
    //    }
    //}


    // Solve the system using LAPACK's sgesv function
    info = LAPACKE_sgesv(LAPACK_ROW_MAJOR, n, nrhs, A->data, lda, ipiv, x.data, ldb);

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
    mat.data = (float *)malloc(rows * cols * sizeof(float));

    if (mat.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    float min = 0.0f;
    float max = 1.0f;
	for (int i = 0; i < rows * cols; i++) {
	    mat.data[i] = min + (max - min) * ((float)rand() / (float)RAND_MAX);
	}

    return mat;
}
///*
// Function to generate a new Matrix with random values
Matrix matrandPara(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (float *)malloc(rows * cols * sizeof(float));

    if (mat.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    float min = 0.0f;
    float max = 1.0f;
    // Parallelize the random number generation
    #pragma omp parallel
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
//*/
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
        C.data[i * C.cols + (i - 1)] = 1.0f;
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

// Calculate any norm of a matrix using BLAS
float normMatrix(const Matrix *mat, float p)
{
    float norm_value = 0.0f;

    if (p == 1) {
        // 1-norm (maximum column sum)
        for (int j = 0; j < mat->cols; j++) 
        {
            float col_sum = 0.0f;
            for (int i = 0; i < mat->rows; i++) 
            {
                col_sum += fabs(mat->data[i * mat->cols + j]);
            }
            if (col_sum > norm_value) 
            {
                norm_value = col_sum;
            }
        }
    } else if (p == 2) 
    {
        // Frobenius norm (similar to Euclidean norm for vectors)
        norm_value = cblas_snrm2(mat->rows * mat->cols, mat->data, 1);
    } else if (p == INFINITY) 
    {
        // Infinity-norm (maximum row sum)
        for (int i = 0; i < mat->rows; i++) 
        {
            float row_sum = 0.0f;
            for (int j = 0; j < mat->cols; j++) 
            {
                row_sum += fabs(mat->data[i * mat->cols + j]);
            }
            if (row_sum > norm_value) 
            {
                norm_value = row_sum;
            }
        }
    } else 
    {
        fprintf(stderr, "Unsupported norm type for matrices.\n");
        exit(EXIT_FAILURE);
    }

    return norm_value;
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
    result.data = (float *)malloc(result.size * sizeof(float));
    if (result.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Perform matrix-vector multiplication using BLAS
    // result = alpha * A * x + beta * result
    float alpha = 1.0f; // Scalar multiplier for A * x
    float beta = 0.0f;  // Scalar multiplier for the initial value of result
    cblas_sgemv(CblasRowMajor, CblasNoTrans, A->rows, A->cols, alpha, A->data, 
                                    A->cols, x->data, 1, beta, result.data, 1);

    return result;
}

// Function to solve the linear system Ax = b
Vector linsolve(const Matrix *A, const Vector *b) {
    // Check for dimension compatibility
    if (A == NULL || b == NULL || A->data == NULL || b->data == NULL) {
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
    x.data = (float *)malloc(x.size * sizeof(float));
    if (x.data == NULL) {
        fprintf(stderr, "Memory allocation failed for result vector.\n");
        exit(EXIT_FAILURE);
    }

    // Copy the input vector b to the result vector x
    memcpy(x.data, b->data, x.size * sizeof(float));

    // Allocate memory for pivot indices
    int *ipiv = (int *)malloc(A->rows * sizeof(int));
    if (ipiv == NULL) {
        fprintf(stderr, "Memory allocation failed for pivot indices.\n");
        free(x.data);
        exit(EXIT_FAILURE);
    }

    // Solve the system using LAPACKE_sgesv
    int info = LAPACKE_sgesv(LAPACK_ROW_MAJOR,  // Row-major order
                             A->rows,          // Number of rows in A
                             1,                // Number of right-hand sides
                             A->data,          // Pointer to matrix A
                             A->cols,          // Leading dimension of A
                             ipiv,             // Pivot indices
                             x.data,           // Pointer to result vector x
                             1);               // Leading dimension of x

    if (info != 0) {
        fprintf(stderr, "Error: LAPACKE_sgesv failed with info = %d.\n", info);
        free(x.data);
        free(ipiv);
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
    C.data = (float *)malloc(C.size * sizeof(float));
    memcpy(C.data, B->data, C.size * sizeof(float));

    // Constant times a vector plus a vector
    // ALPHA * A + C;   ALPHA = 1.0;
    cblas_saxpy(C.size, 1.0f, A->data, 1, C.data, 1);

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

// Create a vector full of random values
Vector vecrand(int dim) {
    Vector vec;
    vec.size = dim;
    vec.data = (float *)malloc(dim * sizeof(float));

    if (vec.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    float min = 0.0f;
    float max = 1.0f;
	for (int i = 0; i < vec.size; i++) {
	    vec.data[i] = min + (max - min) * ((float)rand() / (float)RAND_MAX);
	}

    return vec;
}
///*
// Create a vector full of random values in parallel
Vector vecrandPara(int dim) {
    Vector vec;
    vec.size = dim;
    vec.data = (float *)malloc(dim * sizeof(float));

    if (vec.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    float min = 0.0f;
    float max = 1.0f;
    
    // Parallelize the random number generation
    #pragma omp parallel
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
//*/
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
///*
// Function to perform matrix element-wise multiplication in parallel
Vector vecelemPara(const Vector *A, const Vector *B) {
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
    // Perform element-wise multiplication using OpenMP for parallelization
    #pragma omp parallel for default(none) shared(A, B, C) private(i)	
	for (i = 0; i < C.size; i++) {
        C.data[i] = A->data[i] * B->data[i];
        //printf("Thread %d processing index %d\n", omp_get_thread_num(), i);
    }

    return C;
}

// Calculate any norm of a vector using BLAS
float norm(const Vector *vec, float p)
{
    float norm_value = 0.0f;

    if (p == 1) 
    {
        // 1-norm using BLAS
        norm_value = cblas_sasum(vec->size, vec->data, 1);
    } else if (p == 2) 
    {
        // 2-norm using BLAS
        norm_value = cblas_snrm2(vec->size, vec->data, 1);
    } else if (p == INFINITY) 
    {
        // Infinity-norm (Maximum norm)
        for (int i = 0; i < vec->size; i++) 
        {
            float abs_val = fabs(vec->data[i]);
            if (abs_val > norm_value) 
            {
                norm_value = abs_val;
            }
        }
    } else 
    {
        // p-norm
        for (int i = 0; i < vec->size; i++) 
        {
            norm_value += pow(fabs(vec->data[i]), p);
        }
        norm_value = pow(norm_value, 1.0f / p);
    }

    return norm_value;
}

// Function to create a linearly spaced array
Vector createArray(float start, float end, float step)
{
    // Check if step size is valid
    if (step == 0.0f) {
        fprintf(stderr, "Error: Step size cannot be zero.\n");
        exit(EXIT_FAILURE);
    }
    if ((end - start) / step < 0) {
        fprintf(stderr, "Error: Step size is inconsistent with start and end points.\n");
        exit(EXIT_FAILURE);
    }

    // Calculate the number of elements
    int size = (int)ceilf((end - start) / step) + 1;

    // Allocate memory for the array
    Vector result;
    result.size = size;
    result.data = (float *)malloc(size * sizeof(float));
    if (result.data == NULL) {
        fprintf(stderr, "Memory allocation failed for result array.\n");
        exit(EXIT_FAILURE);
    }

    // Fill the array with values
    float current = start;
    for (int i = 0; i < size; i++) {
        result.data[i] = current;
        current += step;
        if ((step > 0 && current > end) || (step < 0 && current < end)) {
            current = end; // Ensure last value is exactly the endpoint
        }
    }

    return result;
}

// Function to perform weighted polynomial fit
Vector polyfitweighted(const Vector *x, const Vector *y, const Vector *w, int n)
{
    int rows = x->size;
    int cols = n + 1;

    // Allocate memory for the Vandermonde matrix
    Matrix V;
    V.rows = rows;
    V.cols = cols;
    V.data = (float *)calloc(V.rows * V.cols, sizeof(float));
    Vector wy;
    wy.size = rows;
    wy.data = (float *)calloc(wy.size, sizeof(float));

    if (V.data == NULL || wy.data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    // Construct the weighted Vandermonde matrix 
    for (int i = 0; i < rows; i++) { 
        // Set the last column as the weight
        V.data[i * cols + n] = w->data[i]; // Last column is the weight 
        
        // Fill the Vandermonde row starting from the higest power to the lowest
        for (int j = n - 1; j >= 0; j--) { 
            V.data[i * cols + j] = x->data[i] * V.data[i * cols + j + 1];
        }
        // Fill the Vandermonde row starting from the lowest power to the highest
        //for (int j = 0; j <= n; j++) {
        //    if (j == 0) {
        //        V.data[i * cols + j] = 1.0; // First element is 1
        //    } else {
        //        V.data[i * cols + j] = x->data[i] * V.data[i * cols + j - 1];
        //    }
        //}

        // Calculate weighted y values
        wy.data[i] = w->data[i] * y->data[i]; 
    }

    // QR decomposition of V
    int info;
    float *tau = (float *)calloc(cols, sizeof(float));
    if (tau == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    // Perform QR decomposition
    info = LAPACKE_sgeqrf(LAPACK_ROW_MAJOR, rows, cols, V.data, V.cols, tau);
    if (info != 0) {
        fprintf(stderr, "Error in QR decomposition: %d\n", info);
        free(tau);
        freeMatrix(&V);
    }

    // Allocate memory for Q and R matrices
    Matrix Q;
    Q.rows = rows;
    Q.cols = cols;
    Q.data = (float *)calloc(Q.rows * Q.cols, sizeof(float));

    // Copy the contents of V to Q
    memcpy(Q.data, V.data, rows * cols * sizeof(float));

    // Generate the orthogonal matrix Q from the QR decomposition
    info = LAPACKE_sorgqr(LAPACK_ROW_MAJOR, rows, cols, cols, Q.data, Q.cols, tau);
    if (info != 0) {
        fprintf(stderr, "Error in generating Q: %d\n", info);
        free(tau);
        freeMatrix(&V);
        freeMatrix(&Q);
    }

    // Prepare and extract the R matrix from V
    Matrix R;
    R.rows = cols;
    R.cols = cols;
    R.data = (float *)calloc(R.rows * R.cols, sizeof(float));
    if (R.data == NULL) {
        fprintf(stderr, "Memory allocation failed for R.\n");
        free(tau);
        freeMatrix(&V);
        freeMatrix(&Q);
    }

    // Copy the upper triangular part of V to R
    for (int i = 0; i < R.rows; i++) {
        for (int j = 0; j < R.cols; j++) {
            if (i <= j) {
                R.data[i * R.cols + j] = V.data[i * V.cols + j];
            } //else {  // Already handeled by calloc
              //  R.data[i * R.cols + j] = 0.0f; // Fill with zeros
            //}
        }
    }

    // Compute Q^T * (w .* y) using 'sgemv'
    // Q^T is the transpose of Q, we can use DGEMV for this
    // qwy = Q^T * wy
    Vector qwy;
    qwy.size = Q.cols;
    qwy.data = (float *)calloc(qwy.size, sizeof(float));
    float alpha = 1.0f;
    float beta  = 1.0f;
    cblas_sgemv(CblasRowMajor, CblasTrans, Q.rows, Q.cols, alpha, Q.data, 
                                    Q.cols, wy.data, 1, beta, qwy.data, 1);

    // Step 3: Solve R * p = qwy
    // p = R^(-1) * qwy
    Vector p;
    p.size = R.rows; // p should have the same size as the number of rows in R
    p.data = (float *)calloc(p.size, sizeof(float));

    // Solve R * p = qwy using cblas_dtrsv
    cblas_strsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 
                R.rows, R.data, R.cols, qwy.data, 1); // Solve R * p = qwy, store result in qwy

    // Copy the result from qwy to p
    memcpy(p.data, qwy.data, p.size * sizeof(float));

    // Free allocated memory
    freeVector(&wy);
    free(tau);
    freeMatrix(&V);
    freeMatrix(&Q);
    freeMatrix(&R);
    //freeMatrix(&QT);
    //freeVector(&Qwy);
    freeVector(&qwy);

    return p;
}

// Function to create a linearly spaced vector
Vector linspace(float start, float end, int num)
{
    Vector result;
    result.size = num;
    result.data = (float *)malloc(num * sizeof(float));

    if (result.data == NULL) 
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Calculate the step size
    float step = (end - start) / (num - 1);

    for (int i = 0; i < num; i++) 
    {
        result.data[i] = start + i * step;
    }

    return result;

}

// Function to compute the roots of a polynomial given 
// its coefficients using the Durand-Kerner method (or
// the Weierstrass method). Dis some poweful shit here
Vector roots(const Vector *coeff)
{
    int n = coeff->size - 1; // Degree of the polynomial
    Vector result;
    result.size = 2*n; // Size should be 2*n to hold both real and imaginary parts of one polynomial constant
    result.data = (float *)malloc(result.size * sizeof(float));

    if (result.data == NULL) 
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Initial guess for roots (uniformly spaced on the unit circle)
    float pi = 3.14159265358979323846f;
    float angle;
    for (int i = 0; i < n; i++) 
    {
        angle = 2 * pi * i / n; // Angle for unit circle
        result.data[i] = cos(angle); // Real part
        result.data[i + n] = sin(angle); // Imaginary part (stored in the second half of the array)
    }

    // Iterative method to refine the root estimates
    float real, imag, p, dp, denominator;
    for (int iter = 0; iter < 100; iter++) // Number of iterations
    { 
        for (int i = 0; i < n; i++) 
        {
            real = result.data[i];    
            imag = result.data[i + n];

            // Evaluate the polynomial and its derivative at the current estimate
            p  = 0.0f; // Polynomial value
            dp = 0.0f; // Derivative value

            for (int j = 0; j < coeff->size; j++) 
            {
                p += coeff->data[j] * pow(real, n - j);
                if (j < n) {
                    dp += (n - j) * coeff->data[j] * pow(real, n - j - 1);
                }
            }

            // Update the root estimate
            denominator = p / (dp + 1e-10); // Avoid division by zero
            result.data[i] -= denominator; // Update real part
            result.data[i + n] -= denominator; // Update imaginary part
        }
    }

    return result;
}

// Function to convert roots to polynomial coefficients
Vector polycoefs(const Vector *roots)
{
    int n = roots->size; // Number of roots
    Vector coeff;
    coeff.size = n + 1; // Coefficients array size is degree + 1
    coeff.data = (float *)calloc(coeff.size, sizeof(float));

    if (coeff.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Initialize coefficients: coeff[0] = 1, rest = 0
    coeff.data[0] = 1.0f; // Leading coefficient (x^n)

    // UPDATE: Below block already done with calloc() 
    //for (int i = 1; i < coeff.size; i++) {
    //    coeff.data[i] = 0.0f; // Initialize other coefficients to 0
    //}

    // Construct polynomial coefficients from roots
    float r;
    for (int i = 0; i < n; i++) {
        r = roots->data[i]; // Current root
        // Update coefficients by multiplying with (x - r)
        for (int j = i; j >= 0; j--) {
            coeff.data[j + 1] += coeff.data[j] * (-r);
        }
    }

    return coeff;
}

// Function to compute the convolution of two vectors
Vector conv(const Vector *a, const Vector *b)
{
    int n = a->size + b->size - 1; // Size of the result vector
    Vector result;
    result.size = n;
    result.data = (float *)calloc(result.size, sizeof(float));

    if (result.data == NULL) 
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // UPDATE: Below block already done with calloc() 
    //// Initialize the result vector to zero
    //for (int i = 0; i < result.size; i++) {
    //    result.data[i] = 0.0f;
    //}

    // Perform convolution
    for (int i = 0; i < a->size; i++) 
    {
        for (int j = 0; j < b->size; j++) 
        {
            result.data[i + j] += a->data[i] * b->data[j];
        }
    }

    return result;
}

//*/