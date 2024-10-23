// Compile and execute with
//     $ gcc yl1.c -o yl1 -lm -llapacke -lcblas
//     $ ./yl1
//-------------------------------------------------------     

#include <stdio.h>      // for `printf` to work
#include <math.h>       // mathematical functions (exp, cos, log, ...)
#include <stdlib.h>     // For EXIT_FAILURE to work
#include <lapacke.h>    // Linear Algebra PACKage in Fortran
#include <cblas.h>      // Basic Linear Algebra Subprograms in C    
#include <string.h>     // For the memcpy() command
//#include <omp.h>        // OpenMP for parallelization

typedef struct {
    float *data;
    int rows;
    int cols;
} Matrix;

typedef struct {
    float *data;
    int size;
} Vector;

// Functions for Exercise 2.
float u_func(float x, float y, float z)
{
    return (x*y + z/(x-4.0*y*y) + z*z )/(x*x*x + y*y - 3.0*z*x);
}

float v_func(float x, float y, float z){
    return pow(x*y*z + x*x - (2.0*y)*(2.0*y), 5.0/(x+y));
}

// Function for exercise 3.
float y_func(float x)
{
    return (x*x-3.0)*pow(2.0+x,4.0) - 5.0*exp(x) + 2.0*cos(x+1.0);
}

// Function for exercise 4.
float z_func(float u)
{
    return pow(u,5.0) - 3.0*pow(u,4.0) + pow(u,3.0) + 1.0;
}

// Function for exercise 5.
float F(float x, float y){
    return sin(x-y)/(x*x) + pow( cos(2.0*x+y)/pow((x-y),4.0), 1.0/3.0 );
}

//--------------------------------------------------------------------------------
// Function to free the allocated memory
void freeMatrix(Matrix *matrix) {
    free(matrix->data);
    matrix->data = NULL; // Avoid dangling pointer
}

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

// Function to print the matrix
void printM(const Matrix *matrix) {
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            printf("%.2f ", matrix->data[i * matrix->cols + j]);
        }
        printf("\n");
    }
}

// If the function is declared as 
//         void printm(Matrix matvec){...}
// then there are a few things to consider. 
//
// Firstly, this implementation makes a copy of the entire 
// structure, which can be inefficient for large matrices.
//
// Secondly, not using `const`, means the function could 
// potentially modify the matvec structure, even though it 
// doesn't.
//
// Thirdly, both implementations access the matrix data in
// the same way, using the formula data[i * cols + j]. However, 
// the second implementation uses the arrow operator (->) to 
// access members of the structure through a pointer, which 
// is appropriate.
//
// Finally, the current method can be called by
//         printm(matrix_A);
// but, when using pointers, the approach should be
//         printm(&matrix_A);
// 
// Therefore we declare this funtion in the following matter:
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
    // "free(): double free detected in tcache 2
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

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, C.rows, C.cols, A->cols, 
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
    C.data = (float *)malloc(C.rows * C.cols * sizeof(float));
    memcpy(C.data, B->data, C.rows * C.cols * sizeof(float));

    cblas_saxpy(C.rows * C.cols, 1.0, A->data, 1, C.data, 1);

    return C;
}

// Matrix scalar multiplication function
Matrix matscal(const Matrix *A, float scalar) {
    Matrix C;
    C.rows = A->rows;
    C.cols = A->cols;
    C.data = (float *)malloc(C.rows * C.cols * sizeof(float));
    memcpy(C.data, A->data, C.rows * C.cols * sizeof(float));

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
    float determinant = 1.0;

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

    free(ipiv);
    freeMatrix(&matCopy); // Temporary. Not necessary outside the function.
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

    // Perform element-wise multiplication using OpenMP for parallelization in necessary
    //#pragma omp parallel for
    for (int i = 0; i < C.rows; i++) {
        for (int j = 0; j < C.cols; j++) {
            C.data[i * C.cols + j] = A->data[i * A->cols + j] * B->data[i * B->cols + j];
        }
    }

    return C;
}


//---------------------------------------------------------------------------
// Functions for vectors

// Free allocated vector memory
void freeVector(Vector *vector) {
    free(vector->data);
    vector->data = NULL; // Avoid dangling pointer
}

// Function to save values into a vector
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
    
    // Paralellize if necessary
    //#pragma omp parallel for
    //for (int i = 0; i < size; i++) {
    //    vector.data[i] = values[i];
    //}

    return vector;
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


//-------------------------------------------------------

int main()
{
    float x, y, z, w, u, v;

    printf("Exercise 1.\n");

    x = 13.0;
    y =  7.0;
    z = sqrt(x);
    w = 6.0*x*y - x/y;

    printf("z = %.4f\nw = %.4f\n", z, w);

    printf("\nExercise 2.\n");
    x =  2.0;
    y = -1.0;
    z =  5.0;
    u = u_func(x, y, z);
    v = v_func(x, y, z);
    printf("x = %.1f, y = %.1f, z = %.1f\n", x, y, z);
    printf("u = %.4f\nv = %.4f\n\n", u, v);

    x = -6.0;
    y =  5.0;
    z = -2.0;
    u = u_func(x, y, z);
    v = v_func(x, y, z);
    printf("x = %.1f, y = %.1f, z = %.1f\n", x, y, z);
    printf("u = %.4f\nv = %.4e\n", u, v);
    
    printf("\nExercise 3.\n");
    float x_vec[3] = {-1.0, 4.0, 3.0};
    int i;

    for (i = 0; i < 3; i++)
    {
        printf("y(%.0f) = %.4e\n", x_vec[i], y_func(x_vec[i]) );
    }

    printf("\nExercise 4.\n");
    float u_vec[2] = {2.0, 8.0};
    for (i = 0; i < 2; i++)
    {
        printf("z(%.0f) = %.4f\n", u_vec[i], z_func(u_vec[i]) );
    }

    printf("\nExercise 5.\n");
    printf("F(50, -30) = %.4f\n", F(50.0, -30.0));
    

    printf("\nExercise 6.\n");
    int rows, cols;
    rows = 3;
    cols = 3;

    // Declare and print values for the A matrix
    float mat_values[9] = {3, 12, 52, 4, 6, -11, -2, 7, 2};
    Matrix A = createMatrix(rows, cols, mat_values);
    printf("Matrix A:\n");
    //printM(&A);
    printMatrix(&A);

    // Declare and print values for the b vector
    float vec_values[3] = {13.0, -2.0, 5.0};
    int size = sizeof(vec_values) / sizeof(vec_values[0]); // Calculate the size of the vector
    //printf("size = %d\n", size);
    Vector b = createVector(size, vec_values);
    printf("Vector b: ");
    printVector(&b);

    // Find the transposed and inverse matrix of matrix A
    Matrix A_T   = transposeMatrix(&A);
    Matrix inv_A = inverseMatrix(&A);

    
    // B = 2*A^T + A^-1
    Matrix twoTimesAT = matscal(&A_T, 2.0);
    Matrix B = matadd( &twoTimesAT, &inv_A );
    printf("\nMatrix B: \n");
    printMatrix(&B);

    // C = A*B
    Matrix C = matmul(&A, &B);
    printf("\nMatrix C: \n");
    printMatrix(&C);
    
    // det(A)
    printf("\nDeterminant of matrix A: ");
    float det_A = det(&A);
    printf("%.4f.\n\n", det_A);
    
    Vector d = matvec(&A, &b);
    printf("Vector d: \n");
    printVector(&d);

    // D = A.*B
    Matrix D = matelem(&A, &B);
    printf("\nMatrix D:\n");
    printMatrix(&D);

    // Free the allocated memory
    freeMatrix(&A);
    freeMatrix(&B);
    freeMatrix(&C);
    freeMatrix(&D);
    freeMatrix(&A_T);
    freeMatrix(&inv_A);
    freeMatrix(&twoTimesAT);
    freeVector(&b);
    freeVector(&d);

    return 0;
}
