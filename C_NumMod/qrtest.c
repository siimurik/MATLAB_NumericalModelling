/* gcc qrtest.c -o qr -llapacke -lblas && time ./qr */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>

typedef struct {
    float *data;
    int size;
} Vector;

typedef struct {
    float *data;
    int rows;
    int cols;
} Matrix;

void freeVector(Vector *vector);
Vector createVector(int size, float *values);
void printVector(const Vector *vector);
void freeMatrix(Matrix *matrix);
void printMatrix(const Matrix *mat);
Matrix transposeMatrix(const Matrix *mat);
Vector matvec(const Matrix *A, const Vector *x);

int main() {
    float xData[] = {-1.0f, 0.1f, 0.5f, 3.0f, 4.0f, 6.3f, 7.0f, 9.0f, 14.0f, 21.0f};
    float yData[] = {-2.0f, 0.4f, 0.7f, 2.0f, 4.0f, 3.6f, 3.8f, 6.0f, -1.0f, 12.0f};
    int sizeXY = sizeof(xData) / sizeof(xData[0]);
    Vector x = createVector(sizeXY, xData);
    Vector y = createVector(sizeXY, yData);
    Vector w;
    w.size = sizeXY;
    w.data = (float *)malloc(w.size * sizeof(float));
    for (int i = 0; i < sizeXY; i++) {
        w.data[i] = 1.0f;
    }

    int n = 3; // degree
    int len = x.size;

    // Allocate memory for the Vandermonde matrix
    Matrix V;
    V.rows = len;
    V.cols = n + 1;
    V.data = (float *)calloc(V.rows * V.cols, sizeof(float));
    Vector wy;
    wy.size = len;
    wy.data = (float *)calloc(wy.size, sizeof(float));

    if (V.data == NULL || wy.data == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    // Construct the weighted Vandermonde matrix 
    for (int i = 0; i < len; i++) { 
        // Set the last column as the weight
        V.data[i * (n + 1) + n] = w.data[i]; // Last column is the weight 
        
        // Fill the Vandermonde row starting from the higest power to the lowest
        //for (int j = n - 1; j >= 0; j--) { 
        //    V.data[i * (n + 1) + j] = x.data[i] * V.data[i * (n + 1) + j + 1];
        //}
        // Fill the Vandermonde row starting from the lowest power to the highest
        for (int j = 0; j <= n; j++) {
            if (j == 0) {
                V.data[i * (n + 1) + j] = 1.0; // First element is 1
            } else {
                V.data[i * (n + 1) + j] = x.data[i] * V.data[i * (n + 1) + j - 1];
            }
        }

        // Calculate weighted y values
        wy.data[i] = w.data[i] * y.data[i]; 
    }

    printf("V ="); printMatrix(&V);

    // QR decomposition of V
    int info;
    float *tau = (float *)calloc(n + 1, sizeof(float));
    if (tau == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    // Perform QR decomposition
    info = LAPACKE_sgeqrf(LAPACK_ROW_MAJOR, len, n + 1, V.data, V.cols, tau);
    if (info != 0) {
        fprintf(stderr, "Error in QR decomposition: %d\n", info);
        free(tau);
        freeMatrix(&V);
        return -1;
    }

    // Allocate memory for Q and R matrices
    Matrix Q;
    Q.rows = len;
    Q.cols = n + 1;
    Q.data = (float *)calloc(Q.rows * Q.cols, sizeof(float));

    // Copy the contents of V to Q
    memcpy(Q.data, V.data, len * (n + 1) * sizeof(float));

    // Generate the orthogonal matrix Q from the QR decomposition
    info = LAPACKE_sorgqr(LAPACK_ROW_MAJOR, len, n + 1, n + 1, Q.data, Q.cols, tau);
    if (info != 0) {
        fprintf(stderr, "Error in generating Q: %d\n", info);
        free(tau);
        freeMatrix(&V);
        freeMatrix(&Q);
        return -1;
    }

    // Prepare and extract the R matrix from V
    Matrix R;
    R.rows = n + 1;
    R.cols = n + 1;
    R.data = (float *)calloc(R.rows * R.cols, sizeof(float));
    if (R.data == NULL) {
        fprintf(stderr, "Memory allocation failed for R.\n");
        free(tau);
        freeMatrix(&V);
        freeMatrix(&Q);
        return -1;
    }

    // Copy the upper triangular part of V to R
    for (int i = 0; i < R.rows; i++) {
        for (int j = 0; j < R.cols; j++) {
            if (i <= j) {
                R.data[i * R.cols + j] = V.data[i * V.cols + j];
            } else {
                R.data[i * R.cols + j] = 0.0f; // Fill with zeros
            }
        }
    }

    // Print the results
    printf("\nQ = "); printMatrix(&Q);
    printf("\nR = "); printMatrix(&R);
    printf("\nw.*y = ");
    printVector(&wy);

    // Method 1: Step-by-step
    // Transpose  matrix Q
    Matrix QT = transposeMatrix(&Q);
    printf("\nQ^T = ");
    printMatrix(&QT);
    // Perform matrix-vector multiplcation
    Vector Qwy = matvec(&QT, &wy);
    printf("\nQ^T * (w .* y) = ");
    printVector(&Qwy);

    // Method 2: Compute Q^T * (w .* y) using 'sgemv'
    // Q' is the transpose of Q, we can use DGEMV for this
    // b = Q' * wy
    Vector qwy;
    qwy.size = Q.cols;
    qwy.data = (float *)calloc(qwy.size, sizeof(float));
    float alpha = 1.0f;
    float beta  = 1.0f;
    cblas_sgemv(CblasRowMajor, CblasTrans, Q.rows, Q.cols, alpha, Q.data, 
                                    Q.cols, wy.data, 1, beta, qwy.data, 1);
    printf("\nQ^T * (w .* y) = ");
    printVector(&qwy);

    // Step 3: Solve R * p = b
    // p = R^(-1) * b
    Vector p;
    p.size = R.rows; // p should have the same size as the number of rows in R
    p.data = (float *)calloc(p.size, sizeof(float));

    // Solve R * p = qwy using cblas_dtrsv
    cblas_strsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 
                R.rows, R.data, R.cols, qwy.data, 1); // Solve R * p = qwy, store result in qwy

    // Copy the result from qwy to p
    memcpy(p.data, qwy.data, p.size * sizeof(float));

    printf("\np = ");
    printVector(&p);

    // Free allocated memory
    freeVector(&x);
    freeVector(&y);
    freeVector(&w);
    freeVector(&wy);
    free(tau);
    freeMatrix(&V);
    freeMatrix(&Q);
    freeMatrix(&R);
    freeMatrix(&QT);
    freeVector(&Qwy);
    freeVector(&qwy);
    freeVector(&p);

    return 0;

}

// Free allocated vector memory
void freeVector(Vector *vector) 
{
    free(vector->data);
    vector->data = NULL; // Avoid dangling pointer
}

// Function to save values into a vector
Vector createVector(int size, float *values)
{
    Vector vector;
    vector.size = size;
    vector.data = (float *)malloc(size * sizeof(float)); // Allocate memory

    if (vector.data == NULL) 
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Populate the vector with the provided values 
    memcpy(vector.data, values, vector.size * sizeof(float));

    return vector;
}

void freeMatrix(Matrix *matrix) 
{
    free(matrix->data);
    matrix->data = NULL;    // Avoid dangling pointer.
}  

// Function for printing a matrix
void printMatrix(const Matrix *mat) 
{
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

void printVector(const Vector *vector) 
{
    printf("[");
    
    if (vector->size > 10) 
    {
        // Print the first three elements
        for (int i = 0; i < 3; i++) 
        {
            printf("%.4f", vector->data[i]);
            if (i < 2) {
                printf(", ");
            }
        }
        // Print ellipsis
        printf(", ...");
        // Print the last three elements
        for (int i = vector->size - 3; i < vector->size; i++) 
        {
            printf(", %.4f", vector->data[i]);
        }
    } else 
    {
        // Print all elements if size is 10 or less
        for (int i = 0; i < vector->size; i++) 
        {
            printf("%.4f", vector->data[i]);
            if (i < vector->size - 1) 
            {
                printf(", ");
            }
        }
    }
    
    printf("]\n");
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

// Function to perform matrix-vector multiplication
Vector matvec(const Matrix *A, const Vector *x) 
{
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