/* Compile and execute with
    $ gcc yl7.c -o yl7 -lm -lblas -llapacke
    $ ./yl7
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>

#define pi 3.14159265358979323846

typedef struct {
    float *data;
    int size;
} Vector;

typedef struct {
    float *data;
    int rows;
    int cols;
} Matrix;

float  polyval(Vector *coefficients, float x);
Vector polyvalVec(const Vector *coeff, const Vector *points);
void   freeVector(Vector *vector);
Vector createVector(int size, float *values);
void   printVector(const Vector *vector);
Vector linspace(float start, float end, int num);
Vector roots(const Vector *coeff);
Vector polycoefs(const Vector *roots);
Vector conv(const Vector *a, const Vector *b);
void   freeMatrix(Matrix *matrix);
void   printMatrix(const Matrix *mat);
Matrix transposeMatrix(const Matrix *mat);
Vector matvec(const Matrix *A, const Vector *x);
Vector polyfitweighted(const Vector *x, const Vector *y, const Vector *w, int n);

int main() 
{
    printf("Exercise 1:\n");
    // Define polynomial coefficients for p(x) = 4x^5 - 2x^4 + x^3 + 7^x - 9.
    float coefficients[] = {4.0f, -2.0f, 1.0f, 0.0f, 7.0f, -9.0f}; // Coefficients for x^2, x^1, x^0
    int coeff_size = sizeof(coefficients) / sizeof(coefficients[0]);
    Vector cVec = createVector(coeff_size, coefficients);
    printf("coefs: "); printVector(&cVec);

    float vals[2] = {0.0f, 6.0f};
    float p;
    for (int i = 0; i < 2; i++)
    {
        p = polyval(&cVec, vals[i]);
        printf("p(%2.1f) = %.4f\n", vals[i], p);
    }
   
    Vector xVec = linspace(-1.0f, 8.0f, 100);
    Vector pVec = polyvalVec(&cVec, &xVec);
    printf("p(x) = "); printVector(&pVec);

    printf("\nExercise 2:\n");
    float rVal[3] = {1.0f, -1.0f, -12.0f};
    Vector rVec = createVector(3, rVal);
    Vector r = roots(&rVec);
    printf("r = "); printVector(&r);

    printf("\nExercise 3:\n");
    float rootVals[4] = {-1.0f, 3.0f, 2.0f, 5.0f};
    Vector rootVec = createVector(4, rootVals);
    Vector coefVec = polycoefs(&rootVec);
    printf("Coefficients of the polynomial:\n"); 
    printVector(&coefVec);

    printf("\nExercise 4:\n");
    Vector p1, p2;
    p1.size = 3;
    p1.data = (float *)malloc(p1.size * sizeof(float));
    p1.data[0] =  1.0f;
    p1.data[1] = -2.0f;
    p1.data[2] =  1.0f;
    printf("p1 = "); printVector(&p1);

    p2.size = 2;
    p2.data = (float *)malloc(p2.size * sizeof(float));
    p2.data[0] = 1.0f;
    p2.data[1] = 2.0f;
    printf("p2 = "); printVector(&p2);

    Vector pConv = conv(&p1, &p2);
    printf("p1 * p2 = conv(p1, p2) = ");
    printVector(&pConv);

    printf("\nExercise 5:\n");
    float xData[] = {-1.0f, 0.1f, 0.5f, 3.0f, 4.0f, 6.3f, 7.0f, 9.0f, 14.0f, 21.0f};
    float yData[] = {-2.0f, 0.4f, 0.7f, 2.0f, 4.0f, 3.6f, 3.8f, 6.0f, -1.0f, 12.0f};
    int sizeXY = sizeof(xData)/sizeof(xData[0]);
    Vector xDataVec = createVector(sizeXY, xData);
    Vector yDataVec = createVector(sizeXY, yData);
    Vector k;
    k.size = sizeXY;
    k.data = (float *)malloc(k.size * sizeof(float));
    for (int i = 0; i < sizeXY; i++){
        k.data[i] = 1.0f;
    }
    //printf("k = "); printVector(&k);
    
    int degree = 3;
    Vector pPoints = polyfitweighted(&xDataVec, &yDataVec, &k, degree);
    printf("pPoints = "); printVector(&pPoints);


    // Free allocated memory
    freeVector(&cVec);
    freeVector(&xVec);
    freeVector(&pVec);
    freeVector(&rVec);
    freeVector(&r);
    freeVector(&rootVec);
    freeVector(&coefVec);
    freeVector(&p1);
    freeVector(&p2);
    freeVector(&pConv);
    freeVector(&xDataVec);
    freeVector(&yDataVec);
    freeVector(&k);
    freeVector(&pPoints);

    return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Function to evaluate polynomial at a given point x
float polyval(Vector *coefficients, float x) 
{
    float result = 0.0f;
    for (int i = 0; i < coefficients->size; i++) 
    {
        result += coefficients->data[i] * pow(x, coefficients->size - 1 - i);
    }
    return result;
}

// Function to evaluate a polynomial at given points
Vector polyvalVec(const Vector *coeff, const Vector *points) 
{
    Vector result;
    result.size = points->size;
    result.data = (float *)malloc(result.size * sizeof(float));

    if (result.data == NULL) 
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 0; j < points->size; j++) 
    {
        float x = points->data[j];
        float value = 0.0f;

        // Evaluate polynomial using Horner's method
        for (int i = 0; i < coeff->size; i++) 
        {
            value = value * x + coeff->data[i];
        }

        result.data[j] = value;
    }

    return result;
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
    for (int i = 0; i < n; i++) 
    {
        float angle = 2 * pi * i / n; // Angle for unit circle
        result.data[i] = cos(angle); // Real part
        result.data[i + n] = sin(angle); // Imaginary part (stored in the second half of the array)
    }

    // Iterative method to refine the root estimates
    for (int iter = 0; iter < 100; iter++) // Number of iterations
    { 
        for (int i = 0; i < n; i++) 
        {
            float real = result.data[i];        // declaring varibles inside a for loop is crazy
            float imag = result.data[i + n];    // this comes from a Fortran developer

            // Evaluate the polynomial and its derivative at the current estimate
            float  p = 0.0f; // Polynomial value
            float dp = 0.0f; // Derivative value

            for (int j = 0; j < coeff->size; j++) 
            {
                p += coeff->data[j] * pow(real, n - j);
                if (j < n) {
                    dp += (n - j) * coeff->data[j] * pow(real, n - j - 1);
                }
            }

            // Update the root estimate
            float denominator = p / (dp + 1e-10); // Avoid division by zero
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
    for (int i = 0; i < n; i++) {
        float r = roots->data[i]; // Current root
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

// Function to perform weighted polynomial fit
Vector polyfitweighted(const Vector *x, const Vector *y, const Vector *w, int n) {
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
        //for (int j = n - 1; j >= 0; j--) { 
        //    V.data[i * cols + j] = x.data[i] * V.data[i * cols + j + 1];
        //}
        // Fill the Vandermonde row starting from the lowest power to the highest
        for (int j = 0; j <= n; j++) {
            if (j == 0) {
                V.data[i * cols + j] = 1.0; // First element is 1
            } else {
                V.data[i * cols + j] = x->data[i] * V.data[i * cols + j - 1];
            }
        }

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

    // Method 1: Step-by-step. Slower
    // Transpose  matrix Q
    /*
    Matrix QT = transposeMatrix(&Q);
    printf("Matrix Q^T:\n");
    printMatrix(&QT);
    // Perform matrix-vector multiplcation
    Vector Qwy = matvec(&QT, &wy);
    printf("Vector Q *w.*y:\n");
    printVector(&Qwy);
    */

    // Method 2: Compute Q^T * (w .* y) using 'sgemv'
    // Q' is the transpose of Q, we can use DGEMV for this
    // qwy = Q' * wy
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

