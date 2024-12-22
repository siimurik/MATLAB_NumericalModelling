/* Compile and execute with:
    $ gcc -c splash.c -o splash.o
    $ ar rcs libsplash.a splash.o
    $ gcc yl9_10.c -o yl910 -L. -lsplash -llapacke -lblas -lgsl -lm -fopenmp
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include <omp.h>
#include "splash.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

Matrix vander(const Vector *x);
Vector newlinsolve(const Matrix *A, const Vector *b);
Vector interp1(const Vector *s, const Vector *t, const Vector *ss, const char *method);
Vector createArray(float start, float end, float step);
Vector polyfitweighted(const Vector *x, const Vector *y, const Vector *w, int degree);
Matrix zerosMatrix(int rows, int cols);
Vector zerosVector(int size);
float sum(Vector *vec);
void pow_vector(Vector *vec, float power);

int main()
{
    float datavec[] = {1.0f, 2.0f, 3.0f}; // 
    int vander_dim = sizeof(datavec) / sizeof(datavec[0]);
    Vector data = createVector(vander_dim, datavec);
    //printVector(&data);

    Matrix V = vander(&data);
    printf("V ="); printMatrix(&V);
    //Matrix Vinv = inverseMatrix(&V);

    float BvecData[] = {5.0f, 1.0f, -1.0f};
    Vector B = createVector(vander_dim, BvecData);
    printf("B ="); printVector(&B);

    //Vector X = matvec(&Vinv, &B);
    Vector X = linsolve(&V, &B);
    printf("X ="); printVector(&X);


    printf("\nExercise 4.\n");
    float s_data[] = {0.0f, 0.25f,  0.5f, 0.75f, 1.0f, 1.25f};
    float t_data[] = {0.0f, 25.0f, 49.4f, 73.0f, 96.4f, 119.4f};
    int s_size = sizeof(s_data) / sizeof(s_data[0]);
    Vector s = createVector(s_size, s_data);
    Vector t = createVector(s_size, t_data);
    printf("s = "); printVector(&s);
    printf("t = "); printVector(&t);


    Vector S1_0 = interp1(&s, &t, &s, "linear");
    printf("S1_0 = "); printVector(&S1_0);

    Vector ss = createArray(0.0f, 1.25f, 0.001f);
    printf("ss = "); printVector(&ss);

    Vector S3_2 = interp1(&s, &t, &ss, "spline");
    printf("S3_2 = "); printVector(&S3_2);
    
    float s32valdata[] = {0.45f};
    Vector s32val = createVector(1, s32valdata);
    Vector S3_2_val = interp1(&s, &t, &s32val, "linear");
    printf("S3_2(%.2f)= ", s32valdata[0]); printVector(&S3_2_val);

//fprintf("\nÜlesanne 7\n")
//clear
//x = [1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2];
//y = [9.08 10.43 11.9 13.48 15.19 17.03 19.01 21.13 23.39];
//k = [1 1 2 5 1 4 2 2 1];

    printf("\nExercise 7.\n");
    float xData[] = {1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f};
    float yData[] = {9.08f, 10.43f, 11.9f, 13.48f, 15.19f, 17.03f, 19.01f, 21.13f, 23.39f};
    float kData[] = {1.0f, 1.0f, 2.0f, 5.0f, 1.0f, 4.0f, 2.0f, 2.0f, 1.0f};
    //float kData[] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    int ndim = sizeof(xData) / sizeof(xData[0]);
    Vector x = createVector(ndim, xData);
    Vector y = createVector(ndim, yData);
    Vector k = createVector(ndim, kData);

    Vector p = polyfitweighted(&x, &y, &k, 3);
    printf("p = "); printVector(&p);

    // TODO: Finish alternative method
    //printf("\nAlternative method\n");
    //int N = 3;
    //int M = N + 1;
    //
    //Matrix A = zerosMatrix(M, M);
    //Vector b = zerosVector(M);
    //Vector c = zerosVector(M);
    //
    //// Compute B and A
    //for (int i = 0; i < M; i++) {
    //    pow_vector(&x, M - i);
    //    pow_vector(&y, M - i);
    //    b.data[i] = sum((Vector){ k.data, k.size }) * sum((Vector){ y.data, y.size }) * sum((Vector){ c.data, c.size });
    //    for (int j = 0; j < M; j++) {
    //        pow_vector(&x, 2 * M - i - j);
    //        A.data[i * M + j] = sum((Vector){ k.data, k.size }) * sum((Vector){ c.data, c.size });
    //    }
    //}
    //printMatrix(&A);
    //printVector(&b);



    freeVector(&data);
    freeMatrix(&V);
    //freeMatrix(&Vinv);
    freeVector(&B);
    freeVector(&X);
    freeVector(&s);
    freeVector(&t);
    freeVector(&ss);
    freeVector(&S1_0);
    freeVector(&S3_2);
    freeVector(&s32val);
    freeVector(&x);
    freeVector(&y);
    freeVector(&k);
    freeVector(&p);
    //freeMatrix(&A);
    //freeVector(&b);
    //freeVector(&c);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////////

Matrix vander(const Vector *x) {
    int rows = x->size;
    int cols = rows; // The number of columns should be equal to the number of rows

    // Allocate memory for the Vandermonde matrix
    Matrix V;
    V.rows = rows;
    V.cols = cols;
    V.data = (float *)calloc(V.rows * V.cols, sizeof(float));

    // Construct the Vandermonde matrix
    for (int i = 0; i < rows; i++) {
        // Fill the Vandermonde row starting from the highest power to the lowest
        for (int j = cols - 1; j >= 0; j--) {
            if (j == cols - 1) {
                V.data[i * cols + j] = 1.0; // Last element is 1 (lowest power)
            } else {
                V.data[i * cols + j] = x->data[i] * V.data[i * cols + j + 1];
            }
        }
    }

    return V;
}

//Matrix vander(const Vector *x) {
//    int rows = x->size;
//    int cols = rows; // The number of columns should be equal to the number of rows
//
//    // Allocate memory for the Vandermonde matrix
//    Matrix V;
//    V.rows = rows;
//    V.cols = cols;
//    V.data = (float *)calloc(V.rows * V.cols, sizeof(float));
//
//    // Construct the Vandermonde matrix
//    for (int i = 0; i < rows; i++) {
//        // Fill the Vandermonde row starting from the lowest power to the highest
//        for (int j = 0; j < cols; j++) {
//            if (j == 0) {
//                V.data[i * cols + j] = 1.0; // First element is 1
//            } else {
//                V.data[i * cols + j] = x->data[i] * V.data[i * cols + j - 1];
//            }
//        }
//    }
//
//    return V;
//}



// Function to solve the linear system Ax = b
Vector newlinsolve(const Matrix *A, const Vector *b) {
    // Check for dimension compatibility
    if (A->cols != b->size) {
        fprintf(stderr, "Error: Incompatible dimensions for system.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the result vector
    Vector x;
    x.size = A->cols;
    x.data = (float *)malloc(x.size * sizeof(float));

    if (x.data == NULL) {
        fprintf(stderr, "Memory allocation failed for result vector.\n");
        exit(EXIT_FAILURE);
    }

    // Copy the input vector b to the result vector x
    memcpy(x.data, b->data, x.size * sizeof(float));

    // Allocate memory for the pivot indices
    int *ipiv = (int *)malloc(A->rows * sizeof(int));
    if (ipiv == NULL) {
        fprintf(stderr, "Memory allocation failed for pivot indices.\n");
        free(x.data);
        exit(EXIT_FAILURE);
    }

    // LAPACKE_sgesv requires the matrix data in column-major order.
    // LAPACKE_sgesv modifies the input matrix, so we need a copy of A's data.
    float *A_copy = (float *)malloc(A->rows * A->cols * sizeof(float));
    if (A_copy == NULL) {
        fprintf(stderr, "Memory allocation failed for matrix copy.\n");
        free(x.data);
        free(ipiv);
        exit(EXIT_FAILURE);
    }
    memcpy(A_copy, A->data, A->rows * A->cols * sizeof(float));

    // Solve the system using LAPACKE_sgesv
    int info = LAPACKE_sgesv(LAPACK_ROW_MAJOR, A->rows, 1, A_copy, A->cols, ipiv, x.data, 1);
    if (info != 0) {
        fprintf(stderr, "Error: LAPACKE_sgesv failed with info = %d.\n", info);
        free(x.data);
        free(ipiv);
        free(A_copy);
        exit(EXIT_FAILURE);
    }

    // Free allocated memory for pivot indices and matrix copy
    free(ipiv);
    free(A_copy);

    return x; // Return the result vector
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
    result.data = (float *)malloc(result.size * sizeof(float));
    if (result.data == NULL) {
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
        free(result.data);
        exit(EXIT_FAILURE);
    }

    // Convert float arrays to double for GSL compatibility
    double *s_double = (double *)malloc(s->size * sizeof(double));
    double *t_double = (double *)malloc(t->size * sizeof(double));
    double *ss_double = (double *)malloc(ss->size * sizeof(double));
    if (!s_double || !t_double || !ss_double) {
        fprintf(stderr, "Error: Memory allocation failed for temporary arrays.\n");
        free(result.data);
        free(s_double);
        free(t_double);
        free(ss_double);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < s->size; i++) {
        s_double[i] = s->data[i];
        t_double[i] = t->data[i];
    }
    for (int i = 0; i < ss->size; i++) {
        ss_double[i] = ss->data[i];
    }

    // Create a GSL interpolation object and accelerator
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(interp_type, s->size);

    // Initialize the interpolation object
    if (gsl_spline_init(spline, s_double, t_double, s->size) != GSL_SUCCESS) {
        fprintf(stderr, "Error: Failed to initialize interpolation.\n");
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
        free(result.data);
        free(s_double);
        free(t_double);
        free(ss_double);
        exit(EXIT_FAILURE);
    }

    // Perform interpolation for each value in ss
    for (int i = 0; i < ss->size; i++) {
        if (ss->data[i] < s->data[0] || ss->data[i] > s->data[s->size - 1]) {
            fprintf(stderr, "Error: Interpolation point ss[%d] = %f is out of bounds.\n", i, ss->data[i]);
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            free(result.data);
            free(s_double);
            free(t_double);
            free(ss_double);
            exit(EXIT_FAILURE);
        }
        result.data[i] = (float)gsl_spline_eval(spline, ss_double[i], acc);
    }

    // Free GSL resources and temporary arrays
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    free(s_double);
    free(t_double);
    free(ss_double);

    return result; // Return the interpolated vector
}

// Function to initialize a Matrix with zeros
Matrix zerosMatrix(int rows, int cols) 
{
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (float *)calloc(rows * cols, sizeof(float));
    return mat;
}

// Function to initialize a Matrix with zeros
Vector zerosVector(int size) 
{
    Vector vec;
    vec.size = size;
    vec.data = (float *)calloc(size, sizeof(float));
    return vec;
}

// Function to sum the elements of a Vector
float sum(Vector *vec) 
{
    float result = 0.0;
    for (int i = 0; i < vec->size; i++) {
        result += vec->data[i];
    }
    return result;
}

// Function to compute the power of each element of a Vector
void pow_vector(Vector *vec, float power) {
    for (int i = 0; i < vec->size; i++) {
        vec->data[i] = pow(vec->data[i], power);
    }
}
