/* Compile and execute with:
    $ gcc -c dashpack.c -o dashpack.o
    $ ar rcs libdashpack.a dashpack.o
    $ gcc yl8.c -o yl8 -L. -ldashpack -llapacke -lblas -lm -fopenmp
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include <omp.h> 
#include "dashpack.h"

#define DIM 2   // x, y
#define pi 3.14159265358979323846

void    sim(double (*g)(double), double *x_init, int *count, double tol, int max_iter);
double  f(double x);
double  g(double x);
double  y(double x);
Vector  y_vector(const Vector *x);
void    F(double *z, double *result);
void    dF(double *z, double J[DIM][DIM]);
void    newtonVec(void (*F)(double*, double*), void (*dF)(double*, double[DIM][DIM]), 
            double **initial_values, int initialValuesSize, double epsilon, int max_iter);

int main()
{
    printf("Exercise 1.\n");
    double x_init = 1.1;
    double tol = 1.0E-6;
    int max_iter = 1000;
    int count = 0;
    sim(&g, &x_init, &count, tol, max_iter);
    printf("SIM: x = %f, Iterations: %d\n", x_init, count);

    printf("\nExercise 2.\n");
    double cData[] = {-4.0, 3.0, 6.0, 1.0, 0.0, -7.0};
    int size = sizeof(cData)/sizeof(cData[0]);
    Vector c = createVector(size, cData);
    printf("c = "); printVector(&c);
    double c2   = norm(&c, 2);
    double c1   = norm(&c, 1);
    double cinf = norm(&c, INFINITY);
    printf("||c||_2   = %.4f\n", c2  );
    printf("||c||_1   = %.4f\n", c1  );
    printf("||c||_inf = %.4f\n", cinf);

    printf("\nExercise 3.\n");
    double Adata[] = { -5.0,  3.0,
                        2.0, -1.0};
    Matrix A = createMatrix(2, 2, Adata);
    printf("A = "); printMatrix(&A);
    double Bdata[] = {  4.0, -2.0,
                        7.0,  0.0};
    Matrix B = createMatrix(2, 2, Bdata);
    printf("B = "); printMatrix(&B);
    printf("a) Element-wise multiplication:\n");
    Matrix C = matelem(&A, &B);
    printf("C = "); printMatrix(&C);

    printf("b) First element of the second row:\n");
    printf("B[1][0] = %f\n", B.data[1 * B.cols + 0]);

    printf("\nExercise 4.\n");
    printf("y(5.0) = %lf\n", y(5.0));
    double xData[4] = {0.0, 1.0, 4.0, 3.0};
    Vector y_vec;
    y_vec.size = 4;
    y_vec.data = (double *)malloc(y_vec.size * sizeof(double));
    for (int i = 0; i < 4; i++)
    {
        y_vec.data[i] = y(xData[i]);
    }
    //int xSize = sizeof(xData)/sizeof(xData[0]);
    //Vector x = createVector(xSize, xData);
    //Vector y_vec = y_vector(&x);
    printf("y_vec = "); printVector(&y_vec);

    printf("\nExercise 5.\n");
    tol               = 1.0E-6;
    max_iter          = 10;
    double xy[DIM]    = {-1.0, 2.5};
    double *xy_vec[1] = {xy};
    int xyVecSize = sizeof(xy_vec)/sizeof(xy_vec[0]);
    printf("Newton's Method:\n");
    newtonVec(F, dF, xy_vec, xyVecSize, tol, max_iter);

    freeMatrix(&A);
    freeMatrix(&B);
    freeMatrix(&C);
    freeVector(&c);
    //freeVector(&x);
    freeVector(&y_vec);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////

// Simple iteraion method
void sim(double (*g)(double), double *x_init, int *count, double tol, int max_iter ) 
{

    // Initialize variables
    double x_new = *x_init;
    double x_old = *x_init - 1.0; // Initial value different from x_new
    *count = 0;

    // Iteration process
    while (fabs(x_new - x_old) >= tol && *count < max_iter) {
        x_old = x_new;
        x_new = g(x_old);
        (*count)++;
    }
    *x_init = x_new;

    // Check if maximum iteration count was exceeded
    if (*count >= max_iter) {
        printf("Warning: Maximum iteration count reached.\n");
    }

}

double f(double x)
{
    return pow(x,5) - 2.0*x*x*x + 2.0*x - 2.0;
}

double g(double x)
{   // Extremly IMPORTANT: if you write the power of the 
    // root to be 1/5 instead of 1.0/5.0, it WON'T WORK.
    return pow(2.0*x*x*x - 2.0*x + 2.0, 1.0/5.0);
}

double y(double x)
{
    return x*x*x - 2.0*x*x + 7.0;
}

Vector y_vector(const Vector *x)
{
    Vector y;
    y.size = x->size;
    y.data = (double *)malloc(y.size * sizeof(double));

    if (y.size != x->size) 
    {
        fprintf(stderr, "Error: Input and output vectors must have the same size.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < x->size; i++) 
    {
        y.data[i] = x->data[i] * x->data[i] * x->data[i] 
                   - 2.0 * x->data[i] * x->data[i] 
                   + 7.0;
    }
    return y;
}

// System of nonlinear equations for x, y and z (z0, z1, z2)
void F(double *z, double *result)
{
    result[0] = z[0]*z[0] - 2.0*z[1] + 4.0;
    result[1] = z[1]*z[1] - z[0]*z[0]*z[0]*z[0] + 4.0*z[0] - 1.0;
}

// Define the Jacobian dF(z)
void dF(double *z, double J[DIM][DIM])
{
    J[0][0] = 2.0*z[0];                 J[0][1] = -2.0    ; 
    J[1][0] = 4.0 - 4.0*z[0]*z[0]*z[0]; J[1][1] = 2.0*z[1];
}

// Newton's iteration method for vector inputs
void newtonVec(void (*F)(double*, double*), void (*dF)(double*, double[DIM][DIM]), 
                   double **initial_values, int initialValuesSize, double epsilon, int max_iter) 
{
    int count;
    double z[DIM];
    double delta_z[DIM];
    Vector Fz; 
    Fz.size = DIM;
    Fz.data = (double *)malloc(Fz.size * sizeof(double));
    // Set up the Jacobi matrix J
    double J[DIM][DIM]; // in LAPACK/BLAS needs a pointer
    int ipiv[DIM]; int info; // For inverting J using LAPACK
    int i, j, k, l;

    for (i = 0; i < initialValuesSize; i++) 
    {
        for (j = 0; j < DIM; j++)
        {
            z[j] = initial_values[i][j];
        }
        count = 0;
        
        while (1) 
        { // While TRUE, aka go on forever
            F(z, Fz.data);
            //printf("Iteration %d: z = [%f, %f], F(z) = [%f, %f]\n", count, z[0], z[1], Fz.data[0], Fz.data[1]);
            // The main breaking conditions that stop the code from running forever
            if (norm(&Fz, 1) < epsilon || count >= max_iter) break;

            dF(z, J);
            //printf("Jacobian J = [%f, %f; %f, %f]\n", J[0][0], J[0][1], J[1][0], J[1][1]);

            //  DGETRF computes an LU factorization of a general M-by-N matrix J using partial pivoting with row interchanges.
            info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, DIM, DIM, *J, DIM, ipiv);
            if (info != 0) 
            {
                fprintf(stderr, "Error: Matrix factorization failed.\n");
                exit(EXIT_FAILURE);
            }
            // DGETRI computes the inverse of a matrix using the LU factorization computed by DGETRF.
            info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, DIM, *J, DIM, ipiv);
            if (info != 0) 
            {
                fprintf(stderr, "Error: Matrix inversion failed.\n");
                exit(EXIT_FAILURE);
            }

            // Calculate delta_z = J_inv * F(z) using BLAS
            cblas_dgemv(CblasRowMajor, CblasNoTrans, DIM, DIM, 1.0, *J, DIM, Fz.data, 1, 0.0, delta_z, 1);

            // Update z: z = z - delta_z
            for (k = 0; k < DIM; k++) 
            {
                z[k] -= delta_z[k];
            }
            count++;
        }

        if (count >= max_iter) 
        {
            fprintf(stderr, "Warning: Maximum iterations reached without convergence.\n");
        }
        printf("Set %d: Iterations = %d, Solution: [", i + 1, count);
        for (l = 0; l < DIM; l++) 
        {
            printf("%.6lf", z[l]);  // Print the element without a comma
            if (l < DIM - 1)        // Check if it's not the last element
            {
                printf(", ");       // Print the comma and space if it's not the last element
            }
        }
        printf("]\n");
    }
    freeVector(&Fz);
}