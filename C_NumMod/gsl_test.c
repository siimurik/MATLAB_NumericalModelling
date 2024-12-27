/* Compile and execute with:
    $ gcc -c galax.c -o galax.o
    $ ar rcs libgalax.a galax.o
    $ gcc gsl_test.c -o gtest -L. -lopenblas -lgalax -lgsl -lm
*/

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "galax.h"

#define true 1
#define false 0

double fun(double x, void *params){
    (void)params;
    return (x*x)/(1.0 + sin(x) + cos(x));
}

// Define the Van der Pol oscillator equations
int vdp1(double t, const double y[], double dydt[], void *params) {
    (void)(t);  // Avoid unused parameter warning
    double mu = *(double *)params;
    dydt[0] = y[1];
    dydt[1] = (1 - y[0] * y[0]) * y[1] - y[0];
    return GSL_SUCCESS;
}

// Define the new ODE system
// https://math.stackexchange.com/questions/490208/a-difficult-differential-equation-y2x4y-fracdydx-1-4xy2x2
int diff_eq(double x, const double y[], double dydx[], void *params) { 
    (void)(x); // Avoid unused parameter warning 
    dydx[0] = ((1.0 - 4.0*x*y[0]*y[0])*x*x) / (y[0] * (2*x*x*x*x + y[0]));  
    return GSL_SUCCESS;
} 

int main()
{
    printf("\nBasic Linear Algebra.\n");
    int rows, cols;
    rows = 3;
    cols = 3;

    // Declare and print values for the A matrix
    double mat_values[9] = {3.0, 12.0, 52.0, 
                            4.0, 6.0, -11.0, 
                            -2.0, 7.0, 2.0};
    Matrix A = createMatrix(rows, cols, mat_values);
    printf("Matrix A:\n");
    //printM(&A);
    printMatrix(&A);

    // Declare and print values for the b vector
    double vec_values[3] = {13.0, -2.0, 5.0};
    int size = sizeof(vec_values) / sizeof(vec_values[0]); // Calculate the size of the vector
    //printf("size = %d\n", size);
    Vector b = createVector(size, vec_values);
    printf("Vector b: ");
    printVector(&b);

    // Find the transposed and inverse matrix of matrix A
    Matrix A_T   = transposeMatrix(&A);
    printf("\nMatrix A^T:\n");
    printMatrix(&A_T);
    Matrix inv_A = inverseMatrix(&A);
    printf("\nMatrix A^(-1):\n");
    printMatrix(&inv_A);

    // B = 2*A^T + A^-1
    Matrix twoTimesAT = matscal(&A_T, 2.0);
    Matrix B = matadd( &twoTimesAT, &inv_A );
    printf("\nMatrix B = 2*A^T + A^-1: \n");
    printMatrix(&B);

    // C = A*B
    Matrix C = matmul(&A, &B);
    printf("\nMatrix C: \n");
    printMatrix(&C);
    
    // det(A)
    printf("\nDeterminant of matrix A: ");
    double det_A = det(&A);
    printf("%.4f.\n\n", det_A);
    
    Vector d = matvec(&A, &b);
    printf("Vector d: \n");
    printVector(&d);

    // D = A.*B
    Matrix D = matelem(&A, &B);
    printf("\nMatrix D:\n");
    printMatrix(&D);

    Matrix minus_B = matscal(&B, -1.0);
    printf("\nMatrix -B:\n");
    printMatrix(&minus_B);

    Matrix E = matadd(&A, &minus_B);
    printf("\nMatrix E:\n");
    printMatrix(&E);

    Vector x = linsolve(&A, &b);
    printf("\nx = "); printVector(&x);

    printf("\n");
    eig(&A);
    //Matrix eig_D = eigB(&D, false);
    //printf("\nEigenvalues of matrix D:\n");
    //printMatrix(&eig_D);

    printf("\nNorms of a matrix.");
    int n = 4;
    double norms[n];
    norms[0] = normMatrix(&A, INFINITY);
    norms[1] = normMatrix(&B, INFINITY);
    norms[2] = normMatrix(&C, INFINITY);
    norms[3] = normMatrix(&D, INFINITY);
    for (int i = 0; i < n; i++){
        printf("\nnorm[%d] = %.4f", i, norms[i]);
    }
    printf("\n");

    printf("\nVector Scalar Multiplication and Addition.");
    double scal = -1.0;
    Vector minus_d = vecscal(&d, scal);
    Vector x_minus_d = vecadd(&x, &minus_d);
    printf("\nminus_d = "); printVector(&minus_d);
    printf("x_minus_d = "); printVector(&x_minus_d);

    printf("\nVector norms.\n");
    Vector a1 = createArray(4.0, 8.0, 1.0);
    Vector b1 = createArray(88.0, 92.0, 1.0);
    printf("a1 = "); printVector(&a1);
    printf("b1 = "); printVector(&b1);

    Vector elem_prod = vecelem(&a1, &b1);
    printf("a1 * b1 = "); printVector(&elem_prod);

    Vector c1 = createArray(101.0, 105.0, 1.0);
    
    double norm_a1 = norm(&a1, 1);
    double norm_b1 = norm(&a1, 2);
    double norm_c1 = norm(&a1, INFINITY);
    printf("\n|a1|_1 = %.4f",     norm_a1);
    printf("\n|b1|_2 = %.4f",     norm_b1);
    printf("\n|c1|_inf = %.4f\n", norm_c1);

    printf("\n");
    printf("\nRoots of a polynomial.\n");
    Vector r = vecrand(3);
    printf("r = "); printVector(&r);

    Matrix comp = compan(&r);
    eig(&comp);
    
    Matrix coef = roots(&r);
    printf("coef = "); printMatrix(&coef);

    printf("\nPolynomial Coefficients based on roots.");
    double pData[3] = {-3.0, 0.0, 3.0};
    Vector p = createVector(3, pData);
    Vector p_coef = polycoefs(&p);
    printf("\np_coef = "); printVector(&p_coef);

    printf("\nConvolution.\n");
    Vector p2 = createArray(2.0, 5.0, 1.0);
    printf("p2 = "); printVector(&p2);

    Vector convolution = conv(&p_coef, &p2);
    printf("conv = "); printVector(&convolution);
    
    printf("\nIntegration.");
    double pi = 4.0*atan(1.0); 
    double xmin = 0.0;
    double xmax = pi/2.0;
    double result = integral(fun, xmin, xmax);
    printf("\nint_xmin^xmax (x^2)/(1.0 + sin(x) + cos(x)), xmin = 0.0, xmax = pi/2");
    printf("\nresult = %.10f\n", result);

    printf("\nThe van der Pol equation:\ny'' - μ(1-y^2)y' + y = 0\n");
    // ODE parameters
    double mu = 1.0;
    double t_start = 0.0;
    double y[2] = {2.0, 0.0};
    double t_end = 2.0;
    double dt = 0.1;  // Time step for printing results
    size_t dimension = 2;  // Dimension of the ODE system

    // Solve the ODE with default parameters
    ode45(vdp1, &mu, dimension, &t_start, t_end, y, dt);
    
    printf("\n");
    printf("\nDifferential equation:\ny(2x^4+y)dy/dx = (1-4xy^2)x^2\n");
    // Difficult ODE problem
    double t2 = 0.0;        // This value get modified; hence we have to point to it.
    double y2[1] = {1.0};   // Initial condition
    double tend2 = 1.0;     // End time
    double dt2 = 0.1;       // Time step for printing results
    size_t dim = 1;         // Dimension of the ODE system

    // Solve the ODE with default parameters
    ode45(diff_eq, NULL, dim, &t2, tend2, y2, dt2);
    
    printf("\nLeast squares method:\n");
    double s_data[] = {1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
    double y_data[] = {9.08, 10.43, 11.9, 13.48, 15.19, 17.03, 19.01, 21.13, 23.39};
    double w_data[] = {1.0, 1.0, 2.0, 5.0, 1.0, 4.0, 2.0, 2.0, 1.0};
    int m = sizeof(s_data) / sizeof(s_data[0]);
    Vector sM = createVector(m, s_data); 
    Vector yM = createVector(m, y_data); 
    Vector wM = createVector(m, w_data);
    
    printf("s = "); printVector(&sM);
    printf("y = "); printVector(&yM);
    printf("w = "); printVector(&wM);
    
    Vector p3 = polyfitweighted(&sM, &yM, &wM, 3);
    printf("\np3 = "); printVector(&p3);

    printf("\nInterpolation.\n");
    double distData[] = {0.0, 0.25,  0.5, 0.75, 1.0, 1.25};
    double t_data[] = {0.0, 25.0, 49.4, 73.0, 96.4, 119.4};
    int d_size = sizeof(distData) / sizeof(distData[0]);
    Vector di = createVector(d_size, distData);
    Vector t = createVector(d_size, t_data);
    printf("di = "); printVector(&di);
    printf("t = "); printVector(&t);


    Vector S1_0 = interp1(&di, &t, &di, "linear");
    printf("S1_0 = "); printVector(&S1_0);

    Vector dd = createArray(0.0, 1.25, 0.001);
    printf("dd = "); printVector(&dd);

    Vector S3_2 = interp1(&di, &t, &dd, "spline");
    printf("S3_2 = "); printVector(&S3_2);
    
    double s32valdata[] = {0.45};
    Vector s32val = createVector(1, s32valdata);
    Vector S3_2_val = interp1(&di, &t, &s32val, "linear");
    printf("S3_2(%.2f)= ", s32valdata[0]); printVector(&S3_2_val);

    
    freeMatrix(&A);
    freeMatrix(&B);
    freeMatrix(&minus_B);
    freeMatrix(&C);
    freeMatrix(&D);
    freeMatrix(&E);
    freeMatrix(&A_T);
    freeMatrix(&inv_A);
    freeMatrix(&twoTimesAT);
    freeVector(&b);
    freeVector(&d);
    freeVector(&x);
    //freeMatrix(&eig_D);
    freeVector(&minus_d);
    freeVector(&x_minus_d);
    freeVector(&r);
    freeMatrix(&comp);
    //freeMatrix(&eig_comp);
    freeMatrix(&coef);
    freeVector(&a1);
    freeVector(&b1);
    freeVector(&c1);
    freeVector(&elem_prod);
    freeVector(&p);
    freeVector(&p_coef);
    freeVector(&p2);
    freeVector(&convolution);
    freeVector(&sM);
    freeVector(&yM);
    freeVector(&wM);
    freeVector(&p3);
    freeVector(&di);
    freeVector(&t);
    freeVector(&S1_0);
    freeVector(&dd);
    freeVector(&S3_2);
    freeVector(&s32val);
    freeVector(&S3_2_val);


    return 0;
}
