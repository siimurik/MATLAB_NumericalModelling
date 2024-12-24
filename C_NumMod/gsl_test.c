/* Compile and execute with:
    $ gcc -c galax.c -o galax.o
    $ ar rcs libgalax.a galax.o
    $ gcc gsl_test.c -o gtest -L. -lgalax -lgsl -lm
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define true 1
#define false 0

//#include <math.h>
#include "galax.h"

int main()
{
    printf("\nExercise 6.\n");
    int rows, cols;
    rows = 3;
    cols = 3;

    // Declare and print values for the A matrix
    double mat_values[9] = {3, 12, 52, 4, 6, -11, -2, 7, 2};
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
    eigB(&A, true);
    //Matrix eig_D = eigB(&D, false);
    //printf("\nEigenvalues of matrix D:\n");
    //printMatrix(&eig_D);

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

    double scal = -1.0;
    Vector minus_d = vecscal(&d, scal);
    Vector x_minus_d = vecadd(&x, &minus_d);
    printf("\nminus_d = "); printVector(&minus_d);
    printf("x_minus_d = "); printVector(&x_minus_d);

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
    Vector r = vecrand(3);
    printf("r = "); printVector(&r);

    Matrix comp = compan(&r);
    eig(&comp);
    
    Matrix coef = roots(&r);
    printf("coef = "); printMatrix(&coef);


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
    return 0;
}