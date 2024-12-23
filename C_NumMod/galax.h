/******************************************************/
//  GALAX - GSL Advanced Analytics eXtension
/******************************************************/
/*  
    $ gcc -c galax.c -o galax.o
    $ ar rcs libgalax.a galax.o
*/

#ifndef GALAX_H
#define GALAX_H

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_permutation.h>

// Type Definitions
typedef struct {
    gsl_matrix *gsl_matrix_ptr;
    unsigned int rows;
    unsigned int cols;
} Matrix;

typedef struct {
    gsl_vector *gsl_vector_ptr;
    unsigned int size;
} Vector;

// Matrix Functions
void    freeMatrix(Matrix *matrix);
Matrix  zerosMatrix(int rows, int cols);
Matrix  createMatrix(int rows, int cols, double *values);
void    printMatrix(const Matrix *mat);
Matrix  transposeMatrix(const Matrix *mat);
Matrix  inverseMatrix(const Matrix *mat);
Matrix  matmul(const Matrix *A, const Matrix *B);
Matrix  matadd(const Matrix *A, const Matrix *B);
Matrix  matscal(const Matrix *A, double scalar);
double  det(const Matrix *mat);
Matrix  matelem(const Matrix *A, const Matrix *B);
Matrix  matrand(int rows, int cols);
Matrix  compan(const Vector *coefficients);
void    eig(const Matrix *matrix);
Matrix  eigB(const Matrix *matrix, bool save_result);
double  normMatrix(const Matrix *mat, double p);
Matrix  linsolve_overdet(const Matrix *A, const Matrix *F);
Matrix roots(const Vector *coeff);

// Vector Functions
void    freeVector(Vector *vector);
Vector  zerosVector(int size);
Vector  createVector(int size, double *values);
void    printVector(const Vector *vector);
Vector  matvec(const Matrix *A, const Vector *x);
Vector  linsolve(const Matrix *A, const Vector *b);
Vector  vecadd(const Vector *A, const Vector *B);
Vector  vecscal(const Vector *A, double scalar);
Vector  vecrand(int dim);
Vector  vecelem(const Vector *A, const Vector *B);
double  norm(const Vector *vec, double p);
Vector  createArray(double start, double end, double step);
Vector  polyfitweighted(const Vector *x, const Vector *y, const Vector *w, int n);
Vector  linspace(double start, double end, int num);
//Vector  roots(const Vector *coeff);
Vector  polycoefs(const Vector *roots);
Vector  conv(const Vector *a, const Vector *b);

#endif // GALAX_H
