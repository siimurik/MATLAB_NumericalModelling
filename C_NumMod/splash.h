#ifndef SPLASH_H
#define SPLASH_H

/****************************************************************************/
/*  Officaly known as                                                       */
/*    SPLASH - Single Precision for Linear Algebra and Scientific Handling  */
/*--------------------------------------------------------------------------*/
/*  Or unoffically known as                                                 */
/*    SPLASH - Siim’s Package for Linear Algebra and Scientific Handling    */
/****************************************************************************/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>
#include <cblas.h>
#include <omp.h> 

typedef struct {
    float *data;
    unsigned int rows;
    unsigned int cols;
} Matrix;

typedef struct {
    float *data;
    unsigned int size;
} Vector;

// Matrix operations
void    freeMatrix(Matrix *matrix);
Matrix  createMatrix(int rows, int cols, float *values);
void    printMatrix(const Matrix *mat);
Matrix  transposeMatrix(const Matrix *mat);
Matrix  transposeMatrixPara(const Matrix *mat);
Matrix  inverseMatrix(const Matrix *mat);
Matrix  matmul(const Matrix *A, const Matrix *B);
Matrix  matscal(const Matrix *A, float scalar);
float   det(const Matrix *mat);
Matrix  matelem(const Matrix *A, const Matrix *B);
Matrix  matelemPara(const Matrix *A, const Matrix *B);
Matrix  linsolve_overdet(const Matrix *A, const Matrix *F);
Matrix  matrand(int rows, int cols);
Matrix  matrandPara(int rows, int cols);
Matrix  compan(const float *coefficients, int size);
void    eig(const Matrix *matrix);
float   normMatrix(const Matrix *mat, float p);

// Vector operations
void    freeVector(Vector *vector);
Vector  createVector(int size, float *values);
void    printVector(const Vector *vector);
Vector  matvec(const Matrix *A, const Vector *x);
Vector  linsolve(const Matrix *A, const Vector *b);
Vector  vecadd(const Vector *A, const Vector *B);
Vector  vecscal(const Vector *A, float scalar);
Vector  vecrand(int dim);
Vector  vecrandPara(int dim);
Vector  vecelem(const Vector *A, const Vector *B);
Vector  vecelemPara(const Vector *A, const Vector *B);
float   norm(const Vector *vec, float p);
Vector  createArray(float start, float end, float step);
Vector  polyfitweighted(const Vector *x, const Vector *y, const Vector *w, int n);
Vector  linspace(float start, float end, int num);
Vector  roots(const Vector *coeff);
Vector  polycoefs(const Vector *roots);
Vector  conv(const Vector *a, const Vector *b);

#endif