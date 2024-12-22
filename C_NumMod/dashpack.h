#ifndef DASHPACK_H
#define DASHPACK_H

/************************************************************************************/
/*  Officaly known as                                                               */
/*    DASHPACK - Double-precision Algebraic Solutions for High-performance PACKage  */
/************************************************************************************/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>
#include <cblas.h>
#include <omp.h> 

typedef struct {
    double *data;
    int rows;
    int cols;
} Matrix;

typedef struct {
    double *data;
    int size;
} Vector;

// Matrix operations
void    freeMatrix(Matrix *matrix);
Matrix  createMatrix(int rows, int cols, double *values);
void    printMatrix(const Matrix *mat);
Matrix  transposeMatrix(const Matrix *mat);
Matrix  transposeMatrixPara(const Matrix *mat);
Matrix  inverseMatrix(const Matrix *mat);
Matrix  matmul(const Matrix *A, const Matrix *B);
Matrix  matscal(const Matrix *A, double scalar);
double  det(const Matrix *mat);
Matrix  matelem(const Matrix *A, const Matrix *B);
Matrix  matelemPara(const Matrix *A, const Matrix *B);
Matrix  linsolve_overdet(const Matrix *A, const Matrix *F);
Matrix  matrand(int rows, int cols);
Matrix  matrandPara(int rows, int cols);
Matrix  compan(const double *coefficients, int size);
void    eig(const Matrix *matrix);
double  normMatrix(const Matrix *mat, double p);

// Vector operations
void    freeVector(Vector *vector);
Vector  createVector(int size, double *values);
void    printVector(const Vector *vector);
Vector  matvec(const Matrix *A, const Vector *x);
Vector  linsolve(const Matrix *A, const Vector *b);
Vector  vecadd(const Vector *A, const Vector *B);
Vector  vecscal(const Vector *A, double scalar);
Vector  vecrand(int dim);
Vector  vecrandPara(int dim);
Vector  vecelem(const Vector *A, const Vector *B);
Vector  vecelemPara(const Vector *A, const Vector *B);
double  norm(const Vector *vec, double p);
Vector  createArray(double start, double end, double step);


#endif