/* Compile and execute with
    $ gcc yl6.c -o yl6 -lm -lblas -llapacke
    $ ./yl6
*/

#include <stdio.h>  // for printf
#include <stdlib.h> // for stderr
#include <string.h> // for memcpy
#include <math.h>   // for pow()
#include <cblas.h>  // for cblas_*
#include <lapacke.h>// for LAPACKE_*

#define DIM 3   // x, y, z

typedef struct 
{
    double *data;
    int rows;
    int cols;
} Matrix;

double  norm(double *vec, int dim, double p);
double  normMatrix(const Matrix *mat, double p);
Matrix  createMatrix(int rows, int cols, double *values);
void    freeMatrix(Matrix *matrix);
void    F(double *z, double *result);
void    dF(double *z, double J[DIM][DIM]);
void    newtonVec(void (*F)(double*, double*), void (*dF)(double*, double[DIM][DIM]), 
            double **initial_values, int num_values, double epsilon, int max_iter);

int main()
{
    printf("Exercise 1:\n");
    double a[4] = {6.0, -1.0, -3.0, 8};
    double a1 = norm(a, 4, 1);
    double a2 = norm(a, 4, 2);
    double a5 = norm(a, 4, 5);
    double ainf = norm(a, 4, INFINITY);
    printf("||a||_1   = %.4f\n", a1);
    printf("||a||_2   = %.4f\n", a2);
    printf("||a||_5   = %.4f\n", a5);
    printf("||a||_inf = %.4f\n", ainf);

    printf("\nExercise 2:\n");
    double values[9] = { 8.0,  7.0, -2.0,
                        -6.0,  4.0,  3.0,
                         4.0, -1.0,  5.0};
    Matrix A = createMatrix(3, 3, values);

    printf("1-norm: %.4f\n", normMatrix(&A, 1));
    printf("Frobenius norm: %.4f\n", normMatrix(&A, 2));
    printf("Infinity-norm: %.4f\n", normMatrix(&A, INFINITY));

    printf("\nExercise 3:\n");
    double tol = 1.0E-6;
    int max_iter = 1000;
    double xyz1[3] = {0.5,   0.5,  0.5};
    double xyz2[3] = {-1.0, -1.0, -1.0};
    double xyz3[3] = {-0.7,  0.5, -0.2};
    double *xyz_vec[3] = {xyz1, xyz2, xyz3};
    printf("Newton's Method:\n");
    newtonVec(F, dF, xyz_vec, 3, tol, max_iter);

    printf("\nExercise 4:\n");
    double Fout[DIM];
    F(xyz3, Fout);
    for (int i = 0; i < 3; i++)
    {
        printf("F%d(%3.1f, %3.1f, %3.1f) = %.4f\n", i+1, xyz3[0], xyz3[1], xyz3[2], Fout[i]);
    }

    freeMatrix(&A);

    return 0;
}

///////////////////////////////////////////////////////////////////////

// Calculate any norm of a vector using BLAS
double norm(double *vec, int dim, double p) 
{
    double norm_value = 0.0;

    if (p == 1) 
    {
        // 1-norm using BLAS
        norm_value = cblas_dasum(dim, vec, 1);
    } else if (p == 2) 
    {
        // 2-norm using BLAS
        norm_value = cblas_dnrm2(dim, vec, 1);
    } else if (p == INFINITY) 
    {
        // Infinity-norm (Maximum norm)
        for (int i = 0; i < dim; i++) 
        {
            double abs_val = fabs(vec[i]);
            if (abs_val > norm_value) 
            {
                norm_value = abs_val;
            }
        }
    } else 
    {
        // p-norm
        for (int i = 0; i < dim; i++) 
        {
            norm_value += pow(fabs(vec[i]), p);
        }
        norm_value = pow(norm_value, 1.0 / p);
    }

    return norm_value;
}

// Calculate any norm of a matrix using BLAS
double normMatrix(const Matrix *mat, double p)
{
    double norm_value = 0.0;

    if (p == 1) {
        // 1-norm (maximum column sum)
        for (int j = 0; j < mat->cols; j++) 
        {
            double col_sum = 0.0;
            for (int i = 0; i < mat->rows; i++) 
            {
                col_sum += fabs(mat->data[i * mat->cols + j]);
            }
            if (col_sum > norm_value) 
            {
                norm_value = col_sum;
            }
        }
    } else if (p == 2) 
    {
        // Frobenius norm (similar to Euclidean norm for vectors)
        norm_value = cblas_dnrm2(mat->rows * mat->cols, mat->data, 1);
    } else if (p == INFINITY) 
    {
        // Infinity-norm (maximum row sum)
        for (int i = 0; i < mat->rows; i++) 
        {
            double row_sum = 0.0;
            for (int j = 0; j < mat->cols; j++) 
            {
                row_sum += fabs(mat->data[i * mat->cols + j]);
            }
            if (row_sum > norm_value) 
            {
                norm_value = row_sum;
            }
        }
    } else 
    {
        fprintf(stderr, "Unsupported norm type for matrices.\n");
        exit(EXIT_FAILURE);
    }

    return norm_value;
}


// Function to create and initialize a matrix with specific values
Matrix createMatrix(int rows, int cols, double *values)
{
    Matrix matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.data = (double *)malloc(rows * cols * sizeof(double)); // Allocate memory

    if (matrix.data == NULL) 
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Populate the matrix with the provided values
    memcpy(matrix.data, values, matrix.rows * matrix.cols * sizeof(double));

    return matrix;
}

// Free allocated Matrix
void freeMatrix(Matrix *matrix) 
{
    free(matrix->data);
    matrix->data = NULL;
} 

// System of nonlinear equations for x, y and z (z0, z1, z2)
void F(double *z, double *result)
{
    result[0] = z[0]*z[0] + z[1]*z[1] + z[2]*z[2] - 1.0;
    result[1] = exp(z[0]) + z[1]*cos(z[2]) - 1.0;
    result[2] = sin(z[0]*z[1]) + z[1]*z[1]*z[1] - z[2];
}

// Define the Jacobian dF(z)
void dF(double *z, double J[DIM][DIM])
{
    J[0][0] = 2.0*z[0];            J[0][1] = 2.0*z[1];                            J[0][2] = 2.0*z[2];
    J[1][0] = exp(z[0]);           J[1][1] = cos(z[2]);                           J[1][2] = -z[1]*sin(z[2]);
    J[2][0] = z[1]*cos(z[0]*z[1]); J[2][1] = z[0]*cos(z[0]*z[1]) + 3.0*z[1]*z[1]; J[2][2] = -1.0;
}

// NOTE: This is a bad implementation bc values get redefied in the various loops.
// Also, this may cause the calculations to take more iterations as well.
// Newton's iteration method for vector inputs
void newtonVec(void (*F)(double*, double*), void (*dF)(double*, double[DIM][DIM]), 
                   double **initial_values, int num_values, double epsilon, int max_iter) 
{
    for (int i = 0; i < num_values; i++) 
    {
        double z[DIM];
        for (int j = 0; j < DIM; j++)
        {
            z[j] = initial_values[i][j];
        }
        int count = 0;
        double Fz[DIM], delta_z[DIM];
        
        while (1) 
        { // While TRUE, aka go on forever
            F(z, Fz);
            // The main breaking conditions that stop the code from running forever
            if (norm(Fz, DIM, 1) < epsilon || count >= max_iter) break;

            // Set up the Jacobi matrix J
            double J[DIM][DIM]; // in LAPACK/BLAS needs a pointer
            dF(z, J);

            // Invert J using LAPACK
            int ipiv[DIM];
            int info;
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
            cblas_dgemv(CblasRowMajor, CblasNoTrans, DIM, DIM, 1.0, *J, DIM, Fz, 1, 0.0, delta_z, 1);

            // Update z: z = z - delta_z
            for (int k = 0; k < DIM; k++) 
            {
                z[k] -= delta_z[k];
            }
            count++;
        }

        if (count >= max_iter) 
        {
            fprintf(stderr, "Warning: Maximum iterations reached without convergence.\n");
        }
        printf("Set %d: Iterations = %d, Solution = [%.6f, %.6f, %.6f]\n", i+1, count, z[0], z[1], z[2]);
    }
}