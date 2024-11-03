/* Compile and execute with
    $ gcc yl5.c -o yl5 -lm -lblas -llapacke
    $ ./yl5
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <cblas.h> 
#include <lapacke.h>

// DIM defines the number of unknown variables
#define DIM 2  // For 2D systems, "2" marks the variables x and y

// Function prototypes
double normOld(double *vec, int dim, double p);
double norm(double *vec, int dim, double p);
void F(double *z, double *result);
void dF(double *z, double result[DIM][DIM]);
void G1(double *z, double *result);
void G2(double *z, double *result);
void G3(double *z, double *result);
void G4(double *z, double *result);
void SIM_Vec(void (*F)(double*, double*), void (**G_funcs)(double*, double*), double **xy_vec, 
            int vec_count, double epsilon, int max_iter);
void newtonVec(void (*F)(double*, double*), void (*dF)(double*, double[DIM][DIM]), 
            double **initial_values, int num_values, double epsilon, int max_iter);
void modNewtonVec(void (*F)(double*, double*), void (*dF)(double*, double[DIM][DIM]), 
            double **initial_values, int num_values, double epsilon, int max_iter);

int main() {
    double tol = 1.0E-6;
    int max_iter = 1000;

    // Initial guesses
    double xy1[2] = {-3.8, -2.1};
    double xy2[2] = {-1.6, 2.7};
    double xy3[2] = {2.2, 1.4};
    double xy4[2] = {2.8, -0.1};
    double *xy_vec[4] = {xy1, xy2, xy3, xy4};

    // G functions
    void (*G_funcs[4])(double*, double*) = {G1, G2, G3, G4};

    SIM_Vec(F, G_funcs, xy_vec, 3, tol, max_iter); // Use 3 to exclude G4 and xy4

    printf("Newton's Method:\n");
    newtonVec(F, dF, xy_vec, 4, tol, max_iter);

    printf("Modified Newton's Method:\n");
    modNewtonVec(F, dF, xy_vec, 4, tol, max_iter);

    return 0;
}

// Calculate any norm of a vector
double normOld(double *vec, int dim, double p) {
    double norm_value = 0.0;

    if (p == 1) {
        // 1-norm (Manhattan norm)
        for (int i = 0; i < dim; i++) {
            norm_value += fabs(vec[i]);
        }
    } else if (p == 2) {
        // 2-norm (Euclidean norm)
        for (int i = 0; i < dim; i++) {
            norm_value += vec[i] * vec[i];
        }
        norm_value = sqrt(norm_value);
    } else if (p == INFINITY) {
        // Infinity-norm (Maximum norm)
        for (int i = 0; i < dim; i++) {
            double abs_val = fabs(vec[i]);
            if (abs_val > norm_value) {
                norm_value = abs_val;
            }
        }
    } else {
        // p-norm
        for (int i = 0; i < dim; i++) {
            norm_value += pow(fabs(vec[i]), p);
        }
        norm_value = pow(norm_value, 1.0 / p);
    }

    return norm_value;
}

// Calculate any norm of a vector using BLAS
double norm(double *vec, int dim, double p) {
    double norm_value = 0.0;

    if (p == 1) {
        // 1-norm using BLAS
        norm_value = cblas_dasum(dim, vec, 1);
    } else if (p == 2) {
        // 2-norm using BLAS
        norm_value = cblas_dnrm2(dim, vec, 1);
    } else if (p == INFINITY) {
        // Infinity-norm (Maximum norm)
        for (int i = 0; i < dim; i++) {
            double abs_val = fabs(vec[i]);
            if (abs_val > norm_value) {
                norm_value = abs_val;
            }
        }
    } else {
        // p-norm
        for (int i = 0; i < dim; i++) {
            norm_value += pow(fabs(vec[i]), p);
        }
        norm_value = pow(norm_value, 1.0 / p);
    }

    return norm_value;
}

// Function F
void F(double *z, double *result) 
{
    result[0] = z[0] * z[0] + 2.0*z[1] - 8.0;
    result[1] = z[0] + z[1]*z[1] - z[1] - 3.0;
}

// First G function
void G1(double *z, double *result) {
    result[0] = -sqrt(8.0 - 2.0*z[1]);
    result[1] = -sqrt(z[1] - z[0] + 3.0);
}

// Second G function
void G2(double *z, double *result) 
{
    result[0] = -sqrt(8.0 - 2.0*z[1]);
    result[1] = sqrt(z[1] - z[0] + 3.0);
}

// Third G function
void G3(double *z, double *result) 
{
    result[0] = sqrt(8.0 - 2.0*z[1]);
    result[1] = sqrt(z[1] - z[0] + 3.0);
}

// Fourth G function
void G4(double *z, double *result) 
{
    result[0] =  sqrt(8.0 - 2.0*z[1]);
    result[1] = -sqrt(z[1] - z[0] + 3.0);
}

// The main SIM_Vec function
void SIM_Vec(void (*F)(double*, double*), void (**G_funcs)(double*, double*), double **xy_vec, 
                                            int vec_count, double epsilon, int max_iter) 
{
    for (int i = 0; i < vec_count; i++) 
    {
        double z[DIM];       // reserve memory for z0 and z1 values
        z[0] = xy_vec[i][0]; // This is why "xy_vec" is declared as a double pointer.
        z[1] = xy_vec[i][1]; // It points to the "xy" pair values inside the "xy_vec" 

        void (*G)(double*, double*) = G_funcs[i]; // from the array of G arrays, pick out one G function
        int count = 0;

        double Fz[DIM]; // reserve memory for F(z0) and F(z1)
        F(z, Fz);       // input z is saved into output Fz

        while (norm(Fz, DIM, 1) >= epsilon && count < max_iter) 
        {
            double new_z[DIM];
            G(z, new_z); // Update z
            z[0] = new_z[0];
            z[1] = new_z[1];

            F(z, Fz);
            count++;
        }

        if (count == max_iter) 
        {
            fprintf(stderr,"Warning: Maximum number of iterations reached without convergence.\n");
            exit(EXIT_FAILURE);
        }
        printf("Iterations = %d, Solution z = [%f, %f]\n", count, z[0], z[1]);
    }
}

// Define the Jacobian dF(z)
void dF(double *z, double result[DIM][DIM]) {
    result[0][0] = 2.0 * z[0];
    result[0][1] = 2.0;
    result[1][0] = 1.0;
    result[1][1] = 2.0 * z[1] - 1.0;
}

void newtonVec(void (*F)(double*, double*), void (*dF)(double*, double[DIM][DIM]), 
                   double **initial_values, int num_values, double epsilon, int max_iter) {
    for (int i = 0; i < num_values; i++) 
    {
        double z[DIM] = {initial_values[i][0], initial_values[i][1]};   // for x and y
        int count = 0;
        double Fz[DIM], delta_z[DIM];
        
        while (1) { // While TRUE, aka go on forever
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
            if (info != 0) {
                fprintf(stderr, "Error: Matrix factorization failed.\n");
                exit(EXIT_FAILURE);
            }
            // DGETRI computes the inverse of a matrix using the LU factorization computed by DGETRF.
            info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, DIM, *J, DIM, ipiv);
            if (info != 0) {
                fprintf(stderr, "Error: Matrix inversion failed.\n");
                exit(EXIT_FAILURE);
            }

            // Calculate delta_z = J_inv * F(z) using BLAS
            cblas_dgemv(CblasRowMajor, CblasNoTrans, DIM, DIM, 1.0, *J, DIM, Fz, 1, 0.0, delta_z, 1);

            // Update z: z = z - delta_z
            for (int j = 0; j < DIM; j++) {
                z[j] -= delta_z[j];
            }
            count++;
        }

        if (count >= max_iter) {
            fprintf(stderr, "Warning: Maximum iterations reached without convergence.\n");
        }
        printf("Set %d: Iterations = %d, Solution = [%.6f, %.6f]\n", i+1, count, z[0], z[1]);
    }
}

void modNewtonVec(void (*F)(double*, double*), void (*dF)(double*, double[DIM][DIM]), 
                     double **initial_values, int num_values, double epsilon, int max_iter) {
    for (int i = 0; i < num_values; i++) {
        double z[DIM] = {initial_values[i][0], initial_values[i][1]};
        double z0[DIM] = {z[0], z[1]};  // Initial z0 for Jacobian calculation
        int count = 0;
        double Fz[DIM], delta_z[DIM];

        while (1) {
            F(z, Fz);
            if (cblas_dasum(DIM, Fz, 1) < epsilon || count >= max_iter) break;

            double J[DIM][DIM];
            dF(z0, J);  // Only calculate Jacobian at initial z0

            // Invert J using LAPACK
            int ipiv[DIM];
            int info;
            info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, DIM, DIM, *J, DIM, ipiv);
            if (info != 0) {
                fprintf(stderr, "Error: Matrix factorization failed.\n");
                exit(EXIT_FAILURE);
            }
            info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, DIM, *J, DIM, ipiv);
            if (info != 0) {
                fprintf(stderr, "Error: Matrix inversion failed.\n");
                exit(EXIT_FAILURE);
            }

            // Calculate delta_z = J_inv * F(z) using BLAS
            cblas_dgemv(CblasRowMajor, CblasNoTrans, DIM, DIM, 1.0, *J, DIM, Fz, 1, 0.0, delta_z, 1);

            // Update z: z = z - delta_z
            for (int j = 0; j < DIM; j++) {
                z[j] -= delta_z[j];
            }
            count++;
        }

        if (count >= max_iter) {
            fprintf(stderr, "Warning: Maximum iterations reached without convergence.\n");
        }
        printf("Set %d: Iterations = %d, Solution = [%.6f, %.6f]\n", i+1, count, z[0], z[1]);
    }
}

