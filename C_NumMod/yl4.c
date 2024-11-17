/* Compile and execute with
    $ gcc yl4.c -o yl4 -lm
    $ ./yl4
*/
#include <stdio.h>
#include <math.h>

// Simple iteration method
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
// Newton's iteration method
void newtonIter(double (*f)(double), double (*df)(double), double *x0, int *count, double tol, int max_iter)
{
    double x1 = *x0;
    *count = 0;
    while (fabs(f(x1)) >= tol && *count < max_iter)
    {
        *x0 = x1;
        x1 = *x0  - f(*x0)/df(*x0);
        (*count)++;
    }
    *x0 = x1;

    if (*count >= max_iter)
    {
        printf("Warning: Maximum iteration count reached.\n");
    }
}
// Modified Newton's iteration method
void modNewtonIter(double (*f)(double), double (*df)(double), double *x0, int *count, double tol, int max_iter)
{
    double x1    = *x0;
    *count       = 0;
    double df_x0 = df(*x0); 
    while (fabs(f(x1)) >= tol && *count < max_iter)
    {
        *x0 = x1;
        x1 = *x0  - f(*x0)/df_x0;
        (*count)++;
    }
    *x0 = x1;

    if (*count >= max_iter)
    {
        printf("Warning: Maximum iteration count reached.\n");
    }
}
// Secant iteration method
void secantIter(double (*f)(double), double *x, double x1, int *count, double tol, int max_iter)
{
    double x0 = *x; // if not done this way, pointer arithmetic will lead to failure
    double x2;
    *count = 0;
    while (fabs(f(x1)) >= tol && *count < max_iter)
    {
        x2 = x1 - f(x1) * (x1 - x0)/(f(x1) - f(x0));
        x0 = x1;
        x1 = x2;
        (*count)++;
    }
    *x = x1;

    if (*count >= max_iter)
    {
        printf("Warning: Maximum iteration count reached.\n");
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

double f (double x)
{
    return x*x*x*x + 2.0*x*x - x - 3.0;
}

double df (double x)
{
    return 4.0*x*x*x + 4.0*x - 1.0;
}

double g (double x)
{
    return (5.0 - 2.0*x*x*x + 11.7*x*x)/17.7;   
}

int main()
{
    double x0    = -1.0;
    double tol   = 1E-4;
    int max_iter = 1000;
    int count = 0;

    newtonIter(f, df, &x0, &count, tol, max_iter);
    printf("Newton's method: x = %f, Iterations: %d\n", x0, count);

    // x0 got overwritten in the previous code so reassign new value
    x0 = 1.0;
    modNewtonIter(f, df, &x0, &count, tol, max_iter);
    printf("Modified Newton's method: x = %f, Iterations: %d\n", x0, count);

    x0 = 1.0; double x1 = 1.5; count = 0;
    secantIter(f, &x0, x1, &count, tol, max_iter);
    printf("Secant method: x = %f, Iterations: %d\n", x0, count);

    double x_init = 0.1; count = 0;
    tol = 1.E-6;   max_iter = 1000;
    sim(g, &x_init, &count, tol, max_iter);
    printf("SIM: x = %f, Iterations: %d\n", x_init, count);

    return 0;
}
