#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void printMatrixSize(const gsl_matrix *m) {
    size_t rows = m->size1;
    size_t cols = m->size2;
    printf("Matrix dimensions: %zu x %zu\n", rows, cols);
}

void printVectorSize(const gsl_vector *v) {
    size_t size = v->size;
    printf("Vector size: %zu\n", size);
}

int main() {
    gsl_matrix *matrix = gsl_matrix_alloc(3, 3);
    gsl_vector *vector = gsl_vector_alloc(3);

    printMatrixSize(matrix);
    printVectorSize(vector);

    gsl_matrix_free(matrix);
    gsl_vector_free(vector);

    return 0;
}
