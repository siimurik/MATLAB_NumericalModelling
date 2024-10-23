// Compile and execute with:
// gcc -o yl3 yl3.c splash.c -llapacke -lblas
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <lapacke.h>
#include <cblas.h>
#include "splash.h"

int main()
{
	// Finding the roots of a polynomial x^3 − 8^x + 2 = 0
	float r[4] = {1.0, 0.0, 8.0, 2.0};
	int dim = sizeof(r) / sizeof(r[0]); 
	Matrix matCompan = compan(r, dim);
	printf("Companion matrix:\n");
	printMatrix(&matCompan);
	
	// The eigenvalues are the roots to this polynomial
	eig(&matCompan);

	freeMatrix(&matCompan);
	return 0;
}
