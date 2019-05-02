#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main()
{
	srand(time(0));

	int n = 5;
	int m = 7;
	double *a = malloc(n*m*sizeof(double));
	for (int i = 0; i < n*m; i++)
		a[i] = (double) rand() / (double) RAND_MAX * 2. - 1.;
	
	struct matrix *mat = matrix_malloc(n, m);
	matrix_set(mat, a);

	matrix_fprintf(stdout, mat, "%.10f");

	struct matrix *q = matrix_malloc(m, m-n);
	matrix_null(mat, q);

	printf("\n");
	matrix_fprintf(stdout, mat, "%.10f");

	printf("\n");
	matrix_fprintf(stdout, q, "%.10f");

	matrix_free(q);
	matrix_free(mat);

	free(a);

	return 0;
}
