#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main()
{
	srand(time(0));

	struct matrix *mat = matrix_malloc(10, 11);
	
	double a[mat->n*mat->m];
	for (int i = 0; i < mat->n*mat->m; i++) {
		a[i] = (double)rand()/(double)RAND_MAX*2.-1.;
	}
	matrix_set(mat, a);

	printf("a:\n");
	matrix_fprintf(stdout, mat, "%lf");

	struct matrix *q = matrix_malloc(mat->n, mat->n);
	struct matrix *r = matrix_malloc(mat->n, mat->m);
	struct matrix *qr = matrix_malloc(mat->n, mat->m);

	matrix_qr(mat, q, r);
	printf("\nq:\n");
	matrix_fprintf(stdout, q, "%lf");
	printf("\nr:\n");
	matrix_fprintf(stdout, r, "%lf");

	matrix_mul(qr, q, r);
	printf("\nq*r:\n");
	matrix_fprintf(stdout, qr, "%lf");

	matrix_free(mat);
	matrix_free(q);
	matrix_free(r);
	matrix_free(qr);

	return 0;
}
