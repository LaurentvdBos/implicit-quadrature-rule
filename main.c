#include "matrix.h"

#include <stdio.h>

int main()
{
	struct matrix *mat = matrix_malloc(4, 3);
	
	double a[] = {12., 6., -4., -51., 167., 24., 4., -68., -41., 6., -12., 15.};
	matrix_set(mat, a);

	printf("a:\n");
	matrix_fprintf(stdout, mat, "%lf");

	struct matrix *q1 = matrix_malloc(4, 4);
	matrix_householder(mat, q1, 0);

	printf("\nq1:\n");
	matrix_fprintf(stdout, q1, "%lf");

	struct matrix *r1 = matrix_malloc(4, 3);
	matrix_mul(r1, q1, mat);

	printf("\nr1:\n");
	matrix_fprintf(stdout, r1, "%lf");

	struct matrix *q2 = matrix_malloc(4, 4);
	matrix_householder(r1, q2, 1);

	printf("\nq2:\n");
	matrix_fprintf(stdout, q2, "%lf");
	
	struct matrix *r2 = matrix_malloc(4, 3);
	matrix_mul(r2, q2, r1);
	
	printf("\nr2:\n");
	matrix_fprintf(stdout, r2, "%lf");

	struct matrix *q3 = matrix_malloc(4, 4);
	matrix_householder(r2, q3, 2);

	printf("\nq3:\n");
	matrix_fprintf(stdout, q3, "%lf");

	struct matrix *r3 = matrix_malloc(4, 3);
	matrix_mul(r3, q3, r2);

	printf("\nr3 (r):\n");
	matrix_fprintf(stdout, r3, "%lf");

	struct matrix *h = matrix_malloc(4, 4);

	//(Q1^T*Q2^T*Q3^T) = (Q2*Q1)^T*Q3^T = H^T*Q3^T
	matrix_mul(h, q2, q1);
	matrix_multt(q1, h, q3);

	printf("\nq:\n");
	matrix_fprintf(stdout, q1, "%lf");

	struct matrix *qr = matrix_malloc(4, 3);
	matrix_mul(qr, q1, r3);
	printf("\nq*r:\n");
	matrix_fprintf(stdout, qr, "%lf");

	return 0;
}
