#include "matrix.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <lapacke.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

struct matrix *matrix_malloc(int n, int m)
{
	struct matrix *mat = malloc(sizeof(struct matrix));
	mat->n = n;
	mat->m = m;
	mat->lda = m;
	if (n > 0 && m > 0) {
		mat->a = malloc(sizeof(double)*n*m);
	} else {
		mat->a = NULL;
	}
	return mat;
}

void matrix_free(struct matrix *mat)
{
	if (mat) {
		free(mat->a);
		free(mat);
	}
}

void matrix_set(struct matrix *mat, const double *a)
{
	for (int i = 0; i < mat->n*mat->m; i++) {
		mat->a[i] = a[i];
	}
}

void matrix_swap(struct matrix *a, struct matrix *b)
{
	struct matrix tmp = *a;
	*a = *b;
	*b = tmp;
}

void matrix_fprintf(FILE *f, const struct matrix *mat, const char *fmt)
{
	for (int i = 0; i < mat->n; i++) {
		for (int j = 0; j < mat->m; j++) {
			fprintf(f, fmt, mat->a[i*mat->lda + j]);
			fputc(' ', f);
		}
		fputc('\n', f);
	}
}


// mat = b*c
void matrix_mul(struct matrix *mat, const struct matrix *b, const struct matrix *c)
{
	assert(mat->n == b->n);
	assert(mat->m == c->m);
	assert(b->m == c->n);

	for (int i = 0; i < mat->n; i++) {
		for (int j = 0; j < mat->m; j++) {
			mat->a[i*mat->lda + j] = 0;

			for (int k = 0; k < c->n; k++) {
				mat->a[i*mat->lda + j] += b->a[i*b->lda + k]*c->a[k*c->lda + j];
			}
		}
	}
}

// Determine n null vectors of matrix mat and store in q. The current
// implementation trashes mat. Null vectors are the rows of q. If the number of
// rows of q is larger than the number of null vectors, no error is reported.
void matrix_null(struct matrix *mat, struct matrix *q)
{
	assert(q->n == mat->m);
	assert(q->m <= mat->m);

	double *tau = calloc(min(mat->n, mat->m), sizeof(double));
	lapack_int *pvt = calloc(mat->n, sizeof(lapack_int));

	double lwork, *work;
	lapack_int info;

	// Query for ideal work size
	info = LAPACKE_dgeqp3_work(LAPACK_COL_MAJOR, mat->m, mat->n, mat->a, mat->lda, pvt, tau, &lwork, -1);
	assert(info == 0);

	// Determine pivoted QR decomposition
	work = malloc(sizeof(double)*(lapack_int)lwork);
	info = LAPACKE_dgeqp3_work(LAPACK_COL_MAJOR, mat->m, mat->n, mat->a, mat->lda, pvt, tau, work, (lapack_int)lwork);
	assert(info == 0);

	// Query for ideal work size
	info = LAPACKE_dormqr_work(LAPACK_COL_MAJOR, 'R', 'T', q->m, q->n, min(mat->n, mat->m), mat->a, mat->lda, tau, q->a, q->lda, &lwork, -1);
	assert(info == 0);
	
	// Construct identity matrix in q, selecting the last columns from the decomposition
	for (int i = 0; i < q->n; i++) {
		for (int j = 0; j < q->m; j++) {
			q->a[i*q->lda + j] = (q->n-i == q->m-j ? 1. : 0.);
		}
	}

	// Obtain the last rows
	work = realloc(work, sizeof(double)*(lapack_int)lwork);
	info = LAPACKE_dormqr_work(LAPACK_COL_MAJOR, 'R', 'T', q->m, q->n, min(mat->n, mat->m), mat->a, mat->lda, tau, q->a, q->lda, work, (lapack_int)lwork);
	assert(info == 0);

	free(work);
	free(pvt);
	free(tau);
}
