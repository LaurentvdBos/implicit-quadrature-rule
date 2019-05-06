#include "matrix.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lapacke.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

double *matrix_workspace = NULL;
lapack_int nwork = 0;

static inline void workspace_ensure(lapack_int n)
{
	if (n > nwork) {
		nwork = n;
		matrix_workspace = realloc(matrix_workspace, nwork*sizeof(double));
	}
}

struct matrix *matrix_malloc(int n, int m)
{
	struct matrix *mat = malloc(sizeof(struct matrix));
	mat->n = n;
	mat->m = m;
	mat->nrows = n;
	mat->ncols = m;
	mat->tau = NULL;
	mat->pvt = NULL;
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
		free(mat->tau);
		free(mat->pvt);
		free(mat->a);
		free(mat);
	}
}

void matrix_set(struct matrix *mat, const double *a)
{
	for (int i = 0; i < mat->n; i++) {
		for (int j = 0; j < mat->m; j++) {
			mat->a[i*mat->ncols + j] = a[i*mat->m + j];
		}
	}
}

// mat = b
void matrix_copy(struct matrix *mat, const struct matrix *b)
{
	assert(mat->n == b->n);
	assert(mat->m == b->m);

	for (int i = 0; i < mat->n; i++) {
		for (int j = 0; j < mat->m; j++) {
			mat->a[i*mat->ncols + j] = b->a[i*b->ncols + j];
		}
	}

	if (b->tau) {
		mat->tau = realloc(mat->tau, min(mat->n, mat->m)*sizeof(double));
		memcpy(mat->tau, b->tau, min(mat->n, mat->m)*sizeof(double));
	}
	if (b->pvt) {
		mat->pvt = realloc(mat->pvt, mat->n*sizeof(lapack_int));
		memcpy(mat->pvt, b->pvt, mat->n*sizeof(lapack_int));
	}
}

void matrix_resize(struct matrix *mat, const int n, const int m)
{
	if (mat->ncols >= m) {
		if (mat->nrows < n) {
			mat->a = realloc(mat->a, sizeof(double)*mat->ncols*n);
			mat->nrows = n;
		}
	} else {
		if (mat->nrows < n) {
			mat->nrows = n;
		}
		mat->a = realloc(mat->a, sizeof(double)*m*mat->nrows);
		for (int i = min(mat->n, n)-1; i >= 0; i--) {
			for (int j = mat->m-1; j >= 0; j--) {
				mat->a[i*m + j] = mat->a[i*mat->ncols + j];
			}
		}
		mat->ncols = m;
	}
	mat->n = n;
	mat->m = m;

	// A resize always invalidates a QR or LU decomposition
	free(mat->tau);
	free(mat->pvt);
	mat->tau = NULL;
	mat->pvt = NULL;
}

void matrix_shrink_to_fit(struct matrix *mat)
{
	if (mat->ncols > mat->m) {
		for (int i = 0; i < mat->n; i++) {
			for (int j = 0; j < mat->m; j++) {
				mat->a[i*mat->m + j] = mat->a[i*mat->ncols + j];
			}
		}
	}
	mat->ncols = mat->m;
	mat->nrows = mat->n;
	mat->a = realloc(mat->a, mat->n*mat->m*sizeof(double));
}

void matrix_fprintf(FILE *f, const struct matrix *mat, const char *fmt)
{
	for (int i = 0; i < mat->n; i++) {
		for (int j = 0; j < mat->m; j++) {
			fprintf(f, fmt, mat->a[i*mat->ncols + j]);
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
			mat->a[i*mat->ncols + j] = 0;

			for (int k = 0; k < c->n; k++) {
				mat->a[i*mat->ncols + j] += b->a[i*b->ncols + k]*c->a[k*c->ncols + j];
			}
		}
	}
}

// Calculate QR decomposition of *transpose* of the matrix
void matrix_qr(struct matrix *mat)
{
	double lwork;
	lapack_int info;

	mat->tau = realloc(mat->tau, min(mat->n, mat->m)*sizeof(double));
	mat->pvt = realloc(mat->pvt, mat->n*sizeof(lapack_int));

	for (int i = 0; i < mat->n; i++) {
		mat->pvt[i] = 0;
	}

	// Query for ideal work size
	info = LAPACKE_dgeqp3_work(LAPACK_COL_MAJOR, mat->m, mat->n, mat->a, mat->ncols, mat->pvt, mat->tau, &lwork, -1);
	assert(info == 0);

	workspace_ensure((lapack_int)lwork);

	// Determine pivoted QR decomposition
	info = LAPACKE_dgeqp3_work(LAPACK_COL_MAJOR, mat->m, mat->n, mat->a, mat->ncols, mat->pvt, mat->tau, matrix_workspace, (lapack_int)lwork);
	assert(info == 0);
}

// Determine rightmost vectors of the q matrix stored in the QR decomposition
// of mat, which are null vectors if mat has a nullspace. If mat is not a QR
// decomposition, one is constructed. The number of vectors returned is defined
// by the number of columns of q.
void matrix_qr_null(struct matrix *mat, struct matrix *q)
{
	assert(q->n == mat->m);
	assert(q->m <= mat->m);

	if (!mat->tau) {
		matrix_qr(mat);
	}

	double lwork;
	lapack_int info;

	// Construct identity matrix in q, selecting the last columns from the decomposition
	for (int i = 0; i < q->n; i++) {
		for (int j = 0; j < q->m; j++) {
			q->a[i*q->ncols + j] = (q->n-i == q->m-j ? 1. : 0.);
		}
	}

	// Query for ideal work size
	info = LAPACKE_dormqr_work(LAPACK_COL_MAJOR, 'R', 'T', q->m, q->n, min(mat->n, mat->m), mat->a, mat->ncols, mat->tau, q->a, q->ncols, &lwork, -1);
	assert(info == 0);

	workspace_ensure((lapack_int)lwork);

	// Obtain the last columns by matrix multiplication with Q
	info = LAPACKE_dormqr_work(LAPACK_COL_MAJOR, 'R', 'T', q->m, q->n, min(mat->n, mat->m), mat->a, mat->ncols, mat->tau, q->a, q->ncols, matrix_workspace, (lapack_int)lwork);
	assert(info == 0);
}

// Determine LU decomposition of A^T
void matrix_lu(struct matrix *mat)
{
	mat->pvt = realloc(mat->pvt, mat->n*sizeof(lapack_int));
	free(mat->tau);
	mat->tau = NULL;

	lapack_int info = LAPACKE_dgetrf_work(LAPACK_COL_MAJOR, mat->m, mat->n, mat->a, mat->ncols, mat->pvt);
	assert(info == 0);
}

// Use LU decomposition of A^T to solve Ax = b for vector b
void matrix_lu_solve(struct matrix *mat, struct matrix *b)
{
	assert(mat->n == mat->m);
	assert(b->n == mat->m);
	assert(b->m == 1);

	matrix_shrink_to_fit(b);

	lapack_int info = LAPACKE_dgetrs_work(LAPACK_COL_MAJOR, 'T', mat->n, 1, mat->a, mat->ncols, mat->pvt, b->a, b->n);
	assert(info == 0);
}
