#include "matrix.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lapacke.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

// LAPACK requires a small bit of memory to calculate QR decompositions. It is
// stored in this global variable and freed in main (see implquad.c).
// workspace_ensure can be used to ensure that the workspace has a certain
// size. You probably want to change this part if you want to use this software
// in a multithreaded environment
double *matrix_workspace = NULL;
lapack_int nwork = 0;

static inline void workspace_ensure(lapack_int n)
{
	if (n > nwork) {
		nwork = n;
		matrix_workspace = realloc(matrix_workspace, nwork*sizeof(double));
	}
}

// malloc a matrix with number of rows n and columns m. It is uninitialized.
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

// Free a matrix; after this function the pointer is invalid
void matrix_free(struct matrix *mat)
{
	if (mat) {
		free(mat->tau);
		free(mat->pvt);
		free(mat->a);
		free(mat);
	}
}

// Copy matrix b to mat. We cannot use memcpy here, since the number of rows or
// columns available in b can be different from the number in mat.
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

// Resize a matrix such that it has a different number of rows and columns.
// Notice that the number of rows or columns a matrix *can* store is different
// than the number it actually has. This function only expands the storage, or
// reuses it if it is already available.
void matrix_resize(struct matrix *mat, const int n, const int m)
{
	if (mat->ncols >= m) {
		if (mat->nrows < n) {
			mat->a = realloc(mat->a, mat->ncols*n*sizeof(double));
			mat->nrows = n;
		}
	} else {
		if (mat->nrows < n) {
			mat->nrows = n;
		}
		mat->a = realloc(mat->a, m*mat->nrows*sizeof(double));
		
		// Copy data to new memory block, keeping the new number of columns into account
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

// Shrink the array of the matrix such that it *exactly* contains the number of
// rows and columns of the matrix. This is necessary for some LAPACK functions.
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

// Basic implementation of matrix multiplication. The matrix mat cannot be
// equal to b or c; b and c can be equal.
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

// Calculate QR decomposition of *transpose* of the matrix. We do not have to
// actively transpose the matrix, since LAPACK is column major and the
// implementation here is row major.
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

// Determine LU decomposition of the *transpose* of A. We do not have to
// calculate the transpose of A, since LAPACK is column major and the
// implementation here is row major.
void matrix_lu(struct matrix *mat)
{
	mat->pvt = realloc(mat->pvt, mat->n*sizeof(lapack_int));
	free(mat->tau);
	mat->tau = NULL;

	// No workspace allocation is necessary for an LU decomposition
	lapack_int info = LAPACKE_dgetrf_work(LAPACK_COL_MAJOR, mat->m, mat->n, mat->a, mat->ncols, mat->pvt);
	assert(info == 0);
}

// Use LU decomposition of *transpose* of A to solve Ax = b for vector b.
// LAPACK can solve A^T x = b, so we ask LAPACK to do the transposition for us.
// We need to ensure that b is compact, i.e. the array of b should not be
// larger than necessary.
void matrix_lu_solve(struct matrix *mat, struct matrix *b)
{
	assert(mat->n == mat->m);
	assert(b->n == mat->m);
	assert(b->m == 1);

	matrix_shrink_to_fit(b);

	lapack_int info = LAPACKE_dgetrs_work(LAPACK_COL_MAJOR, 'T', mat->n, 1, mat->a, mat->ncols, mat->pvt, b->a, b->n);
	assert(info == 0);
}
