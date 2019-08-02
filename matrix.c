#include "matrix.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

// malloc a matrix with number of rows n and columns m. It is uninitialized.
struct matrix *matrix_malloc(const int n, const int m)
{
	struct matrix *mat = malloc(sizeof(struct matrix));
	mat->n = n;
	mat->m = m;
	mat->nrows = n;
	mat->ncols = m;
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

	if (b->pvt) {
		mat->pvt = realloc(mat->pvt, mat->n*sizeof(int));
		memcpy(mat->pvt, b->pvt, mat->n*sizeof(int));
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
	free(mat->pvt);
	mat->pvt = NULL;
}

// Shrink the array of the matrix such that it *exactly* contains the number of
// rows and columns of the matrix.
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
