#include "lu.h"
#include "matrix.h"

#include <assert.h>
#include <stdlib.h>
#include <tgmath.h>

// Determine an LU decomposition of a matrix in place. The decomposition is
// determined using Gaussian elimination and partial pivoting (i.e. columnwise
// pivoting)
void matrix_lu(struct matrix *mat)
{
	const int n = mat->n;
	const int m = mat->m;

	int *pvt = realloc(mat->pvt, n*sizeof(int));
	mat->pvt = pvt;
	for (int i = 0; i < n; i++) {
		pvt[i] = i;
	}

	for (int i = 0; i < n; i++) {
		// Find pivot in this column
		double amax = mat->a[pvt[i]*mat->ncols + i];
		int imax = i;
		for (int k = i+1; k < n; k++) {
			double newmax = fabs(mat->a[pvt[k]*mat->ncols + i]);
			if (newmax > amax) {
				amax = newmax; imax = k;
			}
		}

		// Store pivot
		if (imax != i) {
			int k = pvt[i];
			pvt[i] = pvt[imax];
			pvt[imax] = k;
		}

		// Do a Gaussian reduction
		for (int k = i+1; k < n; k++) {
			mat->a[pvt[k]*mat->ncols + i] /= mat->a[pvt[i]*mat->ncols + i];

			for (int j = i+1; j < m; j++) {
				mat->a[pvt[k]*mat->ncols + j] -= mat->a[pvt[k]*mat->ncols + i]*mat->a[pvt[i]*mat->ncols + j];
			}
		}
	}
}

// Solve A*x = b using the LU decomposition of A, so first solve U*y = b and
// then solve L*x = y. These solves can be done in place and the result of the
// operation is stored in b.
void matrix_lu_solve(struct matrix *mat, struct matrix *b)
{
	assert(mat->n == mat->m);
	assert(b->n == mat->m);
	assert(mat->pvt);

	const int n = mat->n;
	const int m = b->m;
	const int *pvt = mat->pvt;

	// Solve L*x = b -> store result immediately in b
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < i; k++) {
			for (int j = 0; j < m; j++) {
				b->a[pvt[i]*b->ncols + j] -= mat->a[pvt[i]*mat->ncols + k] * b->a[pvt[k]*b->ncols + j];
			}
		}
	}

	// Solve U*x = b -> store result immediately in b
	for (int i = n-1; i >= 0; i--) {
		for (int k = i+1; k < n; k++) {
			for (int j = 0; j < m; j++) {
				b->a[pvt[i]*b->ncols + j] -= mat->a[pvt[i]*mat->ncols + k] * b->a[pvt[k]*b->ncols + j];
			}
		}

		for (int j = 0; j < m; j++) {
			b->a[pvt[i]*b->ncols + j] /= mat->a[pvt[i]*mat->ncols + i];
		}
	}

	// Unpivot b: look iteratively where the element has been going to
	for (int i = 0; i < n; i++) {
		int k = pvt[i];
		while (k < i) {
			k = pvt[k];
		}

		for (int j = 0; j < m; j++) {
			double tmp = b->a[i*b->ncols + j];
			b->a[i*b->ncols + j] = b->a[k*b->ncols + j];
			b->a[k*b->ncols + j] = tmp;
		}
	}
}

// Determine null vectors of A using the LU decomposition of A. It is only
// necessary to calculate the null vectors of U, since L*0 = 0. The null
// vectors are determined by solving the linear system A*x = C, where A are the
// leftmost columns of U forming a square and C are the rightmost columns of U
// forming the rest. Everything is done in place.
void matrix_lu_null(struct matrix *mat, struct matrix *c)
{
	assert(c->n == mat->m);
	assert(c->m <= mat->m - mat->n);
	assert(mat->pvt);

	const int n = mat->n;
	const int m = mat->m;
	const int mz = c->m;
	const int *pvt = mat->pvt;

	// Copy rightmost part of U into c; neglecting the bottom part
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < mz; j++) {
			c->a[i*c->ncols + j] = mat->a[pvt[i]*mat->ncols + (m-j-1)];
		}
	}

	// Solve A*x = c, where A is formed by the leftmost n columns of U, and
	// immediately store the result in c.
	for (int i = n-1; i >= 0; i--) {
		for (int k = i+1; k < n; k++) {
			for (int j = 0; j < mz; j++) {
				c->a[i*c->ncols + j] -= mat->a[pvt[i]*mat->ncols + k] * c->a[k*c->ncols + j];
			}
		}
		
		for (int j = 0; j < mz; j++) {
			c->a[i*c->ncols + j] /= mat->a[pvt[i]*mat->ncols + i];
		}
	}

	// Put a reverse negative identity matrix in the bottom of c and zeroes elsewhere
	for (int i = n; i < m; i++) {
		for (int j = 0; j < mz; j++) {
			c->a[i*c->ncols + j] = 0.L;
		}

		int j = m-i-1;
		c->a[i*c->ncols + j] = -1.L;
	}
}
