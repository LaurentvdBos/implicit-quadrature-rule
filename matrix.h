#ifndef MATRIX_H
#define MATRIX_H

#include <lapacke.h>

struct matrix
{
	// Number of rows/columns in the matrix
	int n, m;

	// Number of rows/columns in the data structure, which can be larger than n/m
	int nrows, ncols;

	// Array of matrix coefficients in row-major order (nrows*ncols)
	double *a;

	// These are used by LAPACK to store a decomposition
	// QR used both tau and pvt; LU only uses pvt
	double *tau;
	lapack_int *pvt;
};

struct matrix *matrix_malloc(int n, int m);
void matrix_free(struct matrix *mat);
void matrix_copy(struct matrix *mat, const struct matrix *b);
void matrix_resize(struct matrix *mat, const int n, const int m);
void matrix_mul(struct matrix *mat, const struct matrix *b, const struct matrix *c);

void matrix_qr(struct matrix *mat);
void matrix_qr_null(struct matrix *mat, struct matrix *q);

void matrix_lu(struct matrix *mat);
void matrix_lu_solve(struct matrix *mat, struct matrix *b);

#endif
