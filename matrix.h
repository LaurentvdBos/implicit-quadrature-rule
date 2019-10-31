#ifndef MATRIX_H
#define MATRIX_H

struct matrix
{
	// Number of rows/columns in the matrix
	int n, m;

	// Number of rows/columns in the data structure, which can be larger than n/m
	int nrows, ncols;

	// Array of matrix coefficients in row-major order (nrows*ncols)
	double *a;

	// The LU decomposition uses partial pivoting, stored in the following
	// array. It is also used as an indicator that the matrix is an
	// LU-decomposition
	int *pvt;
};

struct matrix *matrix_malloc(const int n, const int m);
void matrix_free(struct matrix *mat);
void matrix_copy(struct matrix *mat, const struct matrix *b);
void matrix_resize(struct matrix *mat, const int n, const int m);
void matrix_mul(struct matrix *mat, const struct matrix *b, const struct matrix *c);

#endif
