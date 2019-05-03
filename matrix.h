#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>

struct matrix
{
	// Size of the matrix
	int n, m;

	// Leading dimension, i.e. number of columns in the data structure
	int lda;

	// Array of matrix coefficients in row-major order (n*lda)
	double *a;
};

struct matrix *matrix_malloc(int n, int m);
void matrix_free(struct matrix *mat);
void matrix_set(struct matrix *mat, const double *a);
void matrix_copy(struct matrix *mat, const struct matrix *b);
void matrix_resize(struct matrix *mat, const int n, const int m);
void matrix_fprintf(FILE *f, const struct matrix *mat, const char *fmt);
void matrix_add(struct matrix *mat, const double alpha, const struct matrix *b, const double beta, const struct matrix *c);
void matrix_mul(struct matrix *mat, const struct matrix *b, const struct matrix *c);
void matrix_null(struct matrix *mat, struct matrix *q);

#endif
