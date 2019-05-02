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
void matrix_fprintf(FILE *f, const struct matrix *mat, const char *fmt);
void matrix_mul(struct matrix *mat, const struct matrix *b, const struct matrix *c);
void matrix_null(struct matrix *mat, struct matrix *q);

#endif
