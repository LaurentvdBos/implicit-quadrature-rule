#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>

struct matrix
{
	int n, m;
	double *a;
};

struct matrix *matrix_malloc(int n, int m);
void matrix_free(struct matrix *mat);
void matrix_set(struct matrix *mat, const double *a);
void matrix_fprintf(FILE *f, const struct matrix *mat, const char *fmt);
void matrix_mul(struct matrix *mat, const struct matrix *b, const struct matrix *c);
void matrix_multt(struct matrix *mat, const struct matrix *b, const struct matrix *c);
void matrix_householder(struct matrix *mat, struct matrix *q, int k);

#endif
