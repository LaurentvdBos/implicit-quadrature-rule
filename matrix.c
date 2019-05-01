#include "matrix.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

struct matrix *matrix_malloc(int n, int m)
{
	struct matrix *mat = malloc(sizeof(struct matrix));
	mat->n = n;
	mat->m = m;
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
			fprintf(f, fmt, mat->a[i + j*mat->n]);
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
			mat->a[i + j*mat->n] = 0;

			for (int k = 0; k < c->n; k++) {
				mat->a[i + j*mat->n] += b->a[i + k*b->n]*c->a[k + j*c->n];
			}
		}
	}
}

// mat = b^t*c^t = (c*b)^t
void matrix_multt(struct matrix *mat, const struct matrix *b, const struct matrix *c)
{
	assert(mat->n == b->m);
	assert(mat->m == c->n);
	assert(b->n == c->m);

	for (int i = 0; i < mat->m; i++) {
		for (int j = 0; j < mat->n; j++) {
			mat->a[j + i*mat->n] = 0;

			for (int k = 0; k < c->n; k++) {
				mat->a[j + i*mat->n] += c->a[i + k*b->n]*b->a[k + j*c->n];
			}
		}
	}
}

