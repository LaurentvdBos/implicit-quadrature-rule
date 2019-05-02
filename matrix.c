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
			fprintf(f, fmt, mat->a[i*mat->m + j]);
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
			mat->a[i*mat->m + j] = 0;

			for (int k = 0; k < c->n; k++) {
				mat->a[i*mat->m + j] += b->a[i*b->m + k]*c->a[k*c->m + j];
			}
		}
	}
}
