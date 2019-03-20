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
		mat->a = calloc(sizeof(double), n*m);
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

// Constructs a Householder transformation of the bottom right part of the
// matrix mat. The bottom right starts at column k and row k.
// TODO: You probably do not want to explictly construct q, but only construct
// u and store that in the QR matrix (i.e. the "LAPACK" way)
void matrix_householder(struct matrix *mat, struct matrix *q, int k)
{
	assert(q->n == mat->n);
	assert(q->m == mat->n);

	// alpha = ||col_k||
	double alpha = 0;

	// u = x - alpha e_1
	struct matrix *u = matrix_malloc(mat->n-k, 1);
	for (int i = k; i < mat->n; i++) {
		u->a[i-k] = mat->a[i + k*mat->n];
		alpha += mat->a[i + k*mat->n]*mat->a[i + k*mat->n];
	}

	// norm = ||u||^2
	double norm = 2*alpha;

	alpha = sqrt(alpha);
	if (mat->a[k + k*mat->n] > 0) {
		norm -= 2*alpha*u->a[0];
		u->a[0] -= alpha;
	} else {
		norm += 2*alpha*u->a[0];
		u->a[0] += alpha;
	}

	// Construct Q
	for (int i = 0; i < q->n; i++) {
		for (int j = 0; j < q->m; j++) {
			q->a[i + j*q->n] = (i == j ? 1 : 0);
			if (i >= k && j >= k) {
				q->a[i + j*q->n] -= 2*u->a[i-k]*u->a[j-k]/norm;
			}
		}
	}

	matrix_free(u);
}

// Determine QR, assumes column-major ordering
void matrix_qr(struct matrix *mat, struct matrix *q, struct matrix *r)
{
	for (int i = 0; i < min(mat->n, mat->m-1); i++) {
	}
}
