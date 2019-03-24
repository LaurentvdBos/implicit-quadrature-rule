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

// Constructs a Householder transformation of the bottom right part of the
// matrix mat. The bottom right starts at column k and row k.
void matrix_householder(struct matrix *mat, struct matrix *q, int k)
{
	assert(q->n == mat->n);
	assert(q->m == mat->n);

	double alpha = 0;
	for (int i = k; i < mat->n; i++) {
		alpha += mat->a[i + k*mat->n]*mat->a[i + k*mat->n];
	}

	double norm = 2.*alpha;

	alpha = sqrt(alpha);
	if (mat->a[k + k*mat->n] < 0) {
		alpha = -alpha;
	}
	norm -= 2.*alpha*mat->a[k + k*mat->n];

	// Construct Q
	for (int i = 0; i < q->n; i++) {
		for (int j = 0; j < q->m; j++) {
			q->a[i + j*q->n] = (i == j ? 1 : 0);
			if (i >= k && j >= k) {
				q->a[i + j*q->n] -= 2. * (i == k ? mat->a[i + k*mat->n]-alpha : mat->a[i + k*mat->n]) * (j == k ? mat->a[j + k*mat->n]-alpha : mat->a[j + k*mat->n]) / norm;
			}
		}
	}
}

// Determine QR, assumes column-major ordering
void matrix_qr(struct matrix *mat, struct matrix *q, struct matrix *r)
{
	assert(q->n == mat->n);
	assert(q->m == mat->n);
	assert(r->n == mat->n);
	assert(r->m == mat->m);

	int n = mat->n, m = mat->m;

	// Copy mat in r
	for (int i = 0; i < n*m; i++) {
		r->a[i] = mat->a[i];
	}

	// Make q identity
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			q->a[i + j*n] = (i == j ? 1 : 0);
		}
	}

	// Store Householder transformations in this matrix
	struct matrix *h = matrix_malloc(n, n);
	struct matrix *tmp = matrix_malloc(n, n);

	// Construct actual qr decomposition
	for (int i = 0; i < min(n-1, m); i++) {
		matrix_householder(r, h, i);
		matrix_mul(tmp, h, q);
		matrix_swap(tmp, q);
		matrix_mul(r, q, mat);
	}

	// Transpose Q
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			double tmp = q->a[i + j*n];
			q->a[i + j*n] = q->a[j + i*n];
			q->a[j + i*n] = tmp;
		}
	}

	matrix_free(h);
	matrix_free(tmp);
}
