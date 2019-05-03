#include "matrix.h"
#include "total_sequence.h"

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include <unistd.h>

int d = 0, n = 0, m = 0;
struct total_sequence *ts;

void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s [-d dim] [-n number of nodes]\n\n", myname);
	fprintf(stderr, "This program start reading samples from standard input until eof and\n"
	                "prints the quadrature consisting on n nodes to standard output.\n");
}


void vdm_col(struct matrix *v, double *y)
{
	int j = v->m-1;

	total_sequence_set_sum(ts, 0);

	for (int i = 0; i < n; i++) {
		v->a[i*v->lda + j] = 1.;

		int l = 0;
		for (struct node *n = ts->root; n; n = n->next) {
			v->a[i*v->lda + j] *= pow(y[l++], (double)n->val);
		}

		total_sequence_next(ts);
	}
}

int main(int argc, char **argv)
{
	if (argc < 1) {
		usage("[this-program]");
		return EXIT_FAILURE;
	} else if (argc == 1) {
		usage(argv[0]);
		return EXIT_FAILURE;
	}

	int opt;
	while ((opt = getopt(argc, argv, "d:n:m:h?")) != -1) {
		switch (opt) {
			case 'd':
				d = atoi(optarg);
				break;
			case 'n':
				n = atoi(optarg);
				break;
			case 'm':
				m = atoi(optarg);
			case 'h':
			case '?':
			default:
				usage(argv[0]);
				return EXIT_FAILURE;
		}
	}

	if (d <= 0) {
		fprintf(stderr, "%s: invalid value for 'd': %d\n", argv[0], d);
		usage(argv[0]);
		return EXIT_FAILURE;
	}
	if (n <= 0) {
		fprintf(stderr, "%s: invalid value for 'n': %d\n", argv[0], n);
		usage(argv[0]);
		return EXIT_FAILURE;
	}

	// Initialize total sequence
	ts = total_sequence_alloc(d);

	// Vandermonde-matrix
	struct matrix *v = matrix_malloc(n, 0);
	struct matrix *vq = matrix_malloc(n, 0);

	// List of nodes
	struct matrix *x = matrix_malloc(0, d);
	
	// List of weights
	struct matrix *w = matrix_malloc(0, 1);

	// Null space
	struct matrix *c = matrix_malloc(0, 0);

	// Current sample
	double *y = malloc(sizeof(double)*d);

	// Enter main loop; q denotes number of nodes, k number of samples
	for (int q = 0, k = 0; ; q++, k++) {
		// Read a node
		for (int i = 0; i < d; i++) {
			if (scanf(" %lf", &y[i]) == EOF) {
				goto out;
			}
		}

		// Copy node into x
		matrix_resize(x, q+1, d);
		for (int i = 0; i < d; i++) {
			x->a[q*x->lda + i] = y[i];
		}

		// Adjust weights
		for (int i = 0; i < q; i++) {
			w->a[i*w->lda] = w->a[i*w->lda] * (double)k / (double)(k+1);
		}
		matrix_resize(w, q+1, 1);
		w->a[q*w->lda] = 1. / (double)(k+1);

		// Add column to v
		matrix_resize(v, n, q+1);
		vdm_col(v, y);

		if (n < q+1) {
			// Determine null vector
			matrix_resize(c, q+1, 1);
			matrix_resize(vq, v->n, v->m);
			matrix_copy(vq, v);
			matrix_null(vq, c);

			// Determine alpha and k0
			double alpha = INFINITY; int k0;
			for (int i = 0; i < q+1; i++) {
				if (c->a[i*c->lda] > 0. && w->a[i*w->lda] / c->a[i*c->lda] < alpha) {
					alpha = w->a[i*w->lda] / c->a[i*c->lda];
					k0 = i;
				}
			}

			// Update weights
			matrix_add(w, 1., w, -alpha, c);

			// Replace with last row
			w->a[k0*w->lda] = w->a[q*w->lda];
			for (int j = 0; j < d; j++) {
				x->a[k0*x->lda + j] = x->a[q*x->lda + j];
			}
			for (int i = 0; i < n; i++) {
				v->a[i*v->lda + k0] = v->a[i*v->lda + q];
			}

			// Resize
			matrix_resize(x, q, d);
			matrix_resize(w, q, 1);
			matrix_resize(v, n, q);

			q--;
		}
	}
out:
	
	// Print list of nodes; these outputs should be cleared up in the future
	printf("Nodes:\n");
	matrix_fprintf(stdout, x, "%.10f");

	printf("\nWeights:\n");
	matrix_fprintf(stdout, w, "%.10f");

	// Clean up
	free(y);
	matrix_free(c);
	matrix_free(w);
	matrix_free(x);
	matrix_free(vq);
	matrix_free(v);
	total_sequence_free(ts);

	return EXIT_SUCCESS;
}
