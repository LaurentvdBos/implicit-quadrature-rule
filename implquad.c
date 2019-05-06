#include "matrix.h"
#include "total_sequence.h"
#include "isort.h"
#include "trie.h"
#include "stack.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include <unistd.h>

int d = 0, n = 0, m = 0;
struct total_sequence *ts;

extern double *matrix_workspace;

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

void implremovals(int *ybest, struct matrix *N, struct matrix *w)
{
	int nz = N->m;
	int sz = N->n;

	// Measure for best removal: removes the largest number of nodes right from m
	int best = 0;

	// Construct initial y to bootstrap procedure
	int y[nz];

	struct matrix *NN = matrix_malloc(N->n, N->m);
	matrix_copy(NN, N);

	struct matrix *ww = matrix_malloc(w->n, w->m);
	matrix_copy(ww, w);
	
	for (int j = 0; j < nz; j++) {
		// Determine alpha and k0
		double alpha = INFINITY; int k0;
		for (int i = 0; i < sz; i++) {
			assert(ww->a[i*ww->lda] >= 0.);

			if (NN->a[i*NN->lda + j] > 0. && ww->a[i*w->lda] / NN->a[i*NN->lda + j] < alpha) {
				alpha = ww->a[i*w->lda] / NN->a[i*NN->lda + j];
				k0 = i;
			}
		}

		y[j] = k0;
		if (k0 >= m) {
			best++;
		}

		// Scale j-th column of NN
		double pivot = NN->a[k0*NN->lda + j];
		double wk0 = ww->a[k0];
		for (int i = 0; i < sz; i++) {
			NN->a[i*NN->lda + j] /= pivot;
			ww->a[i*ww->lda] -= wk0*NN->a[i*NN->lda + j];
		}

		// Update all other null vectors such that k0-th element is zero
		for (int jj = j+1; jj < nz; jj++) {
			double pivot = NN->a[k0*NN->lda + jj];
			for (int i = 0; i < sz; i++) {
				NN->a[i*NN->lda + jj] -= pivot*NN->a[i*NN->lda + j];
			}
		}
	}
	matrix_free(NN);

	isort(y, nz);
	memcpy(ybest, y, nz*sizeof(int));

	struct trie *trie = trie_alloc();
	trie_add(trie, y, nz);
	struct stack *todo = NULL;
	todo = stack_push(todo, best, y, nz);

	// Matrix used to compute all alpha's, it is basically an LU decomposition
	struct matrix *lu = matrix_malloc(nz, nz);

	// Matrix used to store the rhs
	struct matrix *rhs = matrix_malloc(nz, 1);

	// Matrix used to store the null vector under consideration
	struct matrix *c = matrix_malloc(sz, 1);

	while (todo != NULL) {
		memcpy(y, todo->y, nz*sizeof(int));
		int val = todo->val;
		todo = stack_pop(todo);

		// Extract rows y from null space and put as *columns* in lu
		for (int i = 0; i < nz; i++) {
			for (int j = 0; j < nz; j++) {
				lu->a[i*lu->lda + j] = N->a[y[i]*N->lda + j];
			}
			rhs->a[i*rhs->lda] = w->a[y[i]*w->lda];
		}

		// Compute LU
		matrix_lu(lu);

		// Solve for rhs
		matrix_lu_solve(lu, rhs);
		
		// Update ww = w - N*rhs
		matrix_mul(ww, N, rhs);
		for (int i = 0; i < sz; i++) {
			ww->a[i*ww->lda] = w->a[i*w->lda] - ww->a[i*ww->lda];
		}
		for (int i = 0; i < nz; i++) {
			ww->a[y[i]*ww->lda] = 0.;
		}

		// Main loop to explore all other neighboring points on simplex
		for (int i = 0; i < nz; i++) {
			// Prepare rhs
			for (int j = 0; j < nz; j++) {
				rhs->a[j*rhs->lda] = N->a[y[j]*N->lda + i];
			}
			rhs->a[i*rhs->lda] = 0.;

			// Solve for rhs
			matrix_lu_solve(lu, rhs);

			// Determine the null vector
			matrix_mul(c, N, rhs);
			for (int j = 0; j < sz; j++) {
				c->a[j*c->lda] = N->a[j*N->lda + i] - c->a[j*c->lda];
			}
			for (int j = 0; j < nz; j++) {
				if (i == j) {
					continue;
				}
				c->a[y[j]*c->lda] = 0.;
			}

			// Determine kmin, alphamin, kmax, alphamax
			int k1 = -1, k2 = -1;
			double alpha1 = INFINITY, alpha2 = -INFINITY;
			for (int j = 0; j < sz; j++) {
				assert(ww->a[j*ww->lda] >= 0.);

				if (c->a[j*c->lda] > 0. && ww->a[j*ww->lda] / c->a[j*c->lda] < alpha1) {
					alpha1 = ww->a[j*ww->lda] / c->a[j*c->lda];
					k1 = j;
				}
				if (c->a[j*c->lda] < 0. && ww->a[j*ww->lda] / c->a[j*c->lda] > alpha2) {
					alpha2 = ww->a[j*ww->lda] / c->a[j*c->lda];
					k2 = j;
				}
			}
			assert(k1 > -1 && k2 > -1);
			assert(alpha2 <= alpha1);
			assert(0 <= alpha1);
			assert(alpha2 <= 0);
			assert(k1 == y[i] || k2 == y[i]);

			// Determine k0 and new number of nodes that can be removed
			int k0 = (k1 == y[i] ? k2 : k1);
			int newval = val;
			if (k0 >= m && y[i] < m) {
				newval++;
			}
			if (k0 < m && y[i] >= m) {
				newval--;
			}

			// Process y
			int yt[nz];
			memcpy(yt, y, nz*sizeof(int));
			yt[i] = k0;
			isort(yt, nz);

			// Check if we found a better value
			if (newval > best) {
				best = newval;
				memcpy(ybest, yt, nz*sizeof(int));
			}

			if (!trie_contains(trie, yt, nz)) {
				trie_add(trie, yt, nz);
				todo = stack_push(todo, newval, yt, nz);
			}
		}
	}

	matrix_free(c);
	matrix_free(ww);
	matrix_free(rhs);
	matrix_free(lu);
	trie_free(trie);
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
				break;
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

	// Some helper vectors
	struct matrix *lu = matrix_malloc(0, 0);
	struct matrix *rhs = matrix_malloc(0, 0);
	struct matrix *dvec = matrix_malloc(0, 0);

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
			int nz = v->m - v->n;

			// Determine null vector
			matrix_resize(c, q+1, nz);
			matrix_resize(vq, v->n, v->m);
			matrix_copy(vq, v);
			matrix_qr_null(vq, c);

			// Allocate space for best removal on stack
			int yhat[nz];
			implremovals(yhat, c, w);

			// Cool apply that
			matrix_resize(lu, nz, nz);
			matrix_resize(rhs, nz, 1);
			matrix_resize(dvec, q+1, 1);
			for (int i = 0; i < v->m-v->n; i++) {
				for (int j = 0; j < nz; j++) {
					lu->a[i*lu->lda + j] = c->a[yhat[i]*c->lda + j];
				}
				rhs->a[i*rhs->lda] = w->a[yhat[i]*w->lda];
			}
			matrix_lu(lu);
			matrix_lu_solve(lu, rhs);
			
			// Update w = w - N*rhs
			matrix_mul(dvec, c, rhs);
			for (int i = 0; i < q+1; i++) {
				w->a[i*w->lda] -= dvec->a[i*dvec->lda];
			}

			for (int k0 = nz-1; k0 >= 0; k0--) {
				if (yhat[k0] >= m) {
					w->a[yhat[k0]*w->lda] = w->a[q*w->lda];
					for (int j = 0; j < d; j++) {
						x->a[yhat[k0]*x->lda + j] = x->a[q*x->lda + j];
					}
					for (int i = 0; i < n; i++) {
						v->a[i*v->lda + yhat[k0]] = v->a[i*v->lda + q];
					}

					q--;
				} else {
					w->a[yhat[k0]*w->lda] = 0.;
				}
			}

			// Resize
			matrix_resize(x, q+1, d);
			matrix_resize(w, q+1, 1);
			matrix_resize(v, n, q+1);
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
	matrix_free(lu);
	matrix_free(rhs);
	matrix_free(dvec);
	matrix_free(c);
	matrix_free(w);
	matrix_free(x);
	matrix_free(vq);
	matrix_free(v);
	total_sequence_free(ts);

	free(matrix_workspace);

	return EXIT_SUCCESS;
}
