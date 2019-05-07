#include "matrix.h"
#include "total_sequence.h"
#include "isort.h"
#include "trie.h"
#include "stack.h"

#include <assert.h>
#include <float.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include <lapacke.h>
#include <getopt.h>

// DBL_DECIMAL_DIG is C11 and later only
#ifndef DBL_DECIMAL_DIG
  #ifdef DECIMAL_DIG
    #define DBL_DECIMAL_DIG DECIMAL_DIG
  #else
    #define DBL_DECIMAL_DIG 16
  #endif
#endif

// Size of the progressbar
#define PROGRESS_BAR 80

int d = 0, n = 0, m = 0, K = 0;
bool print_nodes = true, print_weights = true, print_index = false, pretty_print = true;
struct total_sequence *ts;

extern double *matrix_workspace;

void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s [-d dim] [-n number of nodes]\n\n", myname);
	fprintf(stderr, "This program starts reading samples from standard input until eof and\n"
	                "prints the quadrature consisting of n nodes to standard output.\n\n");

	fprintf(stderr, "Options with a + expect a number. Order of options is not important.\n\n");
	fprintf(stderr, "Compulsory options:\n");
	fprintf(stderr, "  -d+ Dimension of the sample space; the samples can be provided unstructured\n");
	fprintf(stderr, "  -n+ Number of nodes in the obtained quadrature rule\n\n");

	fprintf(stderr, "Optional options (defaults to -xwqm0):\n");
	fprintf(stderr, "  -m+ Number of nodes that should be preserved in the rule. Providing a\n"
	                "      non-zero integer yields a quadrature rule that at least contains\n"
	                "      the first m samples.\n");
	fprintf(stderr, "  -y+ Total number of samples provided. If provided, prints a progress\n"
	                "      bar to stderr.\n\n");
	fprintf(stderr, "  -x  Print nodes\n");
	fprintf(stderr, "  -w  Print weights\n");
	fprintf(stderr, "  -i  Print indices of samples used. List of samples is zero-indexed.\n");
	fprintf(stderr, "  -q  Print nodes and weights seperately (otherwise as one big matrix)\n\n");

	lapack_int major, minor, patch;
	LAPACKE_ilaver(&major, &minor, &patch);
	fprintf(stderr, "This %s uses LAPACK %d.%d.%d\n", myname, major, minor, patch);
}

void vdm_col(struct matrix *v, double *y)
{
	int j = v->m-1;

	total_sequence_set_sum(ts, 0);

	for (int i = 0; i < n; i++) {
		v->a[i*v->ncols + j] = 1.;

		int l = 0;
		for (struct node *n = ts->root; n; n = n->next) {
			v->a[i*v->ncols + j] *= pow(y[l++], (double)n->val);
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

	// Arrays to keep track of indices that will be made zero
	int *y = malloc(nz*sizeof(int));
	int *yt = malloc(nz*sizeof(int));

	struct matrix *NN = matrix_malloc(N->n, N->m);
	matrix_copy(NN, N);

	struct matrix *ww = matrix_malloc(w->n, w->m);
	matrix_copy(ww, w);
	
	for (int j = 0; j < nz; j++) {
		// Determine alpha and k0
		double alpha = INFINITY; int k0 = -1;
		for (int i = 0; i < sz; i++) {
			assert(ww->a[i*ww->ncols] >= 0.);

			if (NN->a[i*NN->ncols + j] > 0. && ww->a[i*w->ncols] / NN->a[i*NN->ncols + j] < alpha) {
				alpha = ww->a[i*w->ncols] / NN->a[i*NN->ncols + j];
				k0 = i;
			}
		}
		assert(k0 > -1);

		y[j] = k0;
		if (k0 >= m) {
			best++;
		}

		// Scale j-th column of NN
		double pivot = NN->a[k0*NN->ncols + j];
		double wk0 = ww->a[k0];
		for (int i = 0; i < sz; i++) {
			NN->a[i*NN->ncols + j] /= pivot;
			ww->a[i*ww->ncols] -= wk0*NN->a[i*NN->ncols + j];
		}

		// Update all other null vectors such that k0-th element is zero
		for (int jj = j+1; jj < nz; jj++) {
			double pivot = NN->a[k0*NN->ncols + jj];
			for (int i = 0; i < sz; i++) {
				NN->a[i*NN->ncols + jj] -= pivot*NN->a[i*NN->ncols + j];
			}
		}
	}
	matrix_free(NN);

	isort(y, nz);
	memcpy(ybest, y, nz*sizeof(int));

	struct trie *trie = trie_malloc();
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
				lu->a[i*lu->ncols + j] = N->a[y[i]*N->ncols + j];
			}
			rhs->a[i*rhs->ncols] = w->a[y[i]*w->ncols];
		}

		// Compute LU
		matrix_lu(lu);

		// Solve for rhs
		matrix_lu_solve(lu, rhs);
		
		// Update ww = w - N*rhs
		matrix_mul(ww, N, rhs);
		for (int i = 0; i < sz; i++) {
			ww->a[i*ww->ncols] = w->a[i*w->ncols] - ww->a[i*ww->ncols];
		}
		for (int i = 0; i < nz; i++) {
			ww->a[y[i]*ww->ncols] = 0.;
		}

		// Main loop to explore all other neighboring points on simplex
		for (int i = 0; i < nz; i++) {
			// Prepare rhs
			for (int j = 0; j < nz; j++) {
				rhs->a[j*rhs->ncols] = N->a[y[j]*N->ncols + i];
			}
			rhs->a[i*rhs->ncols] = 0.;

			// Solve for rhs
			matrix_lu_solve(lu, rhs);

			// Determine the null vector
			matrix_mul(c, N, rhs);
			for (int j = 0; j < sz; j++) {
				c->a[j*c->ncols] = N->a[j*N->ncols + i] - c->a[j*c->ncols];
			}
			for (int j = 0; j < nz; j++) {
				if (i == j) {
					continue;
				}
				c->a[y[j]*c->ncols] = 0.;
			}

			// Determine kmin, alphamin, kmax, alphamax
			int k1 = -1, k2 = -1;
			double alpha1 = INFINITY, alpha2 = -INFINITY;
			for (int j = 0; j < sz; j++) {
				assert(ww->a[j*ww->ncols] >= 0.);

				if (c->a[j*c->ncols] > 0. && ww->a[j*ww->ncols] / c->a[j*c->ncols] < alpha1) {
					alpha1 = ww->a[j*ww->ncols] / c->a[j*c->ncols];
					k1 = j;
				}
				if (c->a[j*c->ncols] < 0. && ww->a[j*ww->ncols] / c->a[j*c->ncols] > alpha2) {
					alpha2 = ww->a[j*ww->ncols] / c->a[j*c->ncols];
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

	free(y);
	free(yt);
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
	extern char *optarg;
	while ((opt = getopt(argc, argv, "d:n:m:P:y:xwiqh?")) != -1) {
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
			case 'P':
			case 'y':
				K = atoi(optarg);
				break;
			case 'x':
				print_nodes = !print_nodes;
				break;
			case 'w':
				print_weights = !print_weights;
				break;
			case 'i':
				print_index = !print_index;
				break;
			case 'q':
				pretty_print = !pretty_print;
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
	if (!print_nodes && !print_weights && !print_index) {
		fprintf(stderr, "Nothing to do.\n");
		return EXIT_SUCCESS;
	}

	// Initialize total sequence
	ts = total_sequence_malloc(d);

	// Vandermonde-matrix
	struct matrix *v = matrix_malloc(n, 0);
	struct matrix *vq = matrix_malloc(n, 0);

	// List of nodes
	struct matrix *x = matrix_malloc(0, d);
	
	// List of weights
	struct matrix *w = matrix_malloc(0, 1);

	// Index of nodes
	int *index = NULL;

	// Null space
	struct matrix *c = matrix_malloc(0, 0);

	// Some helper vectors
	struct matrix *lu = matrix_malloc(0, 0);
	struct matrix *rhs = matrix_malloc(0, 0);
	struct matrix *dvec = matrix_malloc(0, 0);

	// Storage for the best removal
	int *yhat = NULL;

	// Current sample
	double *y = malloc(d*sizeof(double));

	// Keep track of progress
	int p = 0;

	// Enter main loop; q denotes number of nodes, k number of samples
	for (int q = 0, k = 0; ; q++, k++) {
		if (K > 0 && k*PROGRESS_BAR > p*K) {
			p++;
			fprintf(stderr, ".");
		}
				
		// Read a node
		for (int i = 0; i < d; i++) {
			if (scanf(" %lf", &y[i]) == EOF) {
				goto out;
			}
		}

		// Copy node into x
		matrix_resize(x, q+1, d);
		for (int i = 0; i < d; i++) {
			x->a[q*x->ncols + i] = y[i];
		}

		// Adjust weights
		for (int i = 0; i < q; i++) {
			w->a[i*w->ncols] = w->a[i*w->ncols] * (double)k / (double)(k+1);
		}
		matrix_resize(w, q+1, 1);
		w->a[q*w->ncols] = 1. / (double)(k+1);

		// Add index
		index = realloc(index, (q+1)*sizeof(int));
		index[q] = k;

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

			// Allocate space for the best removal
			yhat = realloc(yhat, nz*sizeof(int));
			implremovals(yhat, c, w);

			// Cool apply that
			matrix_resize(lu, nz, nz);
			matrix_resize(rhs, nz, 1);
			matrix_resize(dvec, q+1, 1);
			for (int i = 0; i < v->m-v->n; i++) {
				for (int j = 0; j < nz; j++) {
					lu->a[i*lu->ncols + j] = c->a[yhat[i]*c->ncols + j];
				}
				rhs->a[i*rhs->ncols] = w->a[yhat[i]*w->ncols];
			}
			matrix_lu(lu);
			matrix_lu_solve(lu, rhs);
			
			// Update w = w - N*rhs
			matrix_mul(dvec, c, rhs);
			for (int i = 0; i < q+1; i++) {
				w->a[i*w->ncols] -= dvec->a[i*dvec->ncols];
			}

			for (int k0 = nz-1; k0 >= 0; k0--) {
				if (yhat[k0] >= m) {
					w->a[yhat[k0]*w->ncols] = w->a[q*w->ncols];
					index[yhat[k0]] = index[q];
					for (int j = 0; j < d; j++) {
						x->a[yhat[k0]*x->ncols + j] = x->a[q*x->ncols + j];
					}
					for (int i = 0; i < n; i++) {
						v->a[i*v->ncols + yhat[k0]] = v->a[i*v->ncols + q];
					}

					q--;
				} else {
					w->a[yhat[k0]*w->ncols] = 0.;
				}
			}

			// Resize
			matrix_resize(x, q+1, d);
			matrix_resize(w, q+1, 1);
			matrix_resize(v, n, q+1);
		}
	}
out:

	if (K > 0) {
		fprintf(stderr, "\n");
	}
	
	if (pretty_print) {
		if (print_nodes) {
			printf("Nodes:\n");
			for (int i = 0; i < x->n; i++) {
				for (int j = 0; j < x->m; j++) {
					printf("%.*e ", DBL_DECIMAL_DIG, x->a[i*x->ncols + j]);
				}
				printf("\n");
			}
			if (print_weights || print_index) {
				printf("\n");
			}
		}

		if (print_weights) {
			printf("Weights:\n");
			for (int i = 0; i < w->n; i++) {
				printf("%.*e\n", DBL_DECIMAL_DIG, w->a[i*w->ncols]);
			}
			if (print_index) {
				printf("\n");
			}
		}

		if (print_index) {
			printf("Index:\n");
			for (int i = 0; i < w->n; i++) {
				printf("%d\n", index[i]);
			}
		}
	} else {
		for (int i = 0; i < x->n; i++) {
			if (print_index) {
				printf("%d ", index[i]);
			}
			if (print_nodes) {
				for (int j = 0; j < x->m; j++) {
					printf("%.*e ", DBL_DECIMAL_DIG, x->a[i*x->ncols + j]);
				}
			}
			if (print_weights) {
				printf("%.*e", DBL_DECIMAL_DIG, w->a[i*w->ncols]);
			}
			printf("\n");
		}
	}

	// Clean up
	free(y);
	free(yhat);
	free(index);
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
