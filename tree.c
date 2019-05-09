#include "tree.h"

#include <stdlib.h>
#include <string.h>

static inline int compar(const int *a, const int *b, const int n)
{
	for (int i = 0; i < n; i++) {
		if (a[i] != b[i]) {
			return (b[i] - a[i]);
		}
	}

	return 0;
}

struct tree *tree_malloc(const int *a, const int n)
{
	struct tree *root = malloc(sizeof(struct tree));

	root->left = NULL;
	root->right = NULL;

	root->a = malloc(n*sizeof(int));
	memcpy(root->a, a, n*sizeof(int));

	return root;
}

void tree_free(struct tree *root)
{
	if (root) {
		tree_free(root->left);
		tree_free(root->right);

		free(root->a);
		free(root);
	}
}

void tree_add(struct tree *root, int *a, int n)
{
	struct tree *ptr = root;

	while (ptr) {
		if (compar(ptr->a, a, n) > 0) {
			if (ptr->left == NULL) {
				ptr->left = tree_malloc(a, n);
				break;
			} else {
				ptr = ptr->left;
			}
		} else {
			if (ptr->right == NULL) {
				ptr->right = tree_malloc(a, n);
				break;
			} else {
				ptr = ptr->right;
			}
		}
	}
}

bool tree_contains(struct tree *root, int *a, int n)
{
	struct tree *ptr = root;

	while (ptr) {
		int c = compar(ptr->a, a, n);
	
		if (c > 0) {
			ptr = ptr->left;
		} else if (c < 0) {
			ptr = ptr->right;
		} else {
			return true;
		}
	}

	return false;
}
