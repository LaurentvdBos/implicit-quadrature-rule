#include "tree.h"

#include <stdlib.h>
#include <string.h>

static inline int compar(const int *a, const int *b, const int n)
{
	for (int i = 0; i < n; i++) {
		if (b[i] != a[i]) {
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
	root->subtree = NULL;

	root->a = malloc(n*sizeof(int));
	memcpy(root->a, a, n*sizeof(int));

	return root;
}

void tree_free(struct tree *root)
{
	if (root) {
		tree_free(root->left);
		tree_free(root->right);
		tree_free(root->subtree);

		free(root->a);
		free(root);
	}
}

// WARNING: Only add something to the tree if you're very sure that it is not
// contained in the tree yet. Otherwise you'll get some nasty errors.
void tree_add(struct tree *root, const int *a, const int n)
{
	struct tree *ptr = root;

	int i = 0;

	while (ptr) {
		if (ptr->a[0] == a[i]) {
			i++;

			if (ptr->subtree == NULL) {
				ptr->subtree = tree_malloc(a+i, n-i);
				break;
			} else {
				ptr = ptr->subtree;
			}
		} else if (a[i] < ptr->a[0]) {
			if (ptr->left == NULL) {
				ptr->left = tree_malloc(a+i, n-i);
				break;
			} else {
				ptr = ptr->left;
			}
		} else {
			if (ptr->right == NULL) {
				ptr->right = tree_malloc(a+i, n-i);
				break;
			} else {
				ptr = ptr->right;
			}
		}
	}
}

bool tree_contains(const struct tree *root, const int *a, const int n)
{
	const struct tree *ptr = root;

	int i = 0;

	while (ptr) {
		int c = compar(ptr->a, a+i, n-i);

		if (ptr->a[0] == a[i] && c != 0) {
			i++;

			ptr = ptr->subtree;
		} else if (c < 0) {
			ptr = ptr->left;
		} else if (c > 0) {
			ptr = ptr->right;
		} else {
			return true;
		}
	}

	return false;
}
