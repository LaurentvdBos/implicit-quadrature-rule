#include "tree.h"

#include <stdlib.h>
#include <string.h>

// Compare two arrays; returns 0 if b == a, a negative number if b < a, and a
// positive number if b > a.
static inline int compar(const int *a, const int *b, const int n)
{
	for (int i = 0; i < n; i++) {
		if (b[i] != a[i]) {
			return (b[i] - a[i]);
		}
	}

	return 0;
}

// Allocate a tree node containing array a of size n
struct tree *tree_malloc(const int *a, const int n, const int val)
{
	// Allocate node and space for array with one malloc
	struct tree *root = malloc(sizeof(struct tree) + n*sizeof(int));

	root->left = NULL;
	root->right = NULL;
	root->subtree = NULL;

	root->a = (int *)(root+1);
	memcpy(root->a, a, n*sizeof(int));

	root->num = 1;
	root->val = val;

	return root;
}

// Free a tree node; recursively calls itself to free its children
void tree_free(struct tree *root)
{
	if (root) {
		tree_free(root->left);
		tree_free(root->right);
		tree_free(root->subtree);

		free(root);
	}
}

// Add a sequence a with length n to the tree. The algorithm is to
// straightforwardly traverse the tree until an empty leaf is found where the
// array should be placed. The index i is used to count the subtree traversals.
//
// WARNING: Only add something to the tree if you're very sure that it is not
// contained in the tree yet. Otherwise you'll get some nasty errors.
void tree_add(struct tree *root, const int *a, const int n, const int val)
{
	struct tree *ptr = root;

	int i = 0;

	while (ptr) {
		ptr->num++;

		if (ptr->a[0] == a[i]) {
			i++;

			if (ptr->subtree == NULL) {
				ptr->subtree = tree_malloc(a+i, n-i, val);
				break;
			} else {
				ptr = ptr->subtree;
			}
		} else if (a[i] < ptr->a[0]) {
			if (ptr->left == NULL) {
				ptr->left = tree_malloc(a+i, n-i, val);
				break;
			} else {
				ptr = ptr->left;
			}
		} else {
			if (ptr->right == NULL) {
				ptr->right = tree_malloc(a+i, n-i, val);
				break;
			} else {
				ptr = ptr->right;
			}
		}
	}
}

// Verify whether an array a of size n is in the tree. The idea is exactly
// similar to tree_add, though we need to compare all elements this time to
// ensure that the current node is the correct one.
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

// Extract a list from the tree. The extracted list will be stored in *a, and
// marked as "processed" by reducing the field num. We need to find the deepest
// unprocessed child, since otherwise the code might construct unreachable
// children.
int tree_extract(struct tree *root, int *a, const int n)
{
	struct tree *ptr = root;

	int i = 0;

	while (i < n) {
		ptr->num--;

		if (ptr->subtree && ptr->subtree->num > 0) {
			a[i++] = ptr->a[0];
			ptr = ptr->subtree;
		} else if (ptr->left && ptr->left->num > 0) {
			ptr = ptr->left;
		} else if (ptr->right && ptr->right->num > 0) {
			ptr = ptr->right;
		} else {
			// No unprocessed children anymore. Process this one.
			memcpy(a+i, ptr->a, (n-i)*sizeof(int));
			return ptr->val;
		}
	}

	return -1;
}
