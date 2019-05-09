#ifndef TREE_H
#define TREE_H

#include <stdbool.h>

struct tree
{
	// Left and right subtree
	struct tree *left, *right;

	// Array with indices, length is known by user
	int *a;
};

struct tree *tree_malloc(const int *a, const int n);
void tree_free(struct tree *root);
void tree_add(struct tree *root, int *a, int n);
bool tree_contains(struct tree *root, int *a, int n);

#endif
