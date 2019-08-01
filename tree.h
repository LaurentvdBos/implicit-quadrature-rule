#ifndef TREE_H
#define TREE_H

#include <stdbool.h>

// This is a struct containing an "extended" binary tree consisting of
// sequences. The sequences are compared from left to right. If the first
// element *matches*, traverse to subtree, which contains another binary tree
// consisting of the remaining parts of the sequence. It is basically a ternary
// tree that combines a trie and a tree.
// 
// Example: A ternary tree that contains the elements
//   (3, 5, 7), (2, 5, 7), (3, 6, 8), (3, 6, 9), (3, 7, 9), (4, 7, 8)
// might look as follows:
//
// (3, 5, 7)->left = (2, 5, 7)
// (3, 5, 7)->right = (4, 7, 8)
// (3, 5, 7)->subtree = (3, 6, 8), stored as (6, 8)
// (6, 8)->right = (7, 9)
// (6, 8)->subtree = (6, 9), stored as (9)
//
// The advantage of this data structure is that the number of alloc'ed structs
// is the same number as the elements in the tree. A straightforward trie has
// too many structs and therefore uses too much memory to be useful.

struct tree
{
	// Left and right subtree
	struct tree *left, *right;

	// Subtree of subsequences
	struct tree *subtree;

	// Array with indices, length is known by user
	int *a;

	// Value of this node
	int val;

	// Number of unprocessed array elements in this part of the tree
	long num;
};

struct tree *tree_malloc(const int *a, const int n, const int val);
void tree_free(struct tree *root);
void tree_add(struct tree *root, const int *a, const int n, const int val);
bool tree_contains(const struct tree *root, const int *a, const int n);
int tree_extract(struct tree *root, int *a, const int n);

#endif
