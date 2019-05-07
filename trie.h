#ifndef TRIE_H
#define TRIE_H

#include <stdbool.h>

struct trie
{
	// Number of children
	int len;

	// Children and their indices
	struct trie **next;
	int *index;

	// Whether this one ends a sequence
	bool end;
};

struct trie *trie_alloc();
void trie_free(struct trie *root);
void trie_add(struct trie *root, int *a, int n);
bool trie_contains(struct trie *root, int *a, int n);

#endif
