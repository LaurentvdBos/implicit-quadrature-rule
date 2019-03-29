#ifndef TRIE_H
#define TRIE_H

#include <stdbool.h>

struct trie
{
	unsigned int len;
	struct trie **next;

	bool end;
};

struct trie *trie_alloc();
void trie_free(struct trie *root);
void trie_add(struct trie *root, unsigned int *a, int n);
bool trie_contains(struct trie *root, unsigned int *a, int n);

#endif
