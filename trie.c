#include "trie.h"

#include <stdlib.h>

static inline int find(const int *index, const int what, const int n)
{
	for (int i = 0; i < n; i++) {
		if (index[i] == what) {
			return i;
		}
	}

	return -1;
}

struct trie *trie_alloc()
{
	struct trie *root = malloc(sizeof(struct trie));
	root->len = 0;
	root->next = NULL;
	root->index = NULL;
	root->end = false;

	return root;
}

void trie_free(struct trie *root)
{
	if (root) {
		for (int i = 0; i < root->len; i++) {
			trie_free(root->next[i]);
		}

		free(root->next);
		free(root->index);
		free(root);
	}
}

void trie_add(struct trie *root, int *a, int n)
{
	struct trie *ptr = root;

	for (int i = 0; i < n; i++) {
		int ind = find(ptr->index, a[i], ptr->len);

		if (ind == -1) {
			ptr->len++;

			ptr->next = realloc(ptr->next, ptr->len*sizeof(struct trie *));
			ptr->index = realloc(ptr->index, ptr->len*sizeof(int));
			
			ptr->next[ptr->len-1] = trie_alloc();
			ptr->index[ptr->len-1] = a[i];

			ind = ptr->len-1;
		}

		ptr = ptr->next[ind];
	}

	ptr->end = true;
}

bool trie_contains(struct trie *root, int *a, int n)
{
	struct trie *ptr = root;

	for (int i = 0; i < n; i++) {
		int ind = find(ptr->index, a[i], ptr->len);

		if (ind == -1) {
			return false;
		}

		ptr = ptr->next[ind];
	}

	return ptr->end;
}
