#include "trie.h"

#include <stdlib.h>

struct trie *trie_alloc()
{
	struct trie *root = malloc(sizeof(struct trie));
	root->len = 0;
	root->next = NULL;
	root->end = false;

	return root;
}

void trie_free(struct trie *root)
{
	if (root) {
		for (unsigned int i = 0; i < root->len; i++) {
			trie_free(root->next[i]);
		}

		free(root->next);
		free(root);
	}
}

void trie_add(struct trie *root, unsigned int *a, int n)
{
	struct trie *ptr = root;

	for (int i = 0; i < n; i++) {
		if (ptr->len <= a[i]) {
			ptr->next = realloc(ptr->next, (a[i]+1)*sizeof(struct trie *));
			for (int j = ptr->len; j <= a[i]; j++) {
				ptr->next[j] = NULL;
			}
			ptr->len = a[i]+1;
		}

		if (!ptr->next[a[i]]) {
			ptr->next[a[i]] = trie_alloc();
		}

		ptr = ptr->next[a[i]];
	}

	ptr->end = true;
}

bool trie_contains(struct trie *root, unsigned int *a, int n)
{
	struct trie *ptr = root;

	for (int i = 0; i < n; i++) {
		if (!ptr || ptr->len <= a[i]) {
			return false;
		}
		ptr = ptr->next[a[i]];
	}

	if (!ptr) {
		return false;
	} else {
		return ptr->end;
	}
}
