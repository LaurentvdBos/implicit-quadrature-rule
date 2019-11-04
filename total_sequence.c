#include "total_sequence.h"
#include <stdlib.h>

// Allocate a node of a linked list
static inline struct node *node_malloc()
{
	struct node *tmp = malloc(sizeof(struct node));
	tmp->next = NULL;
	tmp->prev = NULL;
	tmp->val = 0;

	return tmp;
}

// Puts a after b; it requires some bookkeeping in case a or b is the root or the leaf
static void node_splice_next(struct total_sequence *ptr, struct node *b, struct node *a)
{
	// Remove a
	if (a == ptr->root) {
		ptr->root = ptr->root->next;
	}
	if (a == ptr->leaf) {
		ptr->leaf = ptr->leaf->prev;
	}
	if (a->next) {
		a->next->prev = a->prev;
	}
	if (a->prev) {
		a->prev->next = a->next;
	}
	
	// Insert a
	if (b->next) {
		b->next->prev = a;
	}
	a->next = b->next;
	a->prev = b;
	b->next = a;

	if (b == ptr->leaf) {
		ptr->leaf = a;
	}
}

// Puts a before b; it requires some bookkeeping in case a or b is the root or leaf
static void node_splice_prev(struct total_sequence *ptr, struct node *b, struct node *a)
{
	// Remove a
	if (a == ptr->root) {
		ptr->root = ptr->root->next;
	}
	if (a == ptr->leaf) {
		ptr->leaf = ptr->leaf->prev;
	}
	if (a->next) {
		a->next->prev = a->prev;
	}
	if (a->prev) {
		a->prev->next = a->next;
	}
	
	// Insert a
	if (b->prev) {
		b->prev->next = a;
	}
	a->prev = b->prev;
	a->next = b;
	b->prev = a;

	if (b == ptr->root) {
		ptr->root = a;
	}
}

// Update the total sum of the sequence, which is basically unrolling the value
// of the sum into the linked list and taking the length constraint of the
// linked list into account.
//
// The linked list root is updated such that it contains the first composition
static void next_sum(struct total_sequence *ptr)
{
	int k = ptr->k;
	int d = ptr->d;

	if (k >= d) {
		ptr->leaf->val = k - d + 1;
	} else if (k > 0) {
		ptr->leaf->val = 1;
	} else {
		ptr->leaf->val = 0;
	}

	struct node *it = ptr->leaf->prev;
	int i = k - ptr->leaf->val;
	while (it) {
		if (i > 0) {
			it->val = 1;
			i--;
		} else {
			it->val = 0;
		}
		it = it->prev;
	}

	ptr->curr = ptr->root;
	ptr->k++;
}

// Determine next composition and store in linked list in ptr
// The algorithm is in-place and constant in time
//
// See also Kelleher [2004], slightly adapted for length constraint
// Before this function: ensure that root is in sorted, ascending order
// After this function: root is in sorted, ascending order containing the next
// composition
static void next_composition(struct total_sequence *ptr)
{
	struct node *it = ptr->leaf;

	if (it->val == ptr->k - 1) {
		next_sum(ptr);
	} else {
		int y = it->val - 1;
		it = it->prev;
		int x = it->val + 1;

		while (x <= y) {
			it->val = x;
			it = it->next;

			if (!it) {
				if (ptr->root->val == 0) {
					node_splice_next(ptr, ptr->leaf, ptr->root);
					it = ptr->leaf;
				} else {
					it = ptr->leaf;
					break;
				}
			}
			y -= x;
		}
		it->val = x+y;

		if (it->next) {
			it->next->val = 0;
			node_splice_prev(ptr, ptr->root, it->next);
		}
	}
	ptr->curr = ptr->root;
}

// Determine next permutation and store in linked list in ptr
// The algorithm is in-place and constant in time
//
// The main logic only works for sequences with more than 2 elements, so the
// special case of sequences of length 1 is caught immediately at entry of the
// function.
//
// See also Williams [2009] (this implementation is in *reverse* order)
// Before this function: ensure curr points to the largest index such that
// everything after curr is ascending
// After this function: root is the next permutation and curr has kept its
// property
static void next_permutation(struct total_sequence *ptr)
{
	if (ptr->d == 1) {
		return;
	}

	struct node *curr = ptr->curr;

	if (curr->prev && curr->prev->prev) {
		if (curr->prev->prev->val > curr->val) {
			node_splice_next(ptr, ptr->leaf, curr->prev);
		} else {
			node_splice_next(ptr, ptr->leaf, curr->prev->prev);
		}
	} else {
		node_splice_next(ptr, ptr->leaf, ptr->root);
		curr = ptr->root;
	}

	if (ptr->leaf->val < ptr->leaf->prev->val) {
		curr = ptr->leaf;
	}
	
	ptr->curr = curr;
}

// Allocate a total sequence containing d elements. It is initialized to
// contain a sequence that sums to 0.
struct total_sequence *total_sequence_malloc(const int d)
{
	struct total_sequence *ptr = malloc(sizeof(struct total_sequence));

	ptr->d = d;
	ptr->k = 0;
	ptr->root = node_malloc();
	ptr->leaf = ptr->root;

	for (int i = 1; i < d; i++) {
		struct node *tmp = node_malloc();
		tmp->next = ptr->root;
		ptr->root->prev = tmp;
		ptr->root = tmp;
	}

	next_sum(ptr);
	next_permutation(ptr);

	return ptr;
}

// Free a total sequence and the linked list it contains
void total_sequence_free(struct total_sequence *ptr)
{
	struct node *tmp = ptr->root;
	while (tmp) {
		ptr->root = ptr->root->next;
		free(tmp);
		tmp = ptr->root;
	}

	free(ptr);
}

// Move the total sequence to the next sequence. It consist of going to the
// next permutation. We are at the last permutation if it is sorted, which is
// the case if curr points to the root (see next_permutation for the invariant
// here).
void total_sequence_next(struct total_sequence *ptr)
{
	if (ptr->curr == ptr->root) {
		next_composition(ptr);
	}
	
	next_permutation(ptr);
}

// Determine the total sum of the total sequence
int total_sequence_sum(const struct total_sequence *ptr)
{
	return ptr->k - 1;
}

// Set the sum of the total sequence to k
void total_sequence_set_sum(struct total_sequence *ptr, const int k)
{
	ptr->k = k;
	next_sum(ptr);
	next_permutation(ptr);
}
