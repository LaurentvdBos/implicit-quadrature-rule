#ifndef TOTAL_SEQUENCE_H
#define TOTAL_SEQUENCE_H

struct node
{
	struct node *next, *prev;
	unsigned int val;
};

struct total_sequence
{
	struct node *root, *leaf, *curr;

	unsigned int k;
	unsigned int d;
};

struct total_sequence *total_sequence_alloc(unsigned int d);
void total_sequence_free(struct total_sequence *ptr);
void total_sequence_next(struct total_sequence *ptr);

unsigned int total_sequence_sum(struct total_sequence *ptr);
void total_sequence_set_sum(struct total_sequence *ptr, unsigned int k);

#endif
