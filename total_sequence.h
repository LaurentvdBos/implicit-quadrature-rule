#ifndef TOTAL_SEQUENCE_H
#define TOTAL_SEQUENCE_H

struct node
{
	struct node *next, *prev;
	int val;
};

struct total_sequence
{
	struct node *root, *leaf, *curr;

	int k;
	int d;
};

struct total_sequence *total_sequence_malloc(const int d);
void total_sequence_free(struct total_sequence *ptr);
void total_sequence_next(struct total_sequence *ptr);

int total_sequence_sum(const struct total_sequence *ptr);
void total_sequence_set_sum(struct total_sequence *ptr, const int k);

#endif
