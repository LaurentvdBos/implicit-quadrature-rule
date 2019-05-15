#ifndef STACK_H
#define STACK_H

#include <stdlib.h>
#include <string.h>

struct stack
{
	int val;
	int *y;
	struct stack *next;
};

static struct stack *stack_push(struct stack *root, int val, int *y, int n)
{
	struct stack *node = malloc(sizeof(struct stack) + n*sizeof(int));

	node->val = val;
	node->y = (int *)(node + 1);
	memcpy(node->y, y, n*sizeof(int));
	node->next = root;

	return node;
}

static struct stack *stack_pop(struct stack *root)
{
	struct stack *node = root->next;

	free(root);

	return node;
}

#endif
