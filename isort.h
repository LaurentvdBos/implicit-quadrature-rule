#ifndef ISORT_H
#define ISORT_H

/**
 * The lists that we need to sort are quite short (length d) and often already
 * approximately sorted. Quicksort is relatively slow to work with these lists,
 * so we use an insertion sort.
 */
void isort(int *a, int n);

#endif
