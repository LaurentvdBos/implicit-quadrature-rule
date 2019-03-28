#ifndef GSORT_H
#define GSORT_H

// Gnome-sort is essentially insertion sort, but if the sequence is already
// sorted (or close to sorted) it is much more efficient than many other
// sorting algorithms.
void gsort(int *start, int *end);

#endif
