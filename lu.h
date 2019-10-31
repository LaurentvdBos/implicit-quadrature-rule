#ifndef LU_H
#define LU_H

#include "matrix.h"

void matrix_lu(struct matrix *mat);
void matrix_lu_solve(const struct matrix *mat, struct matrix *b);
void matrix_lu_null(const struct matrix *mat, struct matrix *c);

#endif
