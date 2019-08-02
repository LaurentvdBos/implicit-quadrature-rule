#ifndef LU_H
#define LU_H

#include "matrix.h"

void matrix_lu(struct matrix *mat);
void matrix_lu_solve(struct matrix *mat, struct matrix *b);
void matrix_lu_null(struct matrix *mat, struct matrix *c);

#endif
