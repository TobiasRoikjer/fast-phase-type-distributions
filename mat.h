#ifndef AMAZEPHASE_MAT_H
#define AMAZEPHASE_MAT_H

#include <stdlib.h>
#include "bbst.h"

typedef struct mat {
    avl_flat_tuple_t **rows;
    avl_flat_tuple_t **cols;

    size_t n_rows;
    size_t n_cols;
} mat_t;

int mat_mult(mat_t *left, mat_t *right);

#endif //AMAZEPHASE_MAT_H