#ifndef AMAZEPHASE_MAT_H
#define AMAZEPHASE_MAT_H

#include <stdlib.h>
#include "bbst.h"

typedef struct mat {
    avl_flat_tuple_t **rows;
    avl_flat_tuple_t **cols;

    size_t *max_row_keys;
    size_t *max_col_keys;

    size_t n_rows;
    size_t n_cols;
} mat_t;

int mat_malloc(mat_t **out, size_t rows, size_t cols);
int mat_flatten(mat_t *out, avl_mat_node_t **rows, avl_mat_node_t **cols);
int mat_clone(mat_t **out, mat_t *in);
int mat_scale_rows(mat_t **out, mat_t *in, mat_entry_t *scalars);
int mat_mult(mat_t **out, mat_t *left, mat_t *right);
int mat_inv(mat_t **out, mat_t *in);
int mat_inv_slow(mat_t **out, mat_t *in);
int mat_identity(mat_t **out, size_t size);
int mat_sub(mat_t **out, mat_t *left, mat_t *right);
int mat_mul_scalar(mat_t **out, mat_t *in, mat_entry_t scalar);
int mat_row_sums(mat_t **out, mat_t *in);
void mat_print_as_matrix(mat_t *mat);
void mat_print_as_matrix_col(mat_t *mat);

#endif //AMAZEPHASE_MAT_H