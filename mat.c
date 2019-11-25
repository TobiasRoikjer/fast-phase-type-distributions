#include "mat.h"

int mat_mult(mat_t **out, mat_t *left, mat_t *right) {
    *out = malloc(sizeof(mat_t));
    (*out)->n_rows = left->n_rows;
    (*out)->n_cols = right->n_cols;
    (*out)->rows = malloc(sizeof(avl_flat_tuple_t*)*(*out)->n_rows);
    (*out)->cols = malloc(sizeof(avl_flat_tuple_t*)*(*out)->n_cols);

    avl_node_t **rows = calloc((*out)->n_rows, sizeof(avl_node_t*));
    avl_node_t **cols = calloc((*out)->n_cols, sizeof(avl_node_t*));

    for (size_t r = 0; r < left->n_rows; r++) {
        for (size_t c = 0; c < right->n_cols; c++) {
            mat_entry_t res = 0;
            avl_flat_tuple_t *row_p = left->rows[r];
            avl_flat_tuple_t *col_p = right->cols[c];

            while (row_p->entry != 0 && col_p->entry != 0) {
                if (row_p->key == col_p->key) {
                    res += row_p->entry * col_p->entry;
                    row_p++;
                    col_p++;
                } else if (row_p->key < col_p->key) {
                    row_p++;
                } else {
                    col_p++;
                }
            }

            if (res != 0) {
                avl_insert_or_inc(&(rows[r]), c, res);
                avl_insert_or_inc(&(cols[c]), r, res);
            }
        }
    }

    size_t n;

    for (size_t i = 0; i < (*out)->n_rows; i++) {
        avl_flatten(&((*out)->rows[i]), &n, rows[i]);
    }

    for (size_t i = 0; i < (*out)->n_cols; i++) {
        avl_flatten(&((*out)->cols[i]), &n, cols[i]);
    }

    free(rows);
    free(cols);

    return 0;
}