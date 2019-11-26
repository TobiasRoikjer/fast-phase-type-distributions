#include "mat.h"

#include <string.h>

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

int mat_clone(mat_t **out, mat_t *in) {
    *out = malloc(sizeof(mat_t));
    (*out)->n_rows = in->n_rows;
    (*out)->n_cols = in->n_cols;
    (*out)->rows = malloc(sizeof(avl_flat_tuple_t*) *(*out)->n_rows);
    (*out)->cols = malloc(sizeof(avl_flat_tuple_t*) *(*out)->n_cols);

    for (size_t r = 0; r < in->n_rows; r++) {
        size_t entries = 1;

        for (avl_flat_tuple_t *p = in->rows[r]; p->entry != 0; p++, entries++);

        ((*out)->rows)[r] = malloc(sizeof(avl_flat_tuple_t) * entries);
        memcpy(((*out)->rows)[r], in->rows[r], entries * sizeof(avl_flat_tuple_t));
    }

    for (size_t c = 0; c < in->n_cols; c++) {
        size_t entries = 1;

        for (avl_flat_tuple_t *p = in->cols[c]; p->entry != 0; p++, entries++);

        ((*out)->cols)[c] = malloc(sizeof(avl_flat_tuple_t) * entries);
        memcpy(((*out)->cols)[c], in->cols[c], entries * sizeof(avl_flat_tuple_t));
    }

    return 0;
}

int mat_scale_rows(mat_t **out, mat_t *in, mat_entry_t *scalars) {
    //mat_clone(out, in);
    *out = malloc(sizeof(mat_t));
    avl_node_t **rows = calloc(in->n_rows, sizeof(avl_node_t*));
    avl_node_t **cols = calloc(in->n_cols, sizeof(avl_node_t*));

    // To avoid zero-entries we start "all over"
    for (size_t r = 0; r < in->n_rows; r++) {
        for (avl_flat_tuple_t *p = in->rows[r]; p->entry != 0; p++) {
            mat_entry_t value = scalars[r] * p->entry;

            if (value != 0) {
                avl_insert_or_inc(&(rows[r]), p->key, value);
            }
        }
    }

    for (size_t c = 0; c < in->n_cols; c++) {
        for (avl_flat_tuple_t *p = in->cols[c]; p->entry != 0; p++) {
            mat_entry_t value = scalars[p->key] * p->entry;

            if (value != 0) {
                avl_insert_or_inc(&(cols[c]), p->key, value);
            }
        }
    }

    size_t rn;

    (*out)->rows = malloc(sizeof(avl_flat_tuple_t*)*in->n_rows);
    (*out)->cols = malloc(sizeof(avl_flat_tuple_t*)*in->n_cols);
    (*out)->n_rows = in->n_rows;
    (*out)->n_cols = in->n_cols;

    for (size_t i = 0; i < in->n_rows; i++) {
        avl_flatten(&((*out)->rows[i]), &rn, rows[i]);
    }

    for (size_t i = 0; i < in->n_cols; i++) {
        avl_flatten(&((*out)->cols[i]), &rn, cols[i]);
    }

    return 0;
}

int mat_repeat_col(mat_t **out, size_t *arr, size_t nmemb, size_t rows) {
    *out = malloc(sizeof(mat_t));
    (*out)->n_rows = rows;
    (*out)->n_cols = nmemb;
    (*out)->rows;
    //TODO
    return 0;
}