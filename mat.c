#include "mat.h"

#include <string.h>
#include <stdio.h>

static inline mat_entry_t mat_rc_sum(avl_flat_tuple_t *row, avl_node_t *col) {
    mat_entry_t res = 0;
    avl_node_t *node;

    while (row->entry != 0) {
        node = avl_find(col, row->key);

        if (node != NULL) {
            res += row->entry * node->entry;
        }

        row++;
    }

    return res;
}

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

int mat_identity(mat_t **out, size_t size) {
    *out = malloc(sizeof(mat_t));
    (*out)->n_rows = size;
    (*out)->n_cols = size;
    (*out)->rows = malloc(sizeof(avl_flat_tuple_t*)*(*out)->n_rows);
    (*out)->cols = malloc(sizeof(avl_flat_tuple_t*)*(*out)->n_cols);

    for (size_t r = 0; r < size; r++) {
        (*out)->rows[r] = malloc(sizeof(avl_flat_tuple_t)*2);
        (*out)->cols[r] = malloc(sizeof(avl_flat_tuple_t)*2);

        (*out)->rows[r][0].key = r;
        (*out)->rows[r][0].entry = 1;
        (*out)->rows[r][1].entry = 0;

        (*out)->cols[r][0].key = r;
        (*out)->cols[r][0].entry = 1;
        (*out)->cols[r][1].entry = 0;
    }

    return 0;
}

int mat_sub(mat_t **out, mat_t *left, mat_t *right) {
    *out = malloc(sizeof(mat_t));
    (*out)->n_rows = left->n_rows;
    (*out)->n_cols = left->n_cols;
    (*out)->rows = malloc(sizeof(avl_flat_tuple_t*)*(*out)->n_rows);
    (*out)->cols = malloc(sizeof(avl_flat_tuple_t*)*(*out)->n_cols);

    avl_node_t **rows = calloc((*out)->n_rows, sizeof(avl_node_t*));
    avl_node_t **cols = calloc((*out)->n_cols, sizeof(avl_node_t*));

    for (size_t r = 0; r < left->n_rows; r++) {
        avl_flat_tuple_t *row_p = left->rows[r];
        avl_flat_tuple_t *row2_p = right->rows[r];

        while (row_p->entry != 0 && row2_p->entry != 0) {
            mat_entry_t lv = row_p->entry;
            mat_entry_t rv = row2_p->entry;
            mat_entry_t res = 0;

            if (row_p->key == row2_p->key) {
                res = lv - rv;

                avl_insert_or_inc(&(rows[r]), row_p->key, res);
                avl_insert_or_inc(&(cols[row_p->key]), r, res);

                row_p++;
                row2_p++;
            } else if (row_p->key < row2_p->key) {
                res = lv;

                avl_insert_or_inc(&(rows[r]), row_p->key, res);
                avl_insert_or_inc(&(cols[row_p->key]), r, res);

                row_p++;
            } else {
                res = -rv;

                avl_insert_or_inc(&(rows[r]), row2_p->key, res);
                avl_insert_or_inc(&(cols[row2_p->key]), r, res);

                row2_p++;
            }
        }

        while (row_p->entry != 0) {
            avl_insert_or_inc(&(rows[r]), row_p->key, row_p->entry);
            avl_insert_or_inc(&(cols[row_p->key]), r, row_p->entry);

            row_p++;
        }

        while (row2_p->entry != 0) {
            avl_insert_or_inc(&(rows[r]), row2_p->key, -row2_p->entry);
            avl_insert_or_inc(&(cols[row2_p->key]), r, -row2_p->entry);

            row2_p++;
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

int mat_mul_scalar(mat_t **out, mat_t *in, mat_entry_t scalar) {
    mat_clone(out, in);

    for (size_t r = 0; r < in->n_rows; r++) {
        for (avl_flat_tuple_t *p = (*out)->rows[r]; p->entry != 0; p++) {
            p->entry *= scalar;
        }
    }

    for (size_t c = 0; c < in->n_cols; c++) {
        for (avl_flat_tuple_t *p = (*out)->cols[c]; p->entry != 0; p++) {
            p->entry *= scalar;
        }
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

int mat_row_sums(mat_t **out, mat_t *in) {
    //mat_clone(out, in);
    mat_entry_t *rows = malloc(sizeof(mat_entry_t)*in->n_rows);
    // To avoid zero-entries we start "all over"
    for (size_t r = 0; r < in->n_rows; r++) {
        mat_entry_t sum = 0;

        for (avl_flat_tuple_t *p = in->rows[r]; p->entry != 0; p++) {
            sum += p->entry;
        }

        rows[r] = sum;
    }

    (*out) = malloc(sizeof(mat_t));
    (*out)->rows = malloc(sizeof(avl_flat_tuple_t*)*in->n_rows);

    for (size_t r = 0; r < in->n_rows; r++) {
        (*out)->rows[r] = malloc(sizeof(avl_flat_tuple_t) * 2);
        (*out)->rows[r][0].key = 0;
        (*out)->rows[r][0].entry = rows[r];
        (*out)->rows[r][1].entry = 0;
    }

    (*out)->cols = malloc(sizeof(avl_flat_tuple_t*)*2);
    (*out)->cols[0] = malloc(sizeof(avl_flat_tuple_t)*(1+in->n_rows));
    (*out)->cols[1] = malloc(sizeof(avl_flat_tuple_t));

    size_t index = 0;

    for (size_t r = 0; r < in->n_rows; r++) {
        if (rows[r] == 0) {
            continue;
        }

        (*out)->cols[0][index].key = r;
        (*out)->cols[0][index].entry = rows[r];
        index++;
    }

    (*out)->cols[0][index].entry = -10;

    (*out)->cols[1][0].entry = -9;

    (*out)->n_rows = in->n_rows;
    (*out)->n_cols = 1;

    return 0;
}

int mat_inv(mat_t **out, mat_t *in) {
    *out = malloc(sizeof(mat_t));

    avl_node_t **rows = calloc(in->n_rows, sizeof(avl_node_t*));
    avl_node_t **cols = calloc(in->n_cols, sizeof(avl_node_t*));

    avl_flat_tuple_t **iter_r = malloc(sizeof(avl_flat_tuple_t*)*in->n_rows);

    for (size_t r = 0; r < in->n_rows; r++) {
        iter_r[r] = in->rows[r];
    }

    for (size_t r = 0; r < in->n_rows; r++) {
        // Assume all diagonals non-zero
        avl_insert_or_inc(&(rows[r]), r, 1/iter_r[r]->entry);
        avl_insert_or_inc(&(cols[r]), r, 1/iter_r[r]->entry);
        iter_r[r]++;
    }

    for (size_t d = 1; d < in->n_rows; d++) {
        for (size_t r = 0; r < in->n_rows - d; r++) {
            mat_key_t column_key = d + r;
            mat_entry_t res = -mat_rc_sum(in->rows[r], cols[column_key])/in->rows[r]->entry;

            if (res != 0) {
                avl_insert_or_inc(&(rows[r]), column_key, res);
                avl_insert_or_inc(&(cols[column_key]), r, res);
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