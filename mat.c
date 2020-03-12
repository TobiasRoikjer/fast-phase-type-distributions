#include "mat.h"

#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

inline int mat_malloc(mat_t **out, size_t rows, size_t cols) {
    *out = malloc(sizeof(mat_t));
    (*out)->n_rows = rows;
    (*out)->n_cols = cols;
    (*out)->rows = malloc(sizeof(avl_flat_tuple_t*) * (*out)->n_rows);
    (*out)->cols = malloc(sizeof(avl_flat_tuple_t*) * (*out)->n_cols);

    // TODO: No need for calloc here, use malloc
    (*out)->max_row_keys = calloc((*out)->n_rows, sizeof(size_t));
    (*out)->max_col_keys = calloc((*out)->n_cols, sizeof(size_t));

    return 0;
}

static inline int avl_malloc(avl_mat_node_t ***rows, avl_mat_node_t ***cols, size_t nrows, size_t ncols) {
    *rows = calloc(nrows, sizeof(avl_mat_node_t*));
    *cols = calloc(ncols, sizeof(avl_mat_node_t*));

    return 0;
}

static inline int avl_single_malloc(avl_mat_node_t ***entries, size_t n) {
    *entries = calloc(n, sizeof(avl_mat_node_t*));

    return 0;
}

inline int mat_flatten(mat_t *out, avl_mat_node_t **rows, avl_mat_node_t **cols) {
    for (size_t i = 0; i < (out)->n_rows; i++) {
        avl_flatten(&(out->rows[i]), &(out->max_row_keys[i]), rows[i]);
    }

    for (size_t i = 0; i < (out)->n_cols; i++) {
        avl_flatten(&(out->cols[i]), &(out->max_col_keys[i]), cols[i]);
    }

    free(rows);
    free(cols);

    return 0;
}

static inline int mat_single_flatten(mat_t *out, avl_mat_node_t **entries) {
    for (size_t i = 0; i < (out)->n_rows; i++) {
        avl_flatten(&(out->rows[i]), &(out->max_row_keys[i]), entries[i]);
    }

    free(entries);

    return 0;
}

static inline mat_entry_t mat_rc_sum_slow(avl_flat_tuple_t *row, avl_mat_node_t *col) {
    mat_entry_t res = 0;
    avl_mat_node_t *node;

    while (row->entry != 0) {
        node = avl_find(col, row->key);

        if (node != NULL) {
            res += row->entry * node->entry;
        }

        row++;
    }

    return res;
}


typedef struct avl_flat_meta {
    avl_flat_tuple_t *arr;
    size_t entries;
    size_t size;
} avl_flat_meta_t;



static inline mat_entry_t mat_rc_sum(avl_flat_tuple_t *row, avl_flat_meta_t *col) {
    if (col->entries == 0) {
        return 0;
    }

    mat_entry_t res = 0;
    size_t c = col->entries - 1;

    while (row->entry != 0) {
        avl_flat_tuple_t entry = col->arr[c];

        while (c > 0 && col->arr[c].key < row->key) {
            c--;
        }

        if (c == 0 && col->arr[c].key < row->key) {
            break;
        }

        if (col->arr[c].key == row->key) {
            res += row->entry * col->arr[c].entry;
        }

        row++;
    }

    return res;
}



//#define MAT_MUL_EPSILON 0.0000001
#define MAT_MUL_EPSILON 0

int mat_mult(mat_t **out, mat_t *left, mat_t *right) {
    mat_malloc(out, left->n_rows, right->n_cols);
    avl_mat_node_t **rows;
    avl_mat_node_t **cols;
    avl_malloc(&rows, &cols, (*out)->n_rows, (*out)->n_cols);

    for (size_t r = 0; r < left->n_rows; r++) {
        //for (size_t c = r; c <= left->max_row_keys[r]; c++) {
        for (size_t c = r; c < right->n_cols; c++) {
            // Shortcut. If non-zero row entries match only
            //   zero col entries, break.
            if (left->max_row_keys[r] < right->cols[c]->key) {
                break;
            }

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

            if (fabs(res) > MAT_MUL_EPSILON) {
                avl_insert_or_inc(&(rows[r]), c, res);
                avl_insert_or_inc(&(cols[c]), r, res);
            }
        }
    }

    mat_flatten(*out, rows, cols);

    return 0;
}

int mat_identity(mat_t **out, size_t size) {
    mat_malloc(out, size, size);

    for (size_t r = 0; r < size; r++) {
        (*out)->rows[r] = malloc(sizeof(avl_flat_tuple_t)*2);
        (*out)->cols[r] = malloc(sizeof(avl_flat_tuple_t)*2);

        (*out)->max_row_keys[r] = r;

        (*out)->rows[r][0].key = r;
        (*out)->rows[r][0].entry = 1;
        (*out)->rows[r][1].entry = 0;

        (*out)->max_col_keys[r] = r;

        (*out)->cols[r][0].key = r;
        (*out)->cols[r][0].entry = 1;
        (*out)->cols[r][1].entry = 0;
    }

    return 0;
}

int mat_sub(mat_t **out, mat_t *left, mat_t *right) {
    mat_malloc(out, left->n_rows, left->n_cols);
    avl_mat_node_t **rows;
    avl_mat_node_t **cols;
    avl_malloc(&rows, &cols, (*out)->n_rows, (*out)->n_cols);

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

    mat_flatten(*out, rows, cols);

    return 0;
}


int mat_clone(mat_t **out, mat_t *in) {
    mat_malloc(out, in->n_rows, in->n_cols);

    for (size_t r = 0; r < in->n_rows; r++) {
        (*out)->max_row_keys[r] = in->max_row_keys[r];

        size_t entries = 1;

        for (avl_flat_tuple_t *p = in->rows[r]; p->entry != 0; p++, entries++);

        ((*out)->rows)[r] = malloc(sizeof(avl_flat_tuple_t) * entries);
        memcpy(((*out)->rows)[r], in->rows[r], entries * sizeof(avl_flat_tuple_t));
    }

    for (size_t c = 0; c < in->n_cols; c++) {
        (*out)->max_col_keys[c] = in->max_col_keys[c];

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
    avl_mat_node_t **rows;
    avl_mat_node_t **cols;
    avl_malloc(&rows, &cols, in->n_rows, in->n_cols);

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

    mat_malloc(out, in->n_rows, in->n_cols);
    mat_flatten(*out, rows, cols);

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

    mat_malloc(out, in->n_rows, 1);

    for (size_t r = 0; r < in->n_rows; r++) {
        (*out)->rows[r] = malloc(sizeof(avl_flat_tuple_t) * 2);
        (*out)->max_row_keys[r] = 0;

        (*out)->rows[r][0].key = 0;
        (*out)->rows[r][0].entry = rows[r];
        (*out)->rows[r][1].entry = 0;
    }

    (*out)->cols = malloc(sizeof(avl_flat_tuple_t*)*2);

    (*out)->max_col_keys = malloc(sizeof(size_t)*2);

    (*out)->max_col_keys[0] = in->n_rows - 1;
    (*out)->max_col_keys[1] = 0;

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

    (*out)->cols[0][index].entry = 0;

    (*out)->cols[1][0].entry = 0;

    (*out)->n_rows = in->n_rows;
    (*out)->n_cols = 1;

    return 0;
}

int mat_inv_slow(mat_t **out, mat_t *in) {
    avl_mat_node_t **rows;
    avl_mat_node_t **cols;
    avl_malloc(&rows, &cols, in->n_rows, in->n_cols);

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

            mat_entry_t res = -mat_rc_sum_slow(in->rows[r], cols[column_key])/in->rows[r]->entry;

            if (res != 0) {
                avl_insert_or_inc(&(rows[r]), column_key, res);
                avl_insert_or_inc(&(cols[column_key]), r, res);
            }
        }
    }

    mat_malloc(out, in->n_rows, in->n_cols);
    mat_flatten(*out, rows, cols);

    return 0;
}

static int arr_should_expand(avl_flat_meta_t *value) {
    return (value->size - 1 == value->entries);
}

static int arr_expand(avl_flat_meta_t *value) {
    value->arr = realloc(value->arr, sizeof(avl_flat_tuple_t) * value->size * 2);
    value->size *= 2;
}

static int arr_expand_by(avl_flat_meta_t *value, size_t n) {
    value->arr = realloc(value->arr, sizeof(avl_flat_tuple_t) * (value->size + n));
    value->size += n;
}

int mat_inv(mat_t **out, mat_t *in) {
    avl_mat_node_t **rows;
    avl_flat_meta_t *flat_arr = malloc(sizeof(avl_flat_meta_t) * in->n_cols);

    for (size_t c = 0; c < in->n_cols; c++) {
        flat_arr[c].arr = malloc(sizeof(avl_flat_tuple_t) * 2);
        flat_arr[c].entries = 0;
        flat_arr[c].size = 2;
    }

    avl_single_malloc(&rows, in->n_rows);

    for (size_t r = 0; r < in->n_rows; r++) {
        // Assume all diagonals non-zero
        avl_insert_or_inc(&(rows[r]), r, 1 / in->rows[r]->entry);

        arr_expand_by(&(flat_arr[r]), 1);

        flat_arr[r].arr[flat_arr[r].entries].key = r;
        flat_arr[r].arr[flat_arr[r].entries].entry = 1 / in->rows[r]->entry;
        flat_arr[r].entries += 1;
    }

    for (size_t d = 1; d < in->n_rows; d++) {
        for (size_t r = 0; r < in->n_rows - d; r++) {
            mat_key_t column_key = d + r;

            mat_entry_t res = -mat_rc_sum(in->rows[r], &(flat_arr[column_key]))/in->rows[r]->entry;

            if (res != 0) {
                if (arr_should_expand(&(flat_arr[column_key]))) {
                    arr_expand(&(flat_arr[column_key]));
                }

                avl_insert_or_inc(&(rows[r]), column_key, res);

                // We insert in order

                flat_arr[column_key].arr[flat_arr[column_key].entries].key = r;
                flat_arr[column_key].arr[flat_arr[column_key].entries].entry = res;
                flat_arr[column_key].entries++;
            }
        }
    }

    mat_malloc(out, in->n_rows, in->n_cols);
    mat_single_flatten(*out, rows);

    for (size_t i = 0; i < in->n_cols; i++) {
        (*out)->cols[i] = malloc(sizeof(avl_flat_tuple_t)*(flat_arr[i].entries+1));

        for (size_t c = 0; c < flat_arr[i].entries; c++) {
            (*out)->cols[i][flat_arr[i].entries - c - 1].key = flat_arr[i].arr[c].key;
            (*out)->cols[i][flat_arr[i].entries - c - 1].entry = flat_arr[i].arr[c].entry;
        }

        (*out)->cols[i][flat_arr[i].entries].entry = 0;
    }

    return 0;
}

void mat_print_as_matrix(mat_t *mat) {
    for (size_t r = 0; r < mat->n_rows; r++) {
        avl_flat_tuple_t *row = mat->rows[r];

        for (size_t c = 0; c < mat->n_cols; c++) {
            if (row->key == c && row->entry != 0) {
                printf("%f\t", row->entry);
                row++;
            } else {
                printf("0\t");
            }
        }

        printf("\n");
    }
}

void mat_print_as_matrix_col(mat_t *mat) {
    for (size_t r = 0; r < mat->n_rows; r++) {
        for (size_t c = 0; c < mat->n_cols; c++) {
            avl_flat_tuple_t *col = mat->cols[c];

            while (col->entry != 0) {
                if (col->key == r) {
                    break;
                }

                col++;
            }

            if (col->entry != 0 && col->key == r) {
                printf("%f\t", col->entry);
            } else {
                printf("0\t");
            }
        }

        printf("\n");
    }
}


int mat_repeat_col(mat_t **out, size_t *arr, size_t nmemb, size_t rows) {
    //TODO
    return 0;
}