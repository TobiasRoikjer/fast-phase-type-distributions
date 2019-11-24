#include <stdio.h>

#include "phdist.h"

size_t phdist_count_non_zeros(phdist_t *phdist) {
    size_t non_zeros = 0;

    for (size_t r = 0; r < phdist->si_mat->n_rows; r++) {
        avl_flat_tuple_t *row = phdist->si_mat->rows[r];

        while (row->entry != 0) {
            non_zeros++;
            row++;
        }
    }

    return non_zeros;
}

void phdist_print_as_matrix(phdist_t *phdist) {
    for (size_t r = 0; r < phdist->si_mat->n_rows; r++) {
        avl_flat_tuple_t *row = phdist->si_mat->rows[r];

        for (size_t c = 0; c < phdist->si_mat->n_cols; c++) {
            if (row->key == c) {
                printf("%zd\t", row->entry);
                row++;
            } else {
                printf("0\t");
            }
        }

        printf("\n");
    }

    printf("\n");

    for (size_t i = 0; i < phdist->n_rw_rows; i++) {
        for (size_t j = 0; j < phdist->n_rw_cols; j++) {
            printf("%zu\t", phdist->rw_arr[i][j]);
        }

        printf("\n");
    }
}

void phdist_print_as_matrix_col(phdist_t *phdist) {
    for (size_t r = 0; r < phdist->si_mat->n_cols; r++) {
        avl_flat_tuple_t *row = phdist->si_mat->cols[r];

        for (size_t c = 0; c < phdist->si_mat->n_rows; c++) {
            if (row->key == c) {
                printf("%zd\t", row->entry);
                row++;
            } else {
                printf("0\t");
            }
        }

        printf("\n");
    }

    printf("\n");

    for (size_t i = 0; i < phdist->n_rw_cols; i++) {
        for (size_t j = 0; j < phdist->n_rw_rows; j++) {
            printf("%zu\t", phdist->rw_arr[j][i]);
        }

        printf("\n");
    }
}