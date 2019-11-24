#include <stdio.h>

#include "phdist.h"

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
}