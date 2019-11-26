#include <stdio.h>
#include <string.h>

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
            if (row->key == c && row->entry != 0) {
                printf("%f\t", row->entry);
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
                printf("%f\t", row->entry);
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

int phdist_clone(phdist_t **out, phdist_t *in) {
    *out = malloc(sizeof(phdist_t));
    (*out)->n_rw_rows = in->n_rw_rows;
    (*out)->n_rw_cols = in->n_rw_cols;
    mat_clone(&((*out)->si_mat), in->si_mat);

    (*out)->rw_arr = malloc(in->n_rw_rows * sizeof(size_t*));

    for (size_t r = 0; r < in->n_rw_rows; r++) {
        (*out)->rw_arr[r] = malloc(in->n_rw_cols * sizeof(size_t));
        memcpy((*out)->rw_arr[r], in->rw_arr[r], in->n_rw_cols * sizeof(size_t));
    }

    return 0;
}


// TODO: Move scalars into coal method
int phdist_reward_transform(phdist_t **out, phdist_t *phdist) {
    *out = malloc(sizeof(phdist_t));

    mat_t *scaled;
    mat_entry_t *scalars = calloc(phdist->n_rw_rows, sizeof(mat_entry_t));

    for (size_t r = 0; r < phdist->n_rw_rows - 1; r++) {
        // Do not include the n-ton
        for (size_t c = 0; c < phdist->n_rw_cols - 1; c++) {
            scalars[r] += phdist->rw_arr[r][c];
        }

        scalars[r] = 1/scalars[r];
    }

    mat_scale_rows(&scaled, phdist->si_mat, scalars);

    for (avl_flat_tuple_t *p = scaled->rows[scaled->n_rows-1]; p->entry != 0; p++) {
        p->entry = 0;
    }

    for (avl_flat_tuple_t *p = scaled->cols[scaled->n_cols-1]; p->entry != 0; p++) {
        p->entry = 0;
    }

    for (size_t r = 0; r < scaled->n_rows-1; r++) {
        avl_flat_tuple_t *p = scaled->rows[r];

        while(1) {
            if (p->entry == 0 || (p+1)->entry == 0) {
                break;
            }

            p++;
        }

        if (p->key == scaled->n_cols - 1) {
            p->entry = 0;
        }
    }

    for (size_t c = 0; c < scaled->n_cols-1; c++) {
        avl_flat_tuple_t *p = scaled->cols[c];

        while(1) {
            if (p->entry == 0 || (p+1)->entry == 0) {
                break;
            }

            p++;
        }

        if (p->key == scaled->n_rows - 1) {
            p->entry = 0;
        }
    }

    (*out)->si_mat = scaled;
    (*out)->rw_arr = NULL;
    (*out)->n_rw_rows = 0;
    (*out)->n_rw_cols = 0;
    (*out)->si_mat->n_rows = phdist->si_mat->n_rows - 1;
    (*out)->si_mat->n_cols = phdist->si_mat->n_cols - 1;

    return 0;
}

int d_ph_gen_fun(ssize_t **out, size_t from, size_t to, d_phgen_args_t *args) {
    /*
     *
df <- data.frame(i = seq(0, 60), prob = c(prob(v2$newinitprob, P, p, seq(0, 3)), prob2(v2$newinitprob, D, p, seq(4, 60))))

theta <- 4
P <- (diag(nrow(v2$newsubintensitymatrix)) - 2/theta * v2$newsubintensitymatrix) %>% solve()
p <- rep(1L, nrow(P)) - P %*% rep(1L, nrow(P))

prob2 <- function(pi, dia, p, ns) {
  #unlist(lapply(ns, function(n) {(pi %*% dia$P %*%(diag(diag(dia$D)^ n)) %*% dia$Pi %*% p)[1,1]}))

  unlist(lapply(ns, function(n) {(pi %*% dia$P %*%(diag(exp(log(diag(dia$D)) * n))) %*% dia$Pi %*% p)[1,1]}))
}

prob <- function(pi, Tmat, p, ns) {
  unlist(lapply(ns, function(n) {(pi %*% (Tmat %^% n) %*% p)[1,1]}))
}
     */


}