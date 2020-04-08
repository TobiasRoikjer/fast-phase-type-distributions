#ifndef AMAZEPHASE_PHDIST_H
#define AMAZEPHASE_PHDIST_H

#include "mat.h"

typedef struct phdist {
    mat_t *si_mat;
    size_t **rw_arr;
    size_t n_rw_rows;
    size_t n_rw_cols;
} phdist_t;

int phdist_clone(phdist_t **out, phdist_t *in);
int phdist_reward_transform(phdist_t **out, phdist_t *phdist);
size_t phdist_count_non_zeros(phdist_t *phdist);
void phdist_print_as_matrix(phdist_t *phdist);
void phdist_print_as_matrix_col(phdist_t *phdist);


#endif //AMAZEPHASE_PHDIST_H
