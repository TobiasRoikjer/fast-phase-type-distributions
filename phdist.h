#ifndef AMAZEPHASE_PHDIST_H
#define AMAZEPHASE_PHDIST_H

#include "mat.h"

typedef struct phdist {
    mat_t *si_mat;
    mat_t *rw_mat;
} phdist_t;

void phdist_print_as_matrix(phdist_t *phdist);

#endif //AMAZEPHASE_PHDIST_H
