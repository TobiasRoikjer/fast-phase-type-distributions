#ifndef AMAZEPHASE_COAL_H
#define AMAZEPHASE_COAL_H

#include "phdist.h"
#include "dist.h"
#include <stdlib.h>

int coal_gen_phdist(phdist_t **phdist, size_t state_size);
int coal_seg_sites(d_dist_t **dist, phdist_t *phdist);


typedef struct d_phgen_args {
    mat_t *reward;
    double theta;
} d_phgen_args_t;

int d_ph_gen_fun(double **out, size_t from, size_t to, void *args);

#endif //AMAZEPHASE_COAL_H
