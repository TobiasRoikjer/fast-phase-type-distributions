#ifndef AMAZEPHASE_COAL_H
#define AMAZEPHASE_COAL_H

#include "phdist.h"
#include "dist.h"
#include <stdlib.h>

int coal_gen_phdist(phdist_t **phdist, size_t state_size);
int coal_seg_sites(d_dist_t **dist, phdist_t *phdist);

#endif //AMAZEPHASE_COAL_H
