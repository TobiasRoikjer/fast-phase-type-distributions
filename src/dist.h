#ifndef AMAZEPHASE_DIST_H
#define AMAZEPHASE_DIST_H

#include <stdlib.h>
#include "utils.h"

typedef int (*generator_fun_t)(double **out, size_t from, size_t to, void *args);

typedef struct d_dist {
    generator_fun_t generator_fun;
    void *args;
} d_dist_t;

double dist_rand01();
void dist_sampling_set_random_seed(unsigned int seed);
weight_t dist_sample_exp(weight_t rate);

#endif //AMAZEPHASE_DIST_H
