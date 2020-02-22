#include "dist.h"

#include <math.h>

double dist_rand01() {
    return random() * 1.0f / (1.0f * RAND_MAX);
}

void dist_sampling_set_random_seed(unsigned int seed) {
    srandom(seed);
}

double dist_sample_exp(double rate) {
    return -log(1 - dist_rand01())/rate;
}