#ifndef AMAZEPHASE_DIST_H
#define AMAZEPHASE_DIST_H

#include <stdlib.h>

typedef int (*generator_fun_t)(ssize_t **out, size_t from, size_t to, void *args);

typedef struct d_dist {
    generator_fun_t generator_fun;
} d_dist_t;

#endif //AMAZEPHASE_DIST_H
