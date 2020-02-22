#include "utils.h"
#include "coal.h"

#ifndef AMAZEPHASE_SAMPLING_H
#define AMAZEPHASE_SAMPLING_H

typedef struct {
    double rate;
    double constant;
} pdf_constant_t;

int sampling_graph_iterative(double **out, coal_graph_node_t *graph, size_t reward_size);
int sampling_graph_pfd_constants(pdf_constant_t **out, size_t *out_size, coal_graph_node_t *graph, size_t reward_size);

#endif //AMAZEPHASE_SAMPLING_H
