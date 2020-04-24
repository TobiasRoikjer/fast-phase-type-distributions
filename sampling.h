#include "utils.h"
#include "coal.h"

#ifndef AMAZEPHASE_SAMPLING_H
#define AMAZEPHASE_SAMPLING_H

typedef long double sampling_number_t;

typedef struct {
    sampling_number_t rate;
    sampling_number_t constant;
} pdf_constant_t;

int sampling_graph_iterative(sampling_number_t **out, coal_graph_node_t *graph, size_t reward_size);
int sampling_graph_pfd_constants(pdf_constant_t **out, size_t *out_size, coal_graph_node_t *graph, size_t reward_size);
int sampling_graph_pfd_constants_rec(pdf_constant_t **out, size_t *out_size, coal_graph_node_t *graph, size_t vector_len, size_t reward_index);
int sampling_graph_pfd_constants_rec_rw(pdf_constant_t **out, size_t *out_size, coal_graph_node_t *graph);
#endif //AMAZEPHASE_SAMPLING_H
