#ifndef AMAZEPHASE_COAL_H
#define AMAZEPHASE_COAL_H

#include "phdist.h"
#include "dist.h"
#include "utils.h"
#include <stdlib.h>

typedef size_t vec_entry_t;

typedef struct {
    vec_entry_t *state;
    double alpha;
    ssize_t vertex_index;
} coal_graph_node_data_t;

typedef struct {
    vector_t *edges;
    coal_graph_node_data_t data;
} coal_graph_node_t;

int coal_gen_phdist(phdist_t **phdist, size_t state_size);
int coal_gen_erlang_phdist(phdist_t **phdist, size_t samples);
int coal_seg_sites(d_dist_t **dist, phdist_t *phdist);
int coal_gen_graph_reward(coal_graph_node_t **graph, size_t n, size_t reward_index);
int coal_graph_as_phdist(phdist_t **phdist, coal_graph_node_t *graph);
double coal_mph_expected(coal_graph_node_t *graph, size_t reward_index);

typedef struct d_phgen_args {
    mat_t *reward;
    double theta;
} d_phgen_args_t;

typedef struct coal_args_hobolth_t {
    size_t n;
} coal_args_hobolth_t;

typedef struct coal_args_compress_t {
    size_t n;
    size_t reward_index;
} coal_args_compress_t;

int d_ph_gen_fun(double **out, size_t from, size_t to, void *args);

#endif //AMAZEPHASE_COAL_H
