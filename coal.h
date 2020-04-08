#ifndef AMAZEPHASE_COAL_H
#define AMAZEPHASE_COAL_H

#include "phdist.h"
#include "dist.h"
#include "utils.h"
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix_long_double.h>

typedef size_t vec_entry_t;

typedef struct {
    vec_entry_t **mat1;
    vec_entry_t **mat2;
    bool flag_mig1to2;
    bool flag_mig2to1;
    bool in_iso;
} im_state_t;

int im_state_init(im_state_t **out,
        size_t n1, size_t n2);

int im_state_as_vec(vec_entry_t **out,
                    im_state_t *state,
                    size_t n1, size_t n2);

size_t im_state_length(size_t n1, size_t n2);

typedef struct {
    vec_entry_t *state_vec;
    void *state;
    double full_path_value;
    double prob;
    double vertex_exp;
    double descendants_exp_sum;
    bool visited;
    ssize_t vertex_index;
    size_t reset_int;
    void *pointer;
} coal_graph_node_data_t;

typedef struct {
    vector_t *edges;
    vector_t *reverse_edges;
    coal_graph_node_data_t data;
} coal_graph_node_t;

typedef long double coal_param_real_t;
typedef long double coal_res_real_t;

int coal_gen_phdist(phdist_t **phdist, size_t state_size);
int coal_gen_erlang_phdist(phdist_t **phdist, size_t samples);
int coal_seg_sites(d_dist_t **dist, phdist_t *phdist);
int coal_gen_kingman_graph(coal_graph_node_t **graph, size_t n);

typedef struct {
    size_t n1;
    size_t n2;
    size_t num_iso_coal_events;
    bool allow_back_migrations;
    coal_param_real_t pop_scale1;
    coal_param_real_t pop_scale2;
    coal_param_real_t mig_scale1;
    coal_param_real_t mig_scale2;
} coal_gen_im_graph_args_t;

int coal_gen_im_graph(coal_graph_node_t **graph, coal_gen_im_graph_args_t args);

typedef struct {
    size_t n1;
    size_t n2;
    size_t left;
    bool allow_back_migrations;
    coal_param_real_t pop_scale1;
    coal_param_real_t pop_scale2;
    coal_param_real_t mig_scale1;
    coal_param_real_t mig_scale2;
} coal_gen_im_pure_cutoff_graph_args_t;

int coal_gen_im_pure_cutoff_graph(coal_graph_node_t **graph, coal_gen_im_pure_cutoff_graph_args_t args);

typedef struct {
    size_t n1;
    size_t n2;
    size_t left_n1;
    size_t left_n2;
    bool allow_back_migrations;
    coal_param_real_t pop_scale1;
    coal_param_real_t pop_scale2;
    coal_param_real_t mig_scale1;
    coal_param_real_t mig_scale2;
} coal_gen_im_cutoff_graph_args_t;

int coal_gen_im_cutoff_graph(coal_graph_node_t **graph, coal_gen_im_cutoff_graph_args_t args);
int coal_gen_im_prob_vertex_graph(coal_graph_node_t **graph, coal_graph_node_t **correct_vertex, coal_gen_im_cutoff_graph_args_t args);
int coal_gen_im_ss_graph(coal_graph_node_t **graph, coal_gen_im_graph_args_t args);

int coal_dist_im_exact_coals(coal_res_real_t *probs, coal_gen_im_graph_args_t args);

int coal_graph_as_mat(weight_t ***weights, size_t *out_size, coal_graph_node_t *graph);
int coal_graph_as_gsl_mat(gsl_matrix_long_double **weights, coal_graph_node_t *graph);
int coal_graph_as_phdist_rw(phdist_t **phdist, coal_graph_node_t *graph);
double coal_mph_expected(coal_graph_node_t *graph, size_t reward_index);
double coal_mph_cov(coal_graph_node_t *graph,
                    size_t reward_index_1,
                    size_t reward_index_2);
int coal_label_vertex_index(size_t *largest_index, coal_graph_node_t *graph);
void coal_graph_reset(coal_graph_node_t *graph);

void coal_print_graph_list(FILE *stream, coal_graph_node_t *graph,
                           bool indexed,
                           size_t vec_length, size_t vec_spacing);
void coal_print_graph_list_im(FILE *stream, coal_graph_node_t *graph,
                              bool indexed,
                              size_t vec_length, size_t vec_spacing,
                              size_t n1, size_t n2);
//print_graph_node(start, im_state_length(n1,n2),0);

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
