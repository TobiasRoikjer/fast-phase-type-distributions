#ifndef AMAZEPHASE_COAL_H
#define AMAZEPHASE_COAL_H

#include "dist.h"
#include "utils.h"
#include "bbst.h"
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix_long_double.h>
#include <gsl/gsl_matrix.h>

typedef size_t vec_entry_t;

typedef struct {
    vec_entry_t **mat1;
    vec_entry_t **mat2;
    size_t n1;
    size_t n2;
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
    long double full_path_value;
    long double prob;
    long double vertex_exp;
    long double descendants_exp_sum;
    bool visited;
    size_t visits;
    ssize_t vertex_index;
    size_t reset_int;
    void *pointer;
    double reward;
    size_t coals;
    size_t type;
} coal_graph_node_data_t;

typedef struct {
    vector_t *edges;
    // TODO: This should not be a weighted_edge_t
    vector_t *reverse_edges;
    coal_graph_node_data_t data;
} coal_graph_node_t;

typedef long double coal_param_real_t;
typedef long double coal_res_real_t;

int coal_gen_kingman_graph(coal_graph_node_t **graph, size_t n);

typedef enum {
    MIG_DIR,
    MIG_ALL,
    MIG_ONCE
} coal_migration_param;

typedef struct {
    size_t n1;
    size_t n2;
    size_t num_iso_coal_events;
    coal_migration_param migration_type;
    coal_param_real_t pop_scale1;
    coal_param_real_t pop_scale2;
    coal_param_real_t mig_scale1;
    coal_param_real_t mig_scale2;
} coal_gen_im_graph_args_t;

int coal_gen_im_graph(coal_graph_node_t **graph, avl_vec_node_t **bst, coal_gen_im_graph_args_t args);

typedef struct {
    size_t n1;
    size_t n2;
    size_t left;
    coal_migration_param migration_type;
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
    coal_migration_param migration_type;
    coal_param_real_t pop_scale1;
    coal_param_real_t pop_scale2;
    coal_param_real_t mig_scale1;
    coal_param_real_t mig_scale2;
} coal_gen_im_cutoff_graph_args_t;

int coal_gen_im_cutoff_graph(coal_graph_node_t **graph, coal_gen_im_cutoff_graph_args_t args);
int coal_gen_im_prob_vertex_graph(coal_graph_node_t **graph, coal_graph_node_t **correct_vertex, coal_gen_im_cutoff_graph_args_t args);
int coal_gen_im_ss_graph(coal_graph_node_t **graph, coal_gen_im_graph_args_t args);

int coal_im_get_number_coals_probs(long double **out,
                                   double isolation_time,
                                   const coal_gen_im_graph_args_t *args);
int coal_get_as_mat(gsl_matrix **S, const coal_graph_node_t *graph);
int coal_get_mat_cdf(long double *out,
                        double t,
                        gsl_matrix *S,
                        coal_graph_node_t *start);
int coal_graph_im_redirect_at_coals(coal_graph_node_t *graph, const size_t coals, const avl_vec_node_t *non_iso_bst);
int coal_graph_as_mat(weight_t ***weights, size_t *out_size, coal_graph_node_t *graph);
int coal_graph_as_gsl_mat(gsl_matrix_long_double **weights, coal_graph_node_t *graph, bool include_absorbing);
long double coal_mph_expected(coal_graph_node_t *graph, size_t reward_index);
gsl_matrix_long_double * coal_mph_im_expected(coal_graph_node_t *graph, size_t n1, size_t n2);
long double coal_mph_cov(coal_graph_node_t *graph,
                    size_t reward_index_1,
                    size_t reward_index_2);
int coal_label_vertex_index(size_t *largest_index, coal_graph_node_t *graph);
size_t coal_get_edges(coal_graph_node_t *graph);
void coal_graph_reset(coal_graph_node_t *graph);
void coal_graph_reset_visited(coal_graph_node_t *graph);
int coal_graph_clone(coal_graph_node_t **out, coal_graph_node_t *graph);

void coal_print_graph_list(FILE *stream, coal_graph_node_t *graph,
                           bool indexed,
                           size_t vec_length, size_t vec_spacing);
void coal_print_graph_list_im(FILE *stream, coal_graph_node_t *graph,
                              bool indexed,
                              size_t vec_length, size_t vec_spacing,
                              size_t n1, size_t n2);

int coal_rewards_set(coal_graph_node_t *graph, double(*reward_function)(coal_graph_node_t *node));
int coal_reward_transform(coal_graph_node_t *graph, coal_graph_node_t **start);

#endif //AMAZEPHASE_COAL_H
