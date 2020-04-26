#include <stdio.h>
#include <time.h>
#include "../coal.h"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>


void test_im_mat_utils() {
    im_state_t *state;
    im_state_init(&state, 3, 4);
    state->mat1[1][0]= 3;
    state->mat1[0][1]= 4;
    state->mat1[3][0]= 1;
    state->mat1[0][4]= 2;
    state->mat1[3][4]= 9;

    state->mat2[1][0]= 3;
    state->mat2[0][1]= 4;
    state->mat2[3][0]= 1;
    state->mat2[0][4]= 2;
    state->mat2[0][2]= 3;
    state->mat2[1][2]= 4;
    state->mat2[2][2]= 5;
    state->mat2[3][4]= 8;

    vec_entry_t *vec;
    im_state_as_vec(&vec, state, 3, 4);
    fflush(stdout);

    for (size_t i = 0; i < (3+1)*(4+1)*2; ++i) {
        fprintf(stderr, "%zu,", vec[i]);
        if ((i+1) % 4 == 0) {
            fprintf(stderr, "\n");
        }
    }

    fprintf(stderr, "\n");
    fflush(stderr);
}
void print_type(coal_graph_node_t *node, FILE *f) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;

    fprintf(f, "%zu, %zu\n", node->data.vertex_index, node->data.coals);

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        print_type((coal_graph_node_t*) values[i].node, f);
    }
}

void test_gen_im() {
    coal_graph_node_t *graph;
    coal_gen_im_graph_args_t args = {
            .n1 = 3,
            .n2 = 3,
            .migration_type = MIG_ALL,
            .num_iso_coal_events = 5,
            .pop_scale1 = 0.3f,
            .pop_scale2 = 0.7f,
            .mig_scale1 = 1.0f,
            .mig_scale2 = 1.0f
    };

    coal_gen_im_graph(&graph, NULL, args);
    //coal_print_graph_list_im(stdout, graph, true, (args.n1 + 1)*(args.n2 + 1)*2+3,
    //                      (args.n1+1), args.n1, args.n2);

    gsl_matrix_long_double *mat;
    coal_graph_as_gsl_mat(&mat, graph, false);

    FILE *f = fopen("../t.tsv","w");

    for (size_t i = 0; i < mat->size1; ++i) {
        for (size_t j = 0; j < mat->size2; ++j) {
            fprintf(f, "%Lf ", gsl_matrix_long_double_get(mat, i, j));
        }

        fprintf(f, "\n");
    }

    fclose(f);

    FILE *f2 = fopen("../t2.tsv","w");

    coal_graph_reset_visited(graph);
    print_type(graph, f2);

    fclose(f2);
}

void test_gen_im_ss() {
    coal_graph_node_t *graph;
    coal_gen_im_graph_args_t args = {
            .n1 = 2,
            .n2 = 2,
            .migration_type = false,
            .num_iso_coal_events = 0,
            .pop_scale1 = 10.0f,
            .pop_scale2 = 1000.0f,
            .mig_scale1 = 1.0f,
            .mig_scale2 = 1.0f
    };

    coal_gen_im_ss_graph(&graph, args);
    coal_print_graph_list(stdout, graph, true, 4, 4);

    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, graph);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            fprintf(stdout, "%Lf ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }
}

void test_gen_im_time() {
    for (size_t i = 0; i < 10; ++i) {
        coal_graph_node_t *graph;
        coal_gen_im_graph_args_t args = {
                .n1 = 10,
                .n2 = 10,
                .migration_type = false,
                .num_iso_coal_events = i,
                .pop_scale1 = 10.0f,
                .pop_scale2 = 1000.0f,
                .mig_scale1 = 1.0f,
                .mig_scale2 = 1.0f
        };

        coal_gen_im_graph(&graph, NULL, args);
    }
}

void test_gen_im_cutoff() {
    coal_graph_node_t *graph;

    coal_gen_im_cutoff_graph_args_t args = {
            .n1 = 5,
            .n2 = 3,
            .left_n1 = 3,
            .left_n2 = 3,
            .migration_type = true,
            .pop_scale1 = 1.0f,
            .pop_scale2 = 1.0f,
            .mig_scale1 = 1.0f,
            .mig_scale2 = 1.0f
    };

    coal_gen_im_cutoff_graph(&graph, args);
    coal_print_graph_list(stdout, graph, true, 4, 4);
    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, graph);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            fprintf(stdout, "%Lf ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }
}

void test_gen_im2() {
    coal_graph_node_t *graph;
    coal_gen_im_graph_args_t args = {
            .n1 = 2,
            .n2 = 2,
            .num_iso_coal_events = 3,
            .migration_type = true,
            .pop_scale1 = 1.0f,
            .pop_scale2 = 1.0f,
            .mig_scale1 = 1.0f,
            .mig_scale2 = 1.0f
    };


    coal_gen_im_graph(&graph, NULL, args);
    coal_print_graph_list_im(stdout, graph, true, (args.n1 + 1)*(args.n2 + 1)*2+3,
                             (args.n1+1), args.n1, args.n2);
    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, graph);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            fprintf(stdout, "%Lf ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }
}

void test_gen_im_prob() {
    coal_graph_node_t *graph;
    coal_gen_im_cutoff_graph_args_t args = {
            .n1 = 3,
            .n2 = 3,
            .left_n1 = 1,
            .left_n2 = 4,
            .migration_type = true,
            .pop_scale1 = 1.0f,
            .pop_scale2 = 1.0f,
            .mig_scale1 = 1.0f,
            .mig_scale2 = 0.0f
    };

    coal_graph_node_t *correct_vertex = NULL;

    coal_gen_im_prob_vertex_graph(&graph, &correct_vertex, args);
    coal_print_graph_list(stdout, graph, true, 4, 4);
    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, graph);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            fprintf(stdout, "%Lf ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }

    if (correct_vertex == NULL) {
        fprintf(stdout, "Correct NULL\n");
    } else {
        fprintf(stdout, "Correct %zu\n", correct_vertex->data.vertex_index);
    }
}
//TODO: The mat_prob function should only be called once
// for all iters of N1
void test_num_coals() {
    coal_gen_im_graph_args_t args = {
            .n1 = 3,
            .n2 = 3,
            .migration_type = true,
            .pop_scale1 = 0.3f,
            .pop_scale2 = 0.7f,
            .mig_scale1 = 0.01f,
            .mig_scale2 = 0.01f,
    };

    long double* probs;
    coal_im_get_number_coals_probs(&probs, 1.0f, &args);

    for (size_t coals = 0; coals < args.n1+args.n2; ++coals) {
        fprintf(stdout, "%zu coals: %Lf\n", coals, probs[coals]);
    }
}

double reward_by(coal_graph_node_t *node) {
    //return node->data.state_vec[1] + node->data.state_vec[2];
    if (node->data.state_vec[1] == 2) {
        return 10;
    } else if (node->data.state_vec[2] == 1) {
        return 1;
    }

    return 0;
}

void test_reward_transform() {
    coal_graph_node_t *graph;
    coal_gen_kingman_graph(&graph, 4);

    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, graph);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            //fprintf(stdout, "%Lf ", mat[i][j]);
        }

        //fprintf(stdout, "\n");
    }

    // Doubletons
    coal_rewards_set(graph, reward_by);
    coal_graph_node_t *start;
    coal_reward_transform(graph, &start);
    graph = start;

    coal_graph_as_mat(&mat, &size, graph);

    fprintf(stdout, "\n");
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            fprintf(stdout, "%Lf ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }
}

void test_redir() {
    size_t n1 = 3;
    size_t n2 = 3;

    for (size_t coals = 1; coals <= n1+n2-1; ++coals) {
        coal_graph_node_t *graph;

        coal_gen_im_graph_args_t args3 = {
                .n1 = n1,
                .n2 = n2,
                .num_iso_coal_events = coals,
                .migration_type = MIG_ALL,
                .pop_scale1 = 1,
                .pop_scale2 = 1,
                .mig_scale1 = 1,
                .mig_scale2 = 1
        };

        coal_gen_im_graph(&graph, NULL, args3);

        //coal_print_graph_list_im(stderr, graph, true, (n1 + 1) * (n2 + 1) * 2 + 3, n1+1, n1, n2);
        //fprintf(stderr, "\n========\n\n");

        coal_graph_node_t *no_iso_graph;

        coal_gen_im_graph_args_t args = {
                .n1 = n1,
                .n2 = n2,
                .num_iso_coal_events = 0,
                .migration_type = MIG_ALL,
                .pop_scale1 = 1,
                .pop_scale2 = 1,
                .mig_scale1 = 1,
                .mig_scale2 = 1
        };

        avl_vec_node_t *bst_no_iso;
        coal_gen_im_graph(&no_iso_graph, &bst_no_iso, args);

        coal_graph_node_t *iso_graph;

        coal_gen_im_graph_args_t args2 = {
                .n1 = n1,
                .n2 = n2,
                .num_iso_coal_events = n1 + n1 - 1,
                .migration_type = MIG_ALL,
                .pop_scale1 = 1,
                .pop_scale2 = 1,
                .mig_scale1 = 1,
                .mig_scale2 = 1
        };

        coal_gen_im_graph(&iso_graph, NULL, args2);

        coal_graph_im_redirect_at_coals(iso_graph, coals, bst_no_iso);

        coal_print_graph_list_im(stderr, iso_graph, false, (n1 + 1) * (n2 + 1) * 2 + 3, n1 + 1, n1, n2);


        weight_t **mat;
        weight_t **mat2;
        size_t size;
        size_t size2;
        coal_graph_as_mat(&mat, &size, graph);

        fprintf(stdout, "\n");
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                fprintf(stdout, "%Lf ", mat[i][j]);
            }

            fprintf(stdout, "\n");
        }

        fprintf(stdout, "\n");
        coal_graph_as_mat(&mat2, &size2, iso_graph);

        fprintf(stdout, "\n");
        for (size_t i = 0; i < size2; ++i) {
            for (size_t j = 0; j < size2; ++j) {
                fprintf(stdout, "%Lf ", mat2[i][j]);

                if (fabsl(mat2[i][j] - mat[i][j]) > 0.01) {
                    fflush(stdout);
                    DIE_ERROR(1, "Mat diff!\n");
                }
            }

            fprintf(stdout, "\n");
        }
    }
}

void test_clone_graph() {
    coal_graph_node_t *graph;
    coal_graph_node_t *cloned;

    coal_gen_im_graph_args_t args = {
            .n1 = 3,
            .n2 = 3,
            .num_iso_coal_events = 1,
            .migration_type =  MIG_ALL,
            .pop_scale1 = 1,
            .pop_scale2 = 1,
            .mig_scale1 = 1,
            .mig_scale2 = 1
    };

    coal_gen_im_graph(&graph, NULL, args);

    coal_graph_clone(&cloned, graph);

    weight_t **mat;
    weight_t **mat2;
    weight_t **mat3;
    size_t size;
    size_t size2;
    size_t size3;
    coal_graph_as_mat(&mat, &size, graph);
    coal_graph_as_mat(&mat2, &size2, cloned);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            if (fabsl(mat2[i][j] - mat[i][j]) > 0.01) {
                fflush(stdout);
                DIE_ERROR(1, "Mat diff!\n");
            }
        }
    }

    // Graph should not change modifying cloned
    coal_rewards_set(cloned, reward_by);
    coal_graph_node_t *start;
    coal_reward_transform(cloned, &start);

    coal_graph_as_mat(&mat3, &size3, graph);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            if (fabsl(mat3[i][j] - mat[i][j]) > 0.01) {
                fflush(stdout);
                DIE_ERROR(1, "Mat diff!\n");
            }
        }
    }
}

int main(int argc, char **argv) {
    //test_gen();
    //printf("\n..\n");
    //test_gen_erlang();
    //printf("\n..\n");
    /*//test_clone();
    printf("\n..\n");
    //test_zeroes();
    printf("\n..\n");
    //test_zeroes2();
    printf("\n..\n");
    test_scalar();
    printf("\n..\n");
    test_reward_sites();
    printf("\n..\n");*/
    //test_inverse();
    //printf("\n..\n");
    //test_time_inverse();
    //printf("\n..\n");
    /*test_mat_sub();
    printf("\n..\n");
    test_mat_scalar();
    printf("\n..\n");
    test_mat_rowsum();
    printf("\n..\n");
    test_seg();
    printf("\n..\n");*/
    //test_time_seg();
    //printf("\n..\n");
    //test_time_seg_erlang();
    //printf("\n..\n");
    //test_time_mul();
    //printf("\n..\n");
    //test_mat_mul();
    //printf("\n..\n");
    //test_gen_reward();
    //printf("\n..\n");
    //test_exp_reward();
    //printf("\n..\n");
    //test_cov_reward();
    //printf("\n..\n");
    //test_im_mat_utils();
    //printf("\n..\n");
    //test_gen_im();
    //printf("\n..\n");
    //test_gen_im_time();
    //printf("\n..\n");
    //test_gen_im_cutoff();
    //printf("\n..\n");
    //test_gen_im2();
    //printf("\n..\n");
    //test_gen_im_prob();
    //printf("\n..\n");
    //test_gen_im_ss();
    //printf("\n..\n");
    //test_num_coals();
    //printf("\n..\n");
    //test_reward_transform();
    //printf("\n..\n");
    //test_redir();
    //printf("\n..\n");
    test_clone_graph();
    printf("\n..\n");
    return 0;
}