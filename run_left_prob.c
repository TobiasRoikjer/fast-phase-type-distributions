#include "coal.h"

int main(int argc, char **argv) {
    if (argc != 10) {
        fprintf(stderr, "Expected 9 arguments got %i\n", argc-1);
        exit(1);
    }

    size_t n1 = (size_t)atoi(argv[1]);
    size_t n2 = (size_t)atoi(argv[2]);
    size_t left_n1 = (size_t)atoi(argv[3]);
    size_t left_n2 = (size_t)atoi(argv[4]);
    bool back_migrations = (bool)atoi(argv[5]);
    coal_param_real_t pop_scale1 = (coal_param_real_t)atof(argv[6]);
    coal_param_real_t pop_scale2 = (coal_param_real_t)atof(argv[7]);
    coal_param_real_t mig_scale1 = (coal_param_real_t)atof(argv[8]);
    coal_param_real_t mig_scale2 = (coal_param_real_t)atof(argv[9]);

    fprintf(stderr, "Args:\n\tn1:%zu, n2: %zu\n\tleft1:%zu, left2: %zu,\n\tM: %i,\n\tnu_p1: %Lf, nu_p2: %Lf,\n\tnu_m1: %Lf, nu_m2: %Lf\n\n",
            n1, n2, left_n1, left_n2, back_migrations, pop_scale1, pop_scale2, mig_scale1, mig_scale2);

    coal_graph_node_t *graph;
    
    coal_gen_im_cutoff_graph_args_t args = {
            .n1 = n1,
            .n2 = n2,
            .left_n1 = left_n1,
            .left_n2 = left_n2,
            .allow_back_migrations = back_migrations,
            .pop_scale1 = pop_scale1,
            .pop_scale2 = pop_scale2,
            .mig_scale1 = mig_scale1,
            .mig_scale2 = mig_scale2
    };

    coal_graph_node_t *correct;
    coal_gen_im_prob_vertex_graph(&graph, &correct, args);
    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, graph);

    for (size_t i = 1; i < size; ++i) {
        weight_t rowsum = 0;
        weight_t taken = 0;

        for (size_t j = 1; j < size; ++j) {
            if (i == j) {
                rowsum = -mat[i][j];
            } else {
                taken += mat[i][j];
            }
        }

        for (size_t j = 1; j < size; ++j) {
            if (i == j) {
                fprintf(stdout, "0 ");
            } else {
                fprintf(stdout, "%Lf ", mat[i][j]/rowsum);
            }
        }

        fprintf(stdout, "%Lf", (rowsum-taken)/rowsum);

        fprintf(stdout, "\n");
    }

    for (size_t j = 0; j < size-1; ++j) {
        fprintf(stdout, "0 ");
    }

    fprintf(stdout, "1\n");

    // Print correct vertex
    fprintf(stdout, "%zu\n", correct->data.vertex_index);

    return 0;
}