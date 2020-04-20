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
    coal_migration_param back_migrations = (coal_migration_param)atoi(argv[5]);
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
            .migration_type = back_migrations,
            .pop_scale1 = pop_scale1,
            .pop_scale2 = pop_scale2,
            .mig_scale1 = mig_scale1,
            .mig_scale2 = mig_scale2
    };

    coal_gen_im_cutoff_graph(&graph, args);
    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, graph);

    for (size_t i = 1; i < size; ++i) {
        for (size_t j = 1; j < size; ++j) {
            fprintf(stdout, "%Lf ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }

    return 0;
}