#include "../src/coal.h"

int main(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr, "Expected 8 arguments got %i\n", argc-1);
        exit(1);
    }

    size_t coals = (size_t)atoi(argv[1]);
    size_t n1 = (size_t)atoi(argv[2]);
    size_t n2 = (size_t)atoi(argv[3]);
    coal_migration_param back_migrations = atoi(argv[4]);
    coal_param_real_t pop_scale1 = (coal_param_real_t)atof(argv[5]);
    coal_param_real_t pop_scale2 = (coal_param_real_t)atof(argv[6]);
    coal_param_real_t mig_scale1 = (coal_param_real_t)atof(argv[7]);
    coal_param_real_t mig_scale2 = (coal_param_real_t)atof(argv[8]);

    fprintf(stderr, "Args:\n\tcoals: %zu,\n\tn1:%zu, n2: %zu\n\tM: %i,\n\tnu_p1: %Lf, nu_p2: %Lf,\n\tnu_m1: %Lf, nu_m2: %Lf\n\n",
            coals, n1, n2, back_migrations, pop_scale1, pop_scale2, mig_scale1, mig_scale2);

    coal_graph_node_t *graph;

    coal_gen_im_graph_args_t args = {
            .n1 = n1,
            .n2 = n2,
            .num_iso_coal_events = coals,
            .migration_type = back_migrations,
            .pop_scale1 = pop_scale1,
            .pop_scale2 = pop_scale2,
            .mig_scale1 = mig_scale1,
            .mig_scale2 = mig_scale2
    };

    coal_gen_im_graph(&graph, NULL, args);
    coal_mph_im_expected(graph, n1, n2);
    gsl_matrix_long_double *mat = graph->data.pointer;

    for (size_t i = 0; i <= n1; ++i) {
        for (size_t j = 0; j <= n2; ++j) {
            fprintf(stdout, "%Lf ", gsl_matrix_long_double_get(mat, i, j));
        }

        fprintf(stdout, "\n");
    }

    return 0;
}