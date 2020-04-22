#include "coal.h"
#include <time.h>

int main(int argc, char **argv) {
    if (argc != 9 && argc != 10) {
        fprintf(stderr, "Expected 8 or 9 arguments got %i\n", argc-1);
        exit(1);
    }

    double isolation_time = atof(argv[1]);
    size_t n1 = (size_t)atoi(argv[2]);
    size_t n2 = (size_t)atoi(argv[3]);
    coal_migration_param back_migrations = atoi(argv[4]);
    coal_param_real_t pop_scale1 = (coal_param_real_t)atof(argv[5]);
    coal_param_real_t pop_scale2 = (coal_param_real_t)atof(argv[6]);
    coal_param_real_t mig_scale1 = (coal_param_real_t)atof(argv[7]);
    coal_param_real_t mig_scale2 = (coal_param_real_t)atof(argv[8]);
    coal_param_real_t factor;

    if (argc == 10) {
        factor = (coal_param_real_t)atof(argv[9]);
    } else {
        factor = 1.0f;
    }

    fprintf(stderr, "Args:\n\tiso time: %f,\n\tn1:%zu, n2: %zu\n\tM: %i,\n\tnu_p1: %Lf, nu_p2: %Lf,\n\tnu_m1: %Lf, nu_m2: %Lf\n\tfactor: %Lf\n\n",
            isolation_time, n1, n2, back_migrations, pop_scale1, pop_scale2, mig_scale1, mig_scale2, factor);

    coal_gen_im_graph_args_t argsp = {
            .n1 = n1,
            .n2 = n2,
            .migration_type = MIG_ALL,
            .pop_scale1 = pop_scale1,
            .pop_scale2 = pop_scale2,
            .mig_scale1 = mig_scale1,
            .mig_scale2 = mig_scale2
    };

    long double* probs;
    coal_im_get_number_coals_probs(&probs, isolation_time, &argsp);

    gsl_matrix_long_double *exp_mat = gsl_matrix_long_double_calloc(n1+1,n2+1);
    clock_t start;
    clock_t tot_gen = 0;
    clock_t tot_exp = 0;

    for (size_t coals = 0; coals < n1 + n2; ++coals) {
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

        // TODO: Idea: copy this graph, don't generate it all the time?
        start = clock();
        coal_gen_im_graph(&graph, args);
        tot_gen += clock() - start;
        start = clock();
        coal_mph_im_expected(graph, n1, n2);
        tot_exp += clock() - start;
        gsl_matrix_long_double *mat = graph->data.pointer;
        gsl_matrix_long_double_scale(mat, (double) probs[coals]);
        gsl_matrix_long_double_add(exp_mat, mat);
    }

    fprintf(stderr, "Total gen time %zu, total exp time %zu\n",
            tot_gen, tot_exp);

    for (size_t i = 0; i <= n1; ++i) {
        for (size_t j = 0; j <= n2; ++j) {
            fprintf(stdout, "%Lf ", gsl_matrix_long_double_get(exp_mat, i, j)*factor);
        }

        fprintf(stdout, "\n");
    }

    return 0;
}