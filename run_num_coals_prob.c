#include "coal.h"
#include "sampling.h"
#include <time.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr, "Expected 8 arguments got %i\n", argc-1);
        exit(1);
    }

    double isolation_time = atof(argv[1]);
    size_t n1 = (size_t)atoi(argv[2]);
    size_t n2 = (size_t)atoi(argv[3]);
    coal_migration_param back_migrations = (coal_migration_param)atoi(argv[4]);
    coal_param_real_t pop_scale1 = (coal_param_real_t)atof(argv[5]);
    coal_param_real_t pop_scale2 = (coal_param_real_t)atof(argv[6]);
    coal_param_real_t mig_scale1 = (coal_param_real_t)atof(argv[7]);
    coal_param_real_t mig_scale2 = (coal_param_real_t)atof(argv[8]);

    fprintf(stderr, "Args:\n\tisolation time: %f,\n\tn1:%zu, n2: %zu\n\tM: %i,\n\tnu_p1: %Lf, nu_p2: %Lf,\n\tnu_m1: %Lf, nu_m2: %Lf\n\n",
            isolation_time, n1, n2, back_migrations, pop_scale1, pop_scale2, mig_scale1, mig_scale2);

    coal_gen_im_graph_args_t args = {
            .n1 = n1,
            .n2 = n2,
            .migration_type = back_migrations,
            .pop_scale1 = pop_scale1,
            .pop_scale2 = pop_scale2,
            .mig_scale1 = mig_scale1,
            .mig_scale2 = mig_scale2
    };

    long double* probs;
    coal_im_get_number_coals_probs(&probs, isolation_time, &args);

    fprintf(stdout, "coals\tprob\n");
    for (size_t coals = 0; coals < args.n1+args.n2; ++coals) {
        fprintf(stdout, "%zu\t%Lf\n", coals, probs[coals]);
    }



    return 0;
}