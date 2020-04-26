#include "../src/coal.h"
#include "../src/sampling.h"
#include <string.h>
#include <time.h>
#include <math.h>

size_t reward_i;
size_t reward_j;
size_t n1, n2;

double reward_tot_ss(coal_graph_node_t* node) {
    return node->data.state_vec[0] + node->data.state_vec[1];
}

double reward_tot(coal_graph_node_t* node) {
    double reward = 0;

    for (size_t i = 0; i <= n1; ++i) {
        for (size_t j = 0; j <= n2; ++j) {
            reward += (double)((((im_state_t *) node->data.state)->mat1)[i][j]);
            reward += (double)((((im_state_t *) node->data.state)->mat2)[i][j]);
        }
    }

    return reward;
}

double reward(coal_graph_node_t* node) {
    double reward = (double)((((im_state_t *) node->data.state)->mat1)[reward_i][reward_j]) +
                    (double)((((im_state_t *) node->data.state)->mat2)[reward_i][reward_j]);

    return reward;
}

clock_t total_conv = 0;

double reward_one(coal_graph_node_t *node) {
    return 1;
}


void print_cdf_graph(FILE *out, coal_graph_node_t *start) {
    coal_rewards_set(start,reward_one);
    start->data.reward = 0;

    pdf_constant_t *constants;
    size_t size;

    sampling_graph_pfd_constants_rec_rw(&constants, &size, start);
    long double const_sum = 0;
    for (size_t i = 0; i < size; ++i) {
        const_sum += constants[i].constant;
        fprintf(stderr, "Adding %Lf, sum is now %Lf\n", constants[i].constant, const_sum);
    }

    for (double t = 0; t <= 30; t += 0.5) {
        fprintf(out, "%f ", t);
    }

    fprintf(out, "\n");

    for (double t = 0; t <= 30; t += 0.5) {
        sampling_number_t cdf = 0;

        for (size_t i = 0; i < size; ++i) {
            cdf += constants[i].constant * (1 - expl(-constants[i].rate * t));
        }

        fprintf(out, "%Lf ", cdf);
    }

    fprintf(out, "\n");
}


void print_cdf_mat(FILE *out, coal_graph_node_t *start) {
    long double res;
    gsl_matrix *S;
    coal_get_as_mat(&S, start);

    for (double t = 0; t <= 30; t += 0.5) {
        fprintf(out, "%f ", t);
    }

    fprintf(out, "\n");

    for (double t = 0; t <= 30; t += 0.5) {
        coal_get_mat_cdf(&res, t, S, start);
        fprintf(out, "%Lf ", res);
    }

    free(S);

    fprintf(out, "\n");
}

int main(int argc, char **argv) {
    if (argc != 9 && argc != 10) {
        fprintf(stderr, "Expected 8 or 9 arguments got %i\n", argc - 1);
        exit(1);
    }

    char *filename = argv[1];
    n1 = (size_t) atoi(argv[2]);
    n2 = (size_t) atoi(argv[3]);
    coal_migration_param back_migrations = (coal_migration_param) atoi(argv[4]);
    coal_param_real_t pop_scale1 = (coal_param_real_t) atof(argv[5]);
    coal_param_real_t pop_scale2 = (coal_param_real_t) atof(argv[6]);
    coal_param_real_t mig_scale1 = (coal_param_real_t) atof(argv[7]);
    coal_param_real_t mig_scale2 = (coal_param_real_t) atof(argv[8]);
    bool ss_only = false;

    if (argc == 10) {
        ss_only = (bool) atof(argv[9]);
    }
    clock_t start;
    clock_t tot_gen = 0;
    clock_t tot_rw = 0;
    clock_t tot_write = 0;

    fprintf(stderr, "Args:\n\tfilename: %s\n\tn1: %zu n2: %zu\n\tM: %i,\n\tnu_p1: %Lf, nu_p2: %Lf,\n\tnu_m1: %Lf, nu_m2: %Lf,\n\tss_only: %i\n\n",
            filename, n1, n2, back_migrations, pop_scale1, pop_scale2, mig_scale1, mig_scale2, ss_only);

    coal_graph_node_t *no_iso_graph;
    avl_vec_node_t *no_iso_bst;

    coal_gen_im_graph_args_t args_no_iso = {
            .n1 = n1,
            .n2 = n2,
            .num_iso_coal_events = 0,
            .migration_type = back_migrations,
            .pop_scale1 = pop_scale1,
            .pop_scale2 = pop_scale2,
            .mig_scale1 = mig_scale1,
            .mig_scale2 = mig_scale2
    };

    start = clock();
    coal_gen_im_graph(&no_iso_graph, &no_iso_bst, args_no_iso);
    tot_gen += clock() - start;

    coal_graph_node_t *iso_graph;

    coal_gen_im_graph_args_t args_iso = {
            .n1 = n1,
            .n2 = n2,
            .num_iso_coal_events = n1+n2-1,
            .migration_type = back_migrations,
            .pop_scale1 = pop_scale1,
            .pop_scale2 = pop_scale2,
            .mig_scale1 = mig_scale1,
            .mig_scale2 = mig_scale2
    };

    start = clock();
    coal_gen_im_graph(&iso_graph, NULL, args_iso);
    tot_gen += clock() - start;

    for (ssize_t coals = n1+n2-1; coals >= 0; --coals) {
        coal_graph_node_t *orig_graph;

        if (coals == 0) {
            orig_graph = no_iso_graph;
        } else {
            orig_graph = iso_graph;
            start = clock();
            coal_graph_im_redirect_at_coals(orig_graph, (size_t) coals, no_iso_bst);
            tot_gen += clock() - start;
        }

        if (!ss_only) {
            // TODO: Redirect based on *reward transformed graphs* such that
            // we do not have to reward transform each time...
            for (reward_i = 0; reward_i <= n1; ++reward_i) {
                for (reward_j = 0; reward_j <= n2; ++reward_j) {
                    char file[2048];
                    snprintf(file, sizeof(file),
                             "%s-%zu-%zu.%zu.tsv", filename, coals, reward_i, reward_j);

                    FILE *f = fopen(file, "w");
                    if (f == NULL) {
                        DIE_PERROR(1, "Failed to open file %s\n", file);
                    }

                    coal_graph_node_t *graph;
                    coal_graph_node_t *start_node;
                    coal_graph_clone(&graph, orig_graph);

                    start = clock();
                    coal_rewards_set(graph, reward);
                    coal_reward_transform(graph, &start_node);
                    tot_rw += clock() - start;
                    start = clock();
                    print_cdf_graph(f, start_node);
                    tot_write += clock() - start;
                    fclose(f);
                }
            }
        }

        char file[2048];
        snprintf(file, sizeof(file),
                 "%s-%zu-TOTAL.tsv", filename, coals);

        FILE *f = fopen(file, "w");
        if (f == NULL) {
            DIE_PERROR(1, "Failed to open file %s\n", file);
        }

        coal_graph_node_t *graphss;
        coal_graph_node_t *start_nodess;

        coal_gen_im_graph_args_t args = {
                .n1 = n1,
                .n2 = n2,
                .num_iso_coal_events = (size_t)coals,
                .migration_type = back_migrations,
                .pop_scale1 = pop_scale1,
                .pop_scale2 = pop_scale2,
                .mig_scale1 = mig_scale1,
                .mig_scale2 = mig_scale2
        };

        start = clock();
        coal_gen_im_ss_graph(&graphss, args);
        tot_gen = clock() -start;
        start = clock();
        coal_rewards_set(graphss, reward_tot_ss);
        coal_reward_transform(graphss, &start_nodess);
        tot_rw += clock() - start;
        start = clock();
        print_cdf_graph(f, start_nodess);
        tot_write += clock() - start;
        fclose(f);
    }

    fprintf(stderr, "Total gen time %zu, total reward time %zu, total conv/write time %zu, total conv time %zu\n",
            tot_gen, tot_rw, tot_write,total_conv);


    return 0;
}