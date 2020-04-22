#include "coal.h"
#include <string.h>

size_t reward_i;
size_t reward_j;

double reward(coal_graph_node_t* node) {
    double reward = (double)((((im_state_t *) node->data.state)->mat1)[reward_i][reward_j]) +
                    (double)((((im_state_t *) node->data.state)->mat2)[reward_i][reward_j]);

    return reward;
}

void print_rw_mat(FILE *out, coal_graph_node_t *start) {
    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, start);

    for (size_t i = 2; i < size; ++i) {
        weighted_edge_t *values = vector_get(start->edges);
        weighted_edge_t *edge = NULL;

        for (size_t j = 0; j < vector_length(start->edges); ++j) {
            if (((coal_graph_node_t *) values[j].node)->data.vertex_index == i) {
                edge = &(values[j]);
                break;
            }
        }

        if (edge != NULL) {
            fprintf(out, "%Lf ", edge->weight);
        } else {
            fprintf(out, "0 ");
        }
    }

    fprintf(out, "\n\n");

    for (size_t i = 2; i < size; ++i) {
        for (size_t j = 2; j < size; ++j) {
            fprintf(out, "%Lf ", mat[i][j]);
        }

        fprintf(out, "\n");
    }
}

int main(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr, "Expected 8 arguments got %i\n", argc - 1);
        exit(1);
    }

    char *filename = argv[1];
    size_t n1 = (size_t) atoi(argv[2]);
    size_t n2 = (size_t) atoi(argv[3]);
    coal_migration_param back_migrations = (coal_migration_param) atoi(argv[4]);
    coal_param_real_t pop_scale1 = (coal_param_real_t) atof(argv[5]);
    coal_param_real_t pop_scale2 = (coal_param_real_t) atof(argv[6]);
    coal_param_real_t mig_scale1 = (coal_param_real_t) atof(argv[7]);
    coal_param_real_t mig_scale2 = (coal_param_real_t) atof(argv[8]);

    fprintf(stderr, "Args:\n\tfilename: %s\n\tM: %i,\n\tnu_p1: %Lf, nu_p2: %Lf,\n\tnu_m1: %Lf, nu_m2: %Lf\n\n",
            filename, back_migrations, pop_scale1, pop_scale2, mig_scale1, mig_scale2);

    for (size_t coals = 0; coals < n1 + n2; ++coals) {
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

                coal_graph_node_t *start;
                coal_gen_im_graph(&graph, NULL, args);
                coal_rewards_set(graph, reward);
                coal_reward_transform(graph, &start);
                print_rw_mat(f, start);
                fclose(f);
            }
        }
    }

    return 0;
}