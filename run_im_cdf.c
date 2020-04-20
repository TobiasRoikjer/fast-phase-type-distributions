#include "coal.h"

size_t reward_i;
size_t reward_j;

double reward(coal_graph_node_t* node) {
    double reward = (double)((((im_state_t *) node->data.state)->mat1)[reward_i][reward_j]) +
                    (double)((((im_state_t *) node->data.state)->mat2)[reward_i][reward_j]);

    return reward;
}

int main(int argc, char **argv) {
    if (argc != 11) {
        fprintf(stderr, "Expected 10 arguments got %i\n", argc-1);
        exit(1);
    }

    reward_i = (size_t)atoi(argv[1]);
    reward_j = (size_t)atoi(argv[2]);
    size_t coals = (size_t)atoi(argv[3]);
    size_t n1 = (size_t)atoi(argv[4]);
    size_t n2 = (size_t)atoi(argv[5]);
    coal_migration_param back_migrations = (coal_migration_param)atoi(argv[6]);
    coal_param_real_t pop_scale1 = (coal_param_real_t)atof(argv[7]);
    coal_param_real_t pop_scale2 = (coal_param_real_t)atof(argv[8]);
    coal_param_real_t mig_scale1 = (coal_param_real_t)atof(argv[9]);
    coal_param_real_t mig_scale2 = (coal_param_real_t)atof(argv[10]);

    fprintf(stderr, "Args:\n\treward index: %zu,%zu\n\tcoals: %zu,\n\tn1:%zu, n2: %zu\n\tM: %i,\n\tnu_p1: %Lf, nu_p2: %Lf,\n\tnu_m1: %Lf, nu_m2: %Lf\n\n",
            reward_i, reward_j, coals, n1, n2, back_migrations, pop_scale1, pop_scale2, mig_scale1, mig_scale2);

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
    coal_gen_im_graph(&graph, args);
    coal_rewards_set(graph, reward);
    coal_reward_transform(graph, &start);

    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, start);

    //weight_t multiplier = -1/mat[1][1];

    for (size_t i = 2; i < size; ++i) {
        weighted_edge_t *values = vector_get(start->edges);
        weighted_edge_t *edge = NULL;

        for (size_t j = 0; j < vector_length(start->edges); ++j) {
            if (((coal_graph_node_t*)values[j].node)->data.vertex_index == i) {
                edge = &(values[j]);
                break;
            }
        }

        if (edge != NULL) {
            fprintf(stdout, "%Lf ", edge->weight);
        } else {
            fprintf(stdout, "0 ");
        }
    }

    fprintf(stdout, "\n\n");

    for (size_t i = 2; i < size; ++i) {
        for (size_t j = 2; j < size; ++j) {
            fprintf(stdout, "%Lf ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }

    return 0;
}