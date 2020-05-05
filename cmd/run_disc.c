#include "../src/coal.h"

size_t rewards[1000];
size_t n;

double reward_sum(coal_graph_node_t *node) {
    double reward = 0;
    for (size_t i = 0; i < n; ++i) {
        if (rewards[i] != 0) {
            reward += node->data.state_vec[i];
        }
    }

    return reward;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Expected at least 3 arguments got %i\n", argc-1);
        exit(1);
    }

    n = (size_t)atoi(argv[1]);
    double theta = (double)atof(argv[2]);
    
    for (size_t i = 0; i < n; ++i) {
        rewards[i] = (size_t)atoi(argv[3+i]);
    }

    fprintf(stderr, "Args:\n\tn: %zu,\n\ttheta:%f\n\n",
            n, theta);

    gsl_matrix_long_double *mat;
    coal_graph_node_t *graph;

    coal_gen_kingman_graph(&graph, n);

    coal_rewards_set(graph, reward_sum);
    coal_construct_unshifted_discrete(graph, theta);

    size_t *weights = calloc(n, sizeof(size_t));

    for (size_t i = 0; i < n; ++i) {
        weights[i] = rewards[i];
    }

    coal_unshifted_discrete_apply_weighting(graph, weights, n);

    coal_graph_as_gsl_mat_discrete(&mat, graph, false);

    for (size_t i = 0; i < mat->size1; ++i) {
        for (size_t j = 0; j < mat->size2; ++j) {
            fprintf(stdout, "%Lf ", gsl_matrix_long_double_get(mat, i, j));
        }

        fprintf(stdout, "\n");
    }

    return 0;
}