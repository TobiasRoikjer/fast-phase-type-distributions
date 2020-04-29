#include "../src/coal.h"

size_t reward_index;

double reward(coal_graph_node_t *node) {
    return node->data.state_vec[reward_index];
}

int main(int argc, char **argv) {
    coal_graph_node_t *graph;
    coal_gen_kingman_graph(&graph, (size_t) atoi(argv[1]));

    if (argc == 3) {
        reward_index = (size_t) atoi(argv[2]);
        coal_graph_node_t *start;
        coal_rewards_set(graph, reward);
        coal_reward_transform(graph, &start);
        graph = start;
    }

    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, graph);

    for (size_t i = 2; i < size; ++i) {
        for (size_t j = 2; j < size; ++j) {
            fprintf(stdout, "%Lf ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }

    return 0;
}