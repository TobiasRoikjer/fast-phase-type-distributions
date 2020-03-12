#include "coal.h"

int main(int argc, char **argv) {
    coal_graph_node_t *graph;
    coal_gen_graph_reward(&graph, (size_t)atoi(argv[1]), 0);

    return 0;
}