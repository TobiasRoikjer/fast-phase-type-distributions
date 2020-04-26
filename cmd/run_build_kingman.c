#include "../src/coal.h"

int main(int argc, char **argv) {
    coal_graph_node_t *graph;
    coal_gen_kingman_graph(&graph, (size_t) atoi(argv[1]));

    return 0;
}