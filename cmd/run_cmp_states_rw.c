#include "../src/coal.h"

int main(int argc, char **argv) {
    printf("n,reward,full,reduced\n");

    size_t n = atoi(argv[1]);

    coal_graph_node_t *graph;
    size_t sizeT;
    coal_gen_kingman_graph(&graph, n);
    coal_label_vertex_index(&sizeT, graph);

    for (size_t rw = 0; rw < n; rw++) {
        coal_gen_kingman_graph_rw(&graph, n, rw);
        size_t size;
        coal_label_vertex_index(&size, graph);
        fprintf(stderr, "%zu,%zu,%zu,%zu\n", n, rw, sizeT, size);
    }

    return 0;
}