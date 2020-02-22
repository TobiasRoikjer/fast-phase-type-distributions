

#include "coal.h"

int main(int argc, char **argv) {
    printf("type,n,i,j,value\n");
    for (size_t n = 1; n < atoi(argv[1]); n++) {
        coal_graph_node_t *graph;
        coal_gen_graph_reward(&graph, n, 0);

        for (size_t i = 0; i < n; i++) {
            printf("exp,%zu,%zu,%zu,%f\n", n, i, i, coal_mph_expected(graph, i));
        }

        for (size_t j = 0; j < n; j++) {
            for (size_t i = j; i < n; i++) {
                printf("cov,%zu,%zu,%zu,%f\n", n, j, i, coal_mph_cov(graph, i, j));
            }
        }
    }

    return 0;
}