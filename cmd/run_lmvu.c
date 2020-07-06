#include "../src/coal.h"

int main(int argc, char **argv) {
    printf("type,n,i,j,value\n");

    for (size_t n = (size_t)atoi(argv[1]); n <= atoi(argv[2]); n++) {
        coal_graph_node_t *graph;
        coal_gen_kingman_graph(&graph, n);
        double ***cov;
        cov = coal_mph_cov_all(graph, n);

        for (size_t i = 0; i < n; i++) {
            printf("exp,%zu,%zu,%zu,%Lf\n", n, i, i, coal_mph_expected(graph, i));
        }

        for (size_t j = 0; j < n; j++) {
            for (size_t i = 0; i <= j; i++) {
                printf("cov,%zu,%zu,%zu,%f\n", n, j, i, (*cov)[j][i]);
            }
        }

        fprintf(stderr, "Done %zu\n", n);

        fflush(stdout);
    }

    return 0;
}