#include "coal.h"
#include "sampling.h"
#include <time.h>
#include <string.h>

int main(int argc, char **argv) {
    if (strcmp(argv[1], "pdf_constants") == 0) {
        size_t n = (size_t) atoi(argv[2]);
        coal_graph_node_t *graph;
        coal_gen_graph_reward(&graph, n, 0);
        fprintf(stderr, "Done generating graph\n");

        time_t start;
        start = time(NULL);

        pdf_constant_t *constants;
        size_t size;

        sampling_graph_pfd_constants_rec(&constants, &size, graph, n, (size_t) atoi(argv[3]) - 1);

        fprintf(stdout, "Time elapsed for %zu samples rewarded by %zu'tons: %lu\n",
                (size_t) atoi(argv[2]), (size_t) atoi(argv[3]), (long) (time(NULL)-start));

        return 0;
    }

    return 0;
}