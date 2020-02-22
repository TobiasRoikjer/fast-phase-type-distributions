#include <stdio.h>
#include <inttypes.h>
#include <time.h>

#include "../dist.h"
#include "../coal.h"
#include "../sampling.h"

static void test_rand01() {
    dist_sampling_set_random_seed((unsigned int)time(NULL));
    for (int i = 0; i < 400; i++) {
        printf("%f,", dist_rand01());
    }
}

static void test_randexp() {
    dist_sampling_set_random_seed((unsigned int)time(NULL));
    for (int i = 0; i < 400; i++) {
        printf("%f,", dist_sample_exp(0.1));
    }
}

static void test_sampling_graph() {
    size_t n = 4;
    coal_graph_node_t *graph;
    coal_gen_graph_reward(&graph, n, 0);

    dist_sampling_set_random_seed((unsigned int)time(NULL));

    for (size_t i = 0; i < 10000; ++i) {
        double *times = calloc(n, sizeof(double));
        sampling_graph_iterative(&times, graph, n);

        for (size_t j = 0; j < n; ++j) {
            printf("%f,", times[j]);
        }
        printf("\n");
    }
}

static void test_sampling_constants_slow() {
    pdf_constant_t *constants;

    size_t size;
    size_t n = 4;
    coal_graph_node_t *graph;
    coal_gen_graph_reward(&graph, n, 0);

    sampling_graph_pfd_constants(&constants, &size, graph, 4);

    pdf_constant_t *p = constants;

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            printf("%zu,%f,%f\n",j, p->constant, p->rate);
            p++;
        }
    }
}

int main(int argc, char **argv) {
//    test_randexp();
    //test_sampling_graph();
    test_sampling_constants_slow();

    return 0;
}