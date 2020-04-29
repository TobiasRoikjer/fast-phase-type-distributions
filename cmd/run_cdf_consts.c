#include "../src/coal.h"
#include "../src/sampling.h"
#include <string.h>
#include <math.h>

sampling_number_t cdf(pdf_constant_t *constants, size_t size, sampling_number_t t) {
    sampling_number_t cdf = 0;
    sampling_number_t sum_constants = 0;

    for (size_t i = 0; i < size; ++i) {
        cdf += constants[i].constant*expl(-constants[i].rate*t);
        sum_constants += constants[i].constant;
    }

    return 1-cdf;
}

size_t r;

double reward_by(coal_graph_node_t *node) {
    return node->data.state_vec[r-1];
}

double reward_one(coal_graph_node_t *node) {
    return 1;
}

int main(int argc, char **argv) {
    size_t n = (size_t) atoi(argv[1]);
    r = (size_t) atoi(argv[2]);
    coal_graph_node_t *graph;
    coal_gen_kingman_graph(&graph, n);
    coal_rewards_set(graph, reward_by);
    pdf_constant_t *constants;
    size_t size;

    sampling_graph_pfd_constants_rec_rw(&constants, &size, graph);
    long double const_sum = 0;
    for (size_t i = 0; i < size; ++i) {
        const_sum += constants[i].constant/1024;
    }

    fprintf(stdout, "t,cdf\n");
    for (size_t j = 0; j < 100; ++j) {
        fflush(stderr);
        fprintf(stdout, "%f,%Lf\n", j*0.1f, cdf(constants, size, (sampling_number_t)j*0.1f));
    }

    return 0;
}

