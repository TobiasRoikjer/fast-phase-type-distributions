#include "coal.h"
#include "sampling.h"
#include <string.h>
#include <math.h>

sampling_number_t cdf(pdf_constant_t *constants, size_t size, sampling_number_t t) {
    sampling_number_t cdf = 0;

    for (size_t i = 0; i < size; ++i) {
        cdf += constants[i].constant*constants[i].rate*(1-expl(-constants[i].rate*t));
    }

    return cdf;
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
    coal_graph_node_t *graph, *start;
    coal_gen_kingman_graph(&graph, n);
    coal_rewards_set(graph, reward_by);
    graph->data.reward = 0;
    coal_reward_transform(graph, &start);
    coal_rewards_set(start,reward_one);
    start->data.reward = 0;

    pdf_constant_t *constants;
    size_t size;
    fprintf(stdout, "%Lf\n", ((weighted_edge_t*)vector_get(start->edges))[0].weight);

    sampling_graph_pfd_constants_rec_rw(&constants, &size, start);
    long double const_sum = 0;
    for (size_t i = 0; i < size; ++i) {
        const_sum += constants[i].constant;
        //constants[i].constant *= 0;//((weighted_edge_t*)vector_get(start->edges))[0].weight;
    }
    fprintf(stdout, "%Lf\n", const_sum);

    fprintf(stdout, "t,cdf\n");
    for (size_t j = 0; j < 100; ++j) {
        fflush(stderr);
        fprintf(stdout, "%f,%Lf\n", j*0.1f, cdf(constants, size, (sampling_number_t)j*0.1f));
    }

    return 0;
}