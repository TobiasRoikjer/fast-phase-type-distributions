#include "coal.h"
#include "sampling.h"
#include <string.h>
#include <math.h>

sampling_number_t cdf(pdf_constant_t *constants, size_t size, sampling_number_t t) {
    sampling_number_t cdf = 0;

    for (size_t i = 0; i < size; ++i) {
        //fprintf(stderr, "%Lf, %Lf\n",constants[i].constant, constants[i].rate);
        if (!isnan(constants[i].constant) && !isinf(constants[i].constant) &&
                !isnan(constants[i].rate) && !isinf(constants[i].rate)) {
            if (constants[i].constant != 0) {
                cdf += constants[i].constant*expl(-constants[i].rate*t);
            }
        }
    }

    return 1-cdf;
}

int main(int argc, char **argv) {
    size_t n = (size_t) atoi(argv[1]);
    size_t r = (size_t) atoi(argv[2]);
    coal_graph_node_t *graph;
    coal_gen_kingman_graph(&graph, n);

    pdf_constant_t *constants;
    size_t size;

    sampling_graph_pfd_constants_rec(&constants, &size, graph, n, r-1);
    fprintf(stdout, "t,cdf\n");
    for (size_t j = 0; j < 100; ++j) {
        fflush(stderr);
        fprintf(stdout, "%f,%Lf\n", j*0.1f, cdf(constants, size, (sampling_number_t)j*0.1f+0.000001f));
    }

    return 0;
}