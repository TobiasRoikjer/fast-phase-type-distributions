#include "coal.h"
#include "sampling.h"
#include <string.h>
#include <math.h>

sampling_number_t cdf(pdf_constant_t *constants, size_t size, sampling_number_t t) {
    sampling_number_t cdf = 0;

    for (size_t i = 0; i < size; ++i) {
        if (!isnan(constants[i].rate) && !isinf(constants[i].rate)) {
            bool negative = (constants[i].constant < 0);
            sampling_number_t k = constants[i].constant;
            int sign = 1;
            if (negative) {
                k = -k;
                sign = -1;
            }
            /*
            fprintf(stderr, "const %Lf\n", constants[i].constant);
            fprintf(stderr, "powl %Lf\n", powl(k, 1/t) );
            fprintf(stderr, "expl %Lf\n", expl(-constants[i].rate));
            fprintf(stderr, "tot %Lf\n", powl((powl(k, 1/t)  * expl(-constants[i].rate)), t));
            */
            //cdf += sign*powl((powl(k, -1/t)  * expl(constants[i].rate)), -t);
            cdf += constants[i].constant*expl(-constants[i].rate*t);
        }
    }

    return 1-cdf;
}

int main(int argc, char **argv) {
    size_t n = (size_t) atoi(argv[1]);
    size_t r = (size_t) atoi(argv[2]);
    coal_graph_node_t *graph;
    coal_gen_graph_reward(&graph, n, 0);

    pdf_constant_t *constants;
    size_t size;

    sampling_graph_pfd_constants_rec(&constants, &size, graph, n, r-1);
    fprintf(stdout, "t,cdf\n");
    for (size_t j = 10; j < 100; ++j) {
        fflush(stderr);
        fprintf(stdout, "%f,%Lf\n", j*0.1f, cdf(constants, size, (sampling_number_t)j*0.1f+0.000001f));
        break;
    }

    return 0;
}