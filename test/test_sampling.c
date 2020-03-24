#include <stdio.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>

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
    coal_gen_kingman_graph(&graph, n);

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
    coal_gen_kingman_graph(&graph, n);

    sampling_graph_pfd_constants(&constants, &size, graph, 4);

    pdf_constant_t *p = constants;

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            printf("%zu,%f,%f\n",j, p->constant, p->rate);
            p++;
        }
    }
}

static void test_sampling_constants_fast() {
    pdf_constant_t *constants;

    size_t size;
    size_t n = 4;
    coal_graph_node_t *graph;
    coal_gen_kingman_graph(&graph, n);

    sampling_graph_pfd_constants_rec(&constants, &size, graph, n, 0);

    pdf_constant_t *p = constants;

    for (size_t i = 0; i < size; ++i) {
        printf("%f,%f\n", p->constant, p->rate);
        p++;
    }
}

static void testpdf(double t) {
    pdf_constant_t *constants;

    size_t size;
    size_t n = 4;
    coal_graph_node_t *graph;
    coal_gen_kingman_graph(&graph, n);

    sampling_graph_pfd_constants_rec(&constants, &size, graph, n, 0);

    pdf_constant_t *p = constants;
    double pdf = 0;

    for (size_t i = 0; i < size; ++i) {
        if (!isnan(constants[i].rate) && !isinf(constants[i].constant)) {
            pdf += constants[i].constant * exp(-constants[i].rate * t);
        }
    }

    fprintf(stderr, "PDF (t=%f): %f", t, pdf);
}

int main(int argc, char **argv) {
//    test_randexp();
    //test_sampling_graph();
    //test_sampling_constants_slow();
    test_sampling_constants_fast();

    for (size_t i = 0; i < 10; ++i) {
        testpdf((double)i*0.1f);
    }

    return 0;
}