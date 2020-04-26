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
        printf("%Lf,", dist_sample_exp(0.1));
    }
}

static void test_sampling_graph() {
    size_t n = 4;
    coal_graph_node_t *graph;
    coal_gen_kingman_graph(&graph, n);

    dist_sampling_set_random_seed((unsigned int)time(NULL));

    for (size_t i = 0; i < 10000; ++i) {
        weight_t *times = calloc(n, sizeof(double));
        sampling_graph_iterative(&times, graph, n);

        for (size_t j = 0; j < n; ++j) {
            printf("%Lf,", times[j]);
        }
        printf("\n");
    }
}

static void test_sampling_constants_slow() {
    pdf_constant_t *constants;

    size_t size;
    size_t n = 5;
    coal_graph_node_t *graph;
    coal_gen_kingman_graph(&graph, n);

    sampling_graph_pfd_constants(&constants, &size, graph, 5);

    pdf_constant_t *p = constants;

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            printf("%zu,%Lf,%Lf\n",j, p->constant, p->rate);
            p++;
        }
    }
}


static void test_sampling_constants_slow2() {
    pdf_constant_t *constants;

    size_t size;
    size_t n = 5;
    coal_graph_node_t *graph;
    coal_gen_kingman_graph(&graph, n);

    sampling_graph_pfd_constants(&constants, &size, graph, 5);

    pdf_constant_t *p = constants;

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            if (j == 1) {
                printf("%Lf,%Lf\n", p->constant, p->rate);
            }
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

    sampling_graph_pfd_constants_rec(&constants, &size, graph, n, 2);

    pdf_constant_t *p = constants;

    for (size_t i = 0; i < size; ++i) {
        printf("%Lf,%Lf\n", p->constant, p->rate);
        p++;
    }
}

static void testcdffast(double t) {
    pdf_constant_t *constants;

    size_t size;
    size_t n = 6;
    coal_graph_node_t *graph;
    coal_gen_kingman_graph(&graph, n);

    sampling_graph_pfd_constants_rec(&constants, &size, graph, n, 0);

    pdf_constant_t *p = constants;
    double cdf = 0;

    for (size_t i = 0; i < size; ++i) {
        if (!isnan(constants[i].rate) && !isinf(constants[i].constant)) {
            cdf += constants[i].constant * expl(-constants[i].rate * t);
        }
    }

    fprintf(stderr, "CDF (t=%f): %f\n", t, cdf);
}


size_t n;

double reward_sing(coal_graph_node_t *node) {
    return node->data.state_vec[1];
}

double reward_singinc(coal_graph_node_t *node) {
    return node->data.state_vec[1]+0.001;
}

double reward_length(coal_graph_node_t *node) {
    size_t lineages = 0;

    for (size_t i = 0; i < n; ++i) {
        lineages += node->data.state_vec[i];
    }

    return lineages;
}

double reward_one(coal_graph_node_t *node) {
    return 1;
}

double reward_two(coal_graph_node_t *node) {
    return 2;
}

void test_rw_cdf() {
    pdf_constant_t *constants;
    n = 8;

    coal_graph_node_t *graph2;
    coal_graph_node_t *graph2s;
    coal_gen_kingman_graph(&graph2, n);
    coal_label_vertex_index(NULL, graph2);
    coal_rewards_set(graph2, reward_sing);
    graph2->data.reward = 0;
    coal_reward_transform(graph2, &graph2s);
    coal_rewards_set(graph2s, reward_one);
    graph2s->data.reward = 0;
    size_t constants_size;
    //coal_graph_node_t *first = (coal_graph_node_t *) ((weighted_edge_t*)vector_get(graph2s->edges))[0].node;
    sampling_graph_pfd_constants_rec_rw(&constants, &constants_size, graph2s);
fflush(stderr);
    FILE *f = fopen("/mnt/c/Users/Tobias/Documents/school/speciale/data/consts.tsv","w");
    for (size_t i = 0; i < constants_size; ++i) {
        fprintf(f, "%Lf %Lf\n", constants[i].constant, constants[i].rate);
        fprintf(stdout, "%Lf %Lf\n", constants[i].constant, constants[i].rate);
    }

    fclose(f);

    for (double t = 0; t < 1; t+=0.1f) {
        long double res = 0;
        long double respdf = 0;
        long double constsum = 0;
        for (size_t i = 0; i < constants_size; ++i) {
            res += constants[i].constant * (1 - expl(-constants[i].rate * t));
            respdf += constants[i].constant *constants[i].rate* expl(-constants[i].rate * t);
            constsum += constants[i].constant;
        }

        fprintf(stdout, "CDF %f: %Lf %Lf (%Lf)\n", t, res, respdf, 1 - constsum);
    }
}

void test_rw_cdf2() {
    pdf_constant_t *constants;
    n = 8;

    coal_graph_node_t *graph2;
    coal_gen_kingman_graph(&graph2, n);
    coal_label_vertex_index(NULL, graph2);
    coal_rewards_set(graph2, reward_singinc);
    graph2->data.reward = 0;
    size_t constants_size;
    sampling_graph_pfd_constants_rec_rw(&constants, &constants_size, graph2);
    fflush(stderr);
    FILE *f = fopen("/mnt/c/Users/Tobias/Documents/school/speciale/data/consts.tsv","w");
    for (size_t i = 0; i < constants_size; ++i) {
        fprintf(f, "%Lf %Lf\n", constants[i].constant, constants[i].rate);
        fprintf(stdout, "%Lf %Lf\n", constants[i].constant, constants[i].rate);
    }

    fclose(f);

    for (double t = 0; t < 1; t+=0.1f) {
        long double res = 0;
        long double respdf = 0;
        long double constsum = 0;
        for (size_t i = 0; i < constants_size; ++i) {
            res += constants[i].constant * (1 - expl(-constants[i].rate * t));
            respdf += constants[i].constant *constants[i].rate* expl(-constants[i].rate * t);
            constsum += constants[i].constant;
        }

        fprintf(stdout, "CDF %f: %Lf %Lf (%Lf)\n", t, res, respdf, 1 - constsum);
    }
}

int main(int argc, char **argv) {
//    test_randexp();
    //test_sampling_graph();
    //test_sampling_constants_slow2();
    //printf("\n\n");
    //test_sampling_constants_fast();

    /*for (size_t i = 0; i < 10; ++i) {
        testcdffast((double)i*0.1f);
    }*/
    test_rw_cdf();
    //test_rw_cdf2();
    //test_sampling_constants_slow();

    return 0;
}