#include "../src/coal.h"
#include "../src/sampling.h"
#include <time.h>
#include <string.h>

size_t reward_index;

double reward(coal_graph_node_t *node) {
    return node->data.state_vec[reward_index];
}

double reward_one(coal_graph_node_t *node) {
    return 1;
}

int main(int argc, char **argv) {
    if (strcmp(argv[1], "ss_kingman") == 0) {
        size_t n = (size_t) atoi(argv[2]);
        coal_graph_node_t *graph;

        time_t start;
        start = time(NULL);
        coal_gen_kingman_graph(&graph, n);

        fprintf(stdout, "Time elapsed for %zu samples: %lu\n",
                (size_t) atoi(argv[2]), (long) (time(NULL)-start));

        return 0;
    } else if (strcmp(argv[1], "pdf_constants") == 0) {
        size_t n = (size_t) atoi(argv[2]);
        reward_index = (size_t) atoi(argv[3]) - 1;
        coal_graph_node_t *graph;
        coal_gen_kingman_graph(&graph, n);
        coal_rewards_set(graph, reward);
        fprintf(stderr, "Done generating graph\n");

        time_t start;
        start = time(NULL);

        pdf_constant_t *constants;
        size_t size;

        sampling_graph_pfd_constants_rec_rw(&constants, &size, graph);

        fprintf(stdout, "Time elapsed for %zu samples rewarded by %zu'tons: %lu\n",
                (size_t) atoi(argv[2]), (size_t) atoi(argv[3]), (long) (time(NULL)-start));

        return 0;
    } else if (strcmp(argv[1], "pdf_constants_trans") == 0) {
        size_t n = (size_t) atoi(argv[2]);
        reward_index = (size_t) atoi(argv[3]) - 1;
        coal_graph_node_t *graph;
        coal_graph_node_t *start_node;
        coal_gen_kingman_graph(&graph, n);
        coal_rewards_set(graph, reward);
        fprintf(stderr, "Done generating graph\n");

        time_t start;
        start = time(NULL);

        pdf_constant_t *constants;
        size_t size;

        coal_reward_transform(graph, &start_node);
        coal_rewards_set(start_node, reward_one);
        start_node->data.reward = 0;

        sampling_graph_pfd_constants_rec_rw(&constants, &size, start_node);

        fprintf(stdout, "Time elapsed for %zu samples rewarded by %zu'tons: %lu\n",
                (size_t) atoi(argv[2]), (size_t) atoi(argv[3]), (long) (time(NULL)-start));

        return 0;
    } else if (strcmp(argv[1], "prob_coals") == 0) {
        size_t n = (size_t) atoi(argv[2]);
        size_t allow_back = (bool) atoi(argv[3]);

        time_t start;
        start = time(NULL);

        coal_gen_im_graph_args_t args = {
                .n1 = n,
                .n2 = n,
                .migration_type = allow_back,
                .pop_scale1 = 0.3f,
                .pop_scale2 = 0.7f,
                .mig_scale1 = 0.01f,
                .mig_scale2 = 0.01f,
        };

        long double* probs;
        coal_im_get_number_coals_probs(&probs, 1.0f, &args);

        fprintf(stdout, "Time elapsed for %zu+%zu samples: %lu\n",
                n, n, (long) (time(NULL)-start));

        return 0;
    }

    return 0;
}