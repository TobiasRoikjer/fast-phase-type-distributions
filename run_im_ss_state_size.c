#include "coal.h"

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Expected 1 argument got %i\n", argc-1);
        exit(1);
    }

    size_t n_tot = (size_t)atoi(argv[1]);

    fprintf(stderr, "Args:\n\tn (total):%zu\n\n", n_tot);

    fprintf(stdout, "n\tstate-type\tmig-type\tstates\tedges\n");


    for (size_t n = 0; n <= n_tot; ++n) {
        for (size_t mig_type = 0; mig_type <= 2; ++mig_type) {
            for (size_t sfs_type = 0; sfs_type <= 1; ++sfs_type) {
                size_t max_states = 0;
                size_t max_edges = 0;

                for (size_t n1 = 0; n1 < n; ++n1) {
                    size_t n2 = n - n1;

                    if (sfs_type == 0) {
                        for (size_t left1 = 0; left1 < n1 + n2; ++left1) {
                            for (size_t left2 = 0; left2 < n1 + n2 - left1; ++left2) {
                                if (left1 == 0 && left2 == 0) {
                                    continue;
                                }

                                coal_graph_node_t *graph;

                                coal_gen_im_cutoff_graph_args_t args = {
                                        .n1 = n1,
                                        .n2 = n2,
                                        .left_n1 = left1,
                                        .left_n2 = left2,
                                        .migration_type = mig_type,
                                        .pop_scale1 = 1,
                                        .pop_scale2 = 1,
                                        .mig_scale1 = 1,
                                        .mig_scale2 = 1
                                };

                                coal_graph_node_t *correct;
                                size_t largest_index;
                                coal_gen_im_prob_vertex_graph(&graph, &correct, args);
                                coal_label_vertex_index(&largest_index, graph);
                                max_states = largest_index > max_states ? largest_index : max_states;
                                size_t edges = coal_get_edges(graph);
                                max_edges = edges > max_edges ? edges : max_edges;
                            }
                        }
                    } else {
                        for (size_t coals = 0; coals < n; ++coals) {
                            coal_graph_node_t *graph;

                            coal_gen_im_graph_args_t args = {
                                    .n1 = n1,
                                    .n2 = n2,
                                    .num_iso_coal_events = coals,
                                    .migration_type = mig_type,
                                    .pop_scale1 = 1,
                                    .pop_scale2 = 1,
                                    .mig_scale1 = 1,
                                    .mig_scale2 = 1
                            };

                            size_t largest_index;
                            coal_gen_im_graph(&graph, NULL, args);
                            coal_label_vertex_index(&largest_index, graph);
                            max_states = largest_index > max_states ? largest_index : max_states;
                            size_t edges = coal_get_edges(graph);
                            max_edges = edges > max_edges ? edges : max_edges;
                        }
                    }
                }

                char *mig;

                switch (mig_type) {
                    case MIG_ALL:
                        mig = "ALL";
                        break;
                    case MIG_DIR:
                        mig = "DIRECTION";
                        break;
                    case MIG_ONCE:
                        mig = "ONCE";
                        break;
                    default:
                        mig = "NULL";
                        break;
                }

                fprintf(stdout, "%zu\t%s\t%s\t%zu\t%zu\n", n, sfs_type ? "JSFS" : "SS", mig, max_states, max_edges);
            }
        }
    }

    return 0;
}