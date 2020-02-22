#include <float.h>
#include "sampling.h"

int sampling_graph_iterative(double **out, coal_graph_node_t *graph, size_t reward_size) {
    double waiting_times[vector_length(graph->edges)];

    coal_graph_node_t *new_state = NULL;
    double min_value = FLT_MAX;

    for (size_t i = 0; i < vector_length(graph->edges); i++) {
        waiting_times[i] = dist_sample_exp(((weighted_edge_t*)(vector_get(graph->edges)))[i].weight);

        if (waiting_times[i] < min_value) {
            min_value = waiting_times[i];
            new_state = (coal_graph_node_t*)(((weighted_edge_t*)(vector_get(graph->edges)))[i].node);
        }
    }

    if (new_state != NULL) {
        for (size_t i = 0; i < reward_size; ++i) {
            (*out)[i] += min_value * graph->data.state[i];
        }

        return sampling_graph_iterative(out, new_state, reward_size);
    }

    return 0;
}

static void print_vector(vec_entry_t *v, size_t nmemb) {
    fprintf(stderr, "(");
    for (size_t i = 0; i < nmemb; i++) {
        fprintf(stderr, "%zu", v[i]);
    }
    fprintf(stderr, ")");
}

static void print_vector2(coal_graph_node_t **node, size_t nmemb) {
    fprintf(stderr, "(");
    for (size_t i = 0; i < nmemb; i++) {
        print_vector(node[i]->data.state, 4);
        fprintf(stderr, ", ");
    }
    fprintf(stderr, ")");
}

#define EPSILON 0.000001f

int _sampling_graph_pfd_constants(pdf_constant_t **out, coal_graph_node_t *graph, size_t reward_size,
        vector_t *path, vector_t *path_rates, double prob) {
    if (vector_length(graph->edges) == 0) {
        for (size_t i = 0; i < reward_size; ++i) {
            coal_graph_node_t **path_n = vector_get(path);
            double *path_r = vector_get(path_rates);

            for (size_t j = 0; j < vector_length(path); ++j) {
                coal_graph_node_t *u = path_n[j];

                if (u->data.state[i] == 0.0f) {
                    continue;
                }

                double u_r = path_r[j];
                double prod = 1.0f;

                for (size_t k = 0; k < vector_length(path); ++k) {
                    coal_graph_node_t *z = path_n[k];

                    if (k == j) {
                        continue;
                    }

                    if (z->data.state[i] == 0.0f) {
                        continue;
                    }

                    double z_r = path_r[k];
                    prod *= (z_r/z->data.state[i])/(z_r/z->data.state[i]-u_r/u->data.state[i]);
                }

                pdf_constant_t *entry = (*out + reward_size * u->data.vertex_index + i);
                entry->constant += prob*prod;
            }
        }
    } else {
        double rate = 0;
        weighted_edge_t * edges = vector_get(graph->edges);

        for (size_t i = 0; i < vector_length(graph->edges); i++) {
            rate += edges[i].weight;
        }

        // Our formula cannot handle the same rates,
        // so we add a small unique constant to each rate
        rate += EPSILON*vector_length(path);
        *((double*)vector_add(path_rates)) = rate;
        *((coal_graph_node_t**)vector_add(path)) = graph;

        // We set the rate here as well, it is redundant
        // for all paths, but it works.
        pdf_constant_t *entry = (*out + reward_size * graph->data.vertex_index);

        for (size_t j = 0; j < reward_size; ++j) {
            if (graph->data.state[j] == 0.0f) {
                entry[j].rate = 0.0f;
            } else {
                entry[j].rate = rate/graph->data.state[j];
            }
        }

        for (size_t i = 0; i < vector_length(graph->edges); i++) {
            _sampling_graph_pfd_constants(out,
                    (coal_graph_node_t *) edges[i].node, reward_size,
                    path, path_rates, prob * edges[i].weight / rate);
        }

        vector_remove_head(path);
        vector_remove_head(path_rates);
    }

    return 0;
}


int sampling_graph_pfd_constants(pdf_constant_t **out, size_t *out_size, coal_graph_node_t *graph, size_t reward_size) {
    size_t size;
    coal_label_vertex_index(&size, graph);
    size++;
    *out = calloc(size, sizeof(pdf_constant_t) * reward_size);
    vector_t *path;
    vector_t *path_rates;

    vector_init(&path, sizeof(coal_graph_node_t*), 8);
    vector_init(&path_rates, sizeof(double), 8);

    *out_size = size;

    return _sampling_graph_pfd_constants(out, graph, reward_size,
            path, path_rates, 1.0f);
}