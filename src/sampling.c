#include <float.h>
#include <math.h>
#include "sampling.h"

int sampling_graph_iterative(sampling_number_t **out, coal_graph_node_t *graph, size_t reward_size) {
    sampling_number_t waiting_times[vector_length(graph->edges)];

    coal_graph_node_t *new_state = NULL;
    sampling_number_t min_value = FLT_MAX;

    for (size_t i = 0; i < vector_length(graph->edges); i++) {
        waiting_times[i] = dist_sample_exp(((weighted_edge_t*)(vector_get(graph->edges)))[i].weight);

        if (waiting_times[i] < min_value) {
            min_value = waiting_times[i];
            new_state = (coal_graph_node_t*)(((weighted_edge_t*)(vector_get(graph->edges)))[i].node);
        }
    }

    if (new_state != NULL) {
        for (size_t i = 0; i < reward_size; ++i) {
            (*out)[i] += min_value * graph->data.state_vec[i];
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
        print_vector(node[i]->data.state_vec, 4);
        fprintf(stderr, ", ");
    }
    fprintf(stderr, ")");
}

void print_array(sampling_number_t *array, size_t len) {
    fprintf(stderr, "array: ");
    for (size_t i = 0; i < len; ++i) {
        fprintf(stderr, "%Lf, ", array[i]);
    }
    fprintf(stderr, "\n");
}

#define EPSILON FLT_EPSILON*100
#define REWARD_INCREASE 0.001f;

int _sampling_graph_pfd_constants(pdf_constant_t **out, coal_graph_node_t *graph, size_t reward_size,
        vector_t *path, vector_t *path_rates, sampling_number_t prob) {
    if (vector_length(graph->edges) == 0) {
        for (size_t i = 0; i < reward_size; ++i) {
            coal_graph_node_t **path_n = vector_get(path);
            sampling_number_t *path_r = vector_get(path_rates);

            for (size_t j = 0; j < vector_length(path); ++j) {
                coal_graph_node_t *u = path_n[j];

                if (u->data.state_vec[i] == 0.0f) {
                    continue;
                }

                sampling_number_t u_r = path_r[j];
                sampling_number_t prod = 1.0f;

                for (size_t k = 0; k < vector_length(path); ++k) {
                    coal_graph_node_t *z = path_n[k];

                    if (k == j) {
                        continue;
                    }

                    if (z->data.state_vec[i] == 0.0f) {
                        continue;
                    }

                    sampling_number_t z_r = path_r[k];
                    prod *= (z_r/z->data.state_vec[i])/(z_r/z->data.state_vec[i]-u_r/u->data.state_vec[i]);
                }

                pdf_constant_t *entry = (*out + reward_size * u->data.vertex_index + i);
                entry->constant += prob*prod;
            }
        }
    } else {
        sampling_number_t rate = 0;
        weighted_edge_t * edges = vector_get(graph->edges);

        for (size_t i = 0; i < vector_length(graph->edges); i++) {
            rate += edges[i].weight;
        }

        // Our formula cannot handle the same rates,
        // so we add a small unique constant to each rate
        rate += EPSILON*vector_length(path);
        *((sampling_number_t*)vector_add(path_rates)) = rate;
        *((coal_graph_node_t**)vector_add(path)) = graph;

        // We set the rate here as well, it is redundant
        // for all paths, but it works.
        pdf_constant_t *entry = (*out + reward_size * graph->data.vertex_index);

        for (size_t j = 0; j < reward_size; ++j) {
            if (graph->data.state_vec[j] == 0.0f) {
                entry[j].rate = 0.0f;
            } else {
                entry[j].rate = rate/graph->data.state_vec[j];
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
    vector_init(&path_rates, sizeof(sampling_number_t), 8);

    *out_size = size;

    return _sampling_graph_pfd_constants(out, graph, reward_size,
            path, path_rates, 1.0f);
}

inline static sampling_number_t get_weight(weighted_edge_t value, size_t index) {
    // Add a small value such that all are unique
    return value.weight * (1 + index * 0.00001);
}

inline static sampling_number_t get_rate(coal_graph_node_t *node, size_t n_vertices) {
    sampling_number_t rate = 0;
    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        // Sum the rates
        rate += get_weight(values[i], n_vertices-1-((coal_graph_node_t*)values[i].node)->data.vertex_index);
    }

    return rate;
}

inline static sampling_number_t get_reward_rate(coal_graph_node_t *node, size_t n_vertices, size_t reward_index) {
    return get_rate(node, n_vertices)/node->data.state_vec[reward_index];
}

inline static bool zero_constants(sampling_number_t *array, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        if (fabsl(array[i]) > EPSILON) {
            return false;
        }
    }

    return true;
}

inline static sampling_number_t reward_scaled_diff(size_t a_index, size_t b_index, const sampling_number_t *reward_rates) {
    sampling_number_t rate_a = reward_rates[a_index];
    sampling_number_t rate_b = reward_rates[b_index];
    
    return rate_a / (rate_a - rate_b);
}

int _sampling_graph_pfd_constants_rec(sampling_number_t **k_out,
                                       sampling_number_t *reward_rates,
                                       coal_graph_node_t *vertex,
                                       size_t n_vertices,
                                       size_t reward_index,
                                       size_t vector_len,
                                       sampling_number_t **vertices_k) {
    if (!vertex->data.visited) {
        DIE_ERROR(1, "Should be visited is not \n");
    }
    if (vertices_k[vertex->data.vertex_index] != NULL) {
        *k_out = vertices_k[vertex->data.vertex_index];
        return 0;
    } else {
        *k_out = calloc((size_t)(vertex->data.vertex_index)+1, sizeof(sampling_number_t));
        sampling_number_t rate = get_rate(vertex, n_vertices);

        if (vector_length(vertex->edges) == 0) {
            return 0;
        } else if (vertex->data.state_vec[reward_index] == 0) {
            weighted_edge_t *values = vector_get(vertex->edges);

            for (size_t i = 0; i < vector_length(vertex->edges); i++) {
                coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;
                sampling_number_t *k_child;
                _sampling_graph_pfd_constants_rec(&k_child, reward_rates,
                                                  child, n_vertices, reward_index,
                                                  vector_len, vertices_k);

                sampling_number_t scale = get_weight(values[i], n_vertices - 1 -(((coal_graph_node_t*)values[i].node)->data.vertex_index)) / rate;

                for (size_t j = 0; j < child->data.vertex_index+1; ++j) {
                    (*k_out)[j] += k_child[j] * scale;
                }
            }

            vertices_k[vertex->data.vertex_index] = *k_out;

            return 0;
        } else {
            weighted_edge_t *values = vector_get(vertex->edges);

            for (size_t i = 0; i < vector_length(vertex->edges); i++) {
                coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;
                sampling_number_t scale = get_weight(values[i], n_vertices - 1 -(((coal_graph_node_t*)values[i].node)->data.vertex_index)) / rate;

                for (size_t j = 0; j < child->data.vertex_index+1; ++j) {
                    sampling_number_t *k_child;
                    _sampling_graph_pfd_constants_rec(&k_child, reward_rates,
                                                      (coal_graph_node_t *) values[i].node, n_vertices, reward_index,
                                                      vector_len, vertices_k);
                    if (fabsl(k_child[j]) <= EPSILON) {
                        continue;
                    }

                    (*k_out)[j] += scale * k_child[j] *
                                   reward_scaled_diff((size_t)vertex->data.vertex_index, j, reward_rates);
                }
            }

            (*k_out)[vertex->data.vertex_index] = 0;

            for (size_t i = 0; i < vector_length(vertex->edges); i++) {
                coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;
                sampling_number_t scale = get_weight(values[i], n_vertices - 1 -(((coal_graph_node_t*)values[i].node)->data.vertex_index)) / rate;
                sampling_number_t *k_child;
                _sampling_graph_pfd_constants_rec(&k_child, reward_rates,
                                                  (coal_graph_node_t *) values[i].node, n_vertices, reward_index,
                                                  vector_len, vertices_k);
                if (zero_constants(k_child, (size_t)(child->data.vertex_index)+1)) {
                    (*k_out)[vertex->data.vertex_index] += values[i].weight / rate;
                } else {
                    for (size_t j = 0; j < child->data.vertex_index+1; ++j) {
                        if (fabsl(k_child[j]) <= EPSILON) {
                            continue;
                        }

                        (*k_out)[vertex->data.vertex_index] += scale * k_child[j] *
                                                               reward_scaled_diff(j, (size_t)vertex->data.vertex_index, reward_rates);
                    }
                }
            }
        }

        vertices_k[vertex->data.vertex_index] = *k_out;

        return 0;
    }
}

int cmp_pdf_constant(const pdf_constant_t *a, const pdf_constant_t *b) {
    if (fabsl(a->rate) < fabsl(b->rate)) {
        return -1;
    } else if (fabsl(a->rate) > fabsl(b->rate)) {
        return 1;
    } else {
        return 0;
    }
}

int sampling_graph_pfd_constants_rec(pdf_constant_t **out, size_t *out_size,
                                     coal_graph_node_t *graph, size_t vector_len, size_t reward_index) {
    size_t size;
    coal_label_vertex_index(&size, graph);
    size++;
    coal_graph_reset_visited(graph);

    sampling_number_t *reward_rates = calloc(size, sizeof(sampling_number_t));

    queue_t *queue;
    queue_create(&queue, 8);
    queue_enqueue(queue, graph);
    size_t index = size - 1;

    while(!queue_empty(queue)) {
        coal_graph_node_t *node = queue_dequeue(queue);

        if (node->data.visited) {
            continue;
        }

        node->data.vertex_index = index;
        index--;
        reward_rates[node->data.vertex_index] = get_reward_rate(node, size, reward_index);

        weighted_edge_t *values = vector_get(node->edges);

        node->data.visited = true;

        for (size_t i = 0; i < vector_length(node->edges); i++) {
            queue_enqueue(queue, values[i].node);
        }
    }

    vec_entry_t *empty_vec = calloc(vector_len, sizeof(vec_entry_t));
    empty_vec[0] = (vec_entry_t) -1;


    sampling_number_t *k_out;
    sampling_number_t **vertices_k = calloc(size, sizeof(sampling_number_t*));

    _sampling_graph_pfd_constants_rec(&k_out,
                                      reward_rates, graph, size, reward_index,
                                      vector_len, vertices_k);

    *out = calloc(size, sizeof(pdf_constant_t));

    size_t count = 0;
    for (size_t j = 0; j < size; ++j) {
        if (!isinf(reward_rates[j]) && !isnan(reward_rates[j])) {
            (*out)[count].rate = reward_rates[j];
            (*out)[count].constant = k_out[j];
            count++;
        }
    }

    *out_size = count;

    qsort(*out, count, sizeof(pdf_constant_t),
          (int(*)(const void *, const void *)) &cmp_pdf_constant);

    return 0;
}


inline static sampling_number_t get_direct_rate(const coal_graph_node_t *node) {
    sampling_number_t rate = 0;
    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        // Sum the rates
        rate += values[i].weight;
    }

    return rate;
}

int _sampling_graph_pfd_constants_rec_rw(sampling_number_t **k_out,
                                      sampling_number_t *reward_rates,
                                      coal_graph_node_t *vertex,
                                      size_t n_vertices,
                                      sampling_number_t **vertices_k) {
    if (vertices_k[vertex->data.vertex_index] != NULL) {
        *k_out = vertices_k[vertex->data.vertex_index];
        return 0;
    } else {
        //fprintf(stderr, "NODE: %zu\n",
        //        vertex->data.vertex_index);
        *k_out = calloc((size_t) vertex->data.vertex_index+1, sizeof(sampling_number_t));
        sampling_number_t rate = get_direct_rate(vertex);

        if (vector_length(vertex->edges) == 0) {
            vertices_k[vertex->data.vertex_index] = *k_out;
            return 0;
        } else if (vertex->data.reward == 0) {
            //fprintf(stderr,"Reward zero (%zu):\n",
            //        vertex->data.vertex_index);
            weighted_edge_t *values = vector_get(vertex->edges);

            for (size_t i = 0; i < vector_length(vertex->edges); i++) {
                coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;
                sampling_number_t *k_child;
                _sampling_graph_pfd_constants_rec_rw(&k_child, reward_rates,
                                                  child, n_vertices,
                                                  vertices_k);

                sampling_number_t prob = values[i].weight / rate;

                for (size_t j = 0; j <= child->data.vertex_index; ++j) {
                    //fprintf(stderr,"RC (%zu): Incrementing my k[%zu] by childk=%Lf * prob %Lf\n",
                    //        vertex->data.vertex_index, j, k_child[j], prob);
                    (*k_out)[j] += k_child[j] * prob;
                }
            }

            vertices_k[vertex->data.vertex_index] = *k_out;

            /*fprintf(stderr, "I am %zu and my k is:\n", vertex->data.vertex_index);

            for (size_t j = 0; j < n_vertices; ++j) {
                fprintf(stderr, "%zu: %Lf\n", j, (*k_out)[j]);
                if (isnan((*k_out)[j]) || isinf((*k_out)[j])) {
                    DIE_ERROR(1, "Infinity or NaN\n");
                }
            }

            fprintf(stderr, "\n\n");*/

            return 0;
        } else {
            weighted_edge_t *values = vector_get(vertex->edges);

            //fprintf(stderr,"C (%zu): I have %zu edges\n",
            //        vertex->data.vertex_index, vector_length(vertex->edges));
            for (size_t i = 0; i < vector_length(vertex->edges); i++) {
                coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;
                sampling_number_t prob = values[i].weight / rate;

                for (size_t j = 0; j <= child->data.vertex_index; ++j) {
                    sampling_number_t *k_child;
                    _sampling_graph_pfd_constants_rec_rw(&k_child, reward_rates,
                                                      child, n_vertices,
                                                      vertices_k);
                    //fprintf(stderr,"C (%zu): For the edge to child %zu: Considering vertex %zu has k %Lf\n",
                    //        vertex->data.vertex_index, child->data.vertex_index, j, k_child[j]);
                    if (fabsl(k_child[j]) <= EPSILON) {
                        continue;
                    }

                    //fprintf(stderr, "C (%zu): Adding to %zu, += %Lf * %Lf * %Lf\n",
                    //        vertex->data.vertex_index,j, prob, k_child[j], reward_scaled_diff((size_t)vertex->data.vertex_index, j, reward_rates));
                    (*k_out)[j] += prob * k_child[j] *
                                   reward_scaled_diff((size_t)vertex->data.vertex_index, j, reward_rates);
                }
            }

            (*k_out)[vertex->data.vertex_index] = 0;

            for (size_t i = 0; i < vector_length(vertex->edges); i++) {
                coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;
                sampling_number_t prob = values[i].weight / rate;
                sampling_number_t *k_child;
                _sampling_graph_pfd_constants_rec_rw(&k_child, reward_rates,
                                                  (coal_graph_node_t *) values[i].node, n_vertices,
                                                  vertices_k);
                long double sum = 0;

                for (size_t j = 0; j <= child->data.vertex_index; ++j) {
                    if (fabsl(k_child[j]) <= EPSILON) {
                        continue;
                    }

                    //fprintf(stderr, "C (%zu): My OWN! Adding to %zu, += %Lf * %Lf * %Lf (which is %Lf)\n",
                    //       vertex->data.vertex_index, vertex->data.vertex_index, prob, k_child[j],
                    //        reward_scaled_diff(j, (size_t) vertex->data.vertex_index, reward_rates),
                    //        prob * k_child[j] *
                    //        reward_scaled_diff(j, (size_t) vertex->data.vertex_index, reward_rates));
                    (*k_out)[vertex->data.vertex_index] += prob * k_child[j] *
                                                           reward_scaled_diff(j,
                                                                              (size_t) vertex->data.vertex_index,
                                                                              reward_rates);

                    //fprintf(stderr, "after %Lf\n", (*k_out)[vertex->data.vertex_index]);
                    sum += k_child[j];
                }


                (*k_out)[vertex->data.vertex_index] += (1-sum) * prob;
            }
        }

        vertices_k[vertex->data.vertex_index] = *k_out;

        /*
        fprintf(stderr, "I am %zu and my k is:\n", vertex->data.vertex_index);

        for (size_t j = 0; j < n_vertices; ++j) {
            fprintf(stderr, "%zu: %Lf\n", j, (*k_out)[j]);
            if (isnan((*k_out)[j]) || isinf((*k_out)[j])) {
                DIE_ERROR(1, "Infinity or NaN\n");
            }
        }

        fprintf(stderr, "\n\n");
        */

        return 0;
    }
}

int dir = 1;

static bool node_is_descendant(coal_graph_node_t *node, coal_graph_node_t *find) {
    if (node == find) {
        return true;
    }

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        if (node_is_descendant((coal_graph_node_t *) values[i].node, find)) {
            return true;
        }
    }

    return false;
}

static long double get_reward_rate_unique(coal_graph_node_t *node, avl_double_node_t **bst) {
    //fprintf(stderr, "The node %zu  has a direct rate of %Lf and a reward of %f\n",
    //        node->data.vertex_index, get_direct_rate(node), node->data.reward);
    if (node->data.reward == 0) {
        return INFINITY;
    }

    long double rw_rate = get_direct_rate(node) / node->data.reward;
    avl_double_node_t *found;
    found = avl_double_find(*bst, rw_rate, EPSILON);

    if (found != NULL) {
        //fprintf(stderr, "Already found %Lf\n", rw_rate);
        node->data.reward += dir * REWARD_INCREASE;
        return get_reward_rate_unique(node, bst);
    }

    avl_double_insert(bst, rw_rate, node, EPSILON);

    return rw_rate;
}

void assign(coal_graph_node_t *node, size_t *index, avl_double_node_t **bst,
            sampling_number_t *reward_rates) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;

    node->data.vertex_index = *index;
    (*index)--;

    long double reward_rate = get_reward_rate_unique(node, bst);
    reward_rates[node->data.vertex_index] = reward_rate;
    dir = dir * -1;

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        assign((coal_graph_node_t *) values[i].node, index, bst, reward_rates);
    }
}

int sampling_graph_pfd_constants_rec_rw(pdf_constant_t **out, size_t *out_size,
                                     coal_graph_node_t *graph) {
    size_t max_index;
    coal_label_topological_index(&max_index, graph);
    size_t length = max_index+1;
    coal_graph_reset_visited(graph);

    sampling_number_t *reward_rates = calloc(length, sizeof(sampling_number_t));
    coal_graph_node_t **vertices = calloc(length, sizeof(coal_graph_node_t *));

    avl_double_node_t *bst = NULL;

    queue_t *queue;
    queue_create(&queue, 8);
    queue_enqueue(queue, graph);

    while(!queue_empty(queue)) {
        coal_graph_node_t *node = queue_dequeue(queue);
        if (node->data.visited) {
            continue;
        }
        node->data.vertex_index = max_index - node->data.vertex_index;

        reward_rates[node->data.vertex_index] = get_reward_rate_unique(node, &bst);
        vertices[node->data.vertex_index] = node;
        dir = dir * -1;

        weighted_edge_t *values = vector_get(node->edges);
        node->data.visited = true;
        for (size_t i = 0; i < vector_length(node->edges); i++) {
            queue_enqueue(queue, values[i].node);
        }
    }

    sampling_number_t *k_out;
    sampling_number_t **vertices_k = calloc(length, sizeof(sampling_number_t*));

    _sampling_graph_pfd_constants_rec_rw(&k_out,
                                      reward_rates, graph, length,
                                      vertices_k);

    *out = calloc(length, sizeof(pdf_constant_t));

    size_t count = 0;
    for (size_t j = 0; j < length; ++j) {
        if (!isinf(reward_rates[j]) && !isnan(reward_rates[j])) {
            (*out)[count].rate = reward_rates[j];
            (*out)[count].constant = k_out[j];
            (*out)[count].node = vertices[j];
            count++;
        }
    }

    *out_size = count;

    qsort(*out, count, sizeof(pdf_constant_t),
          (int(*)(const void *, const void *)) &cmp_pdf_constant);

    return 0;
}