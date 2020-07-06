#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <float.h>
#include <gsl/gsl_matrix_long_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector_double.h>
#include "coal.h"

int gsl_vector_axpby(const double alpha, const gsl_vector * x, const double beta, gsl_vector * y) {
    gsl_vector *clone = gsl_vector_alloc(x->size);
    gsl_vector_memcpy(clone, x);
    gsl_vector_scale(clone, alpha);
    gsl_vector_add(y, clone);
    gsl_vector_free(clone);

    return 0;
}

typedef enum  {
  IM_VERTEX_TYPE_MIG,
  IM_VERTEX_TYPE_COAL,
  IM_VERTEX_TYPE_START
} im_vertex_type_t;

static void reset_graph_visited(coal_graph_node_t *node);

static void print_vector(FILE *stream,vec_entry_t *v, size_t nmemb) {
    fprintf(stream, "(");
    for (size_t i = 0; i < nmemb; i++) {
        fprintf(stream, "%zu", v[i]);
    }
    fprintf(stream, ")");
}

static void print_vector_spacing(FILE *stream,vec_entry_t *v, size_t nmemb, size_t spacing) {
    if (v == NULL) {
        fprintf(stream, "(NULL)");
        return;
    }

    fprintf(stream, "(");

    for (size_t i = 0; i < nmemb; i++) {
        if (i % spacing == 0) {
            fprintf(stream, "-");
        }

        fprintf(stream, "%zu", v[i]);
    }
    fprintf(stream, ")");
}

static void _get_abs_vertex(coal_graph_node_t **abs_vertex, coal_graph_node_t *graph) {
    if (graph->data.visited) {
        return;
    }

    graph->data.visited = true;

    weighted_edge_t *values = vector_get(graph->edges);

    if (vector_length(graph->edges) == 0) {
        if (*abs_vertex == NULL || *abs_vertex == graph) {
            *abs_vertex = graph;
        } else {
            DIE_ERROR(1, "Found multiple absorbing vertices\n");
        }
    }

    for (size_t i = 0; i < vector_length(graph->edges); i++) {
        _get_abs_vertex(abs_vertex, (coal_graph_node_t*)values[i].node);
    }
}

static coal_graph_node_t *get_abs_vertex(coal_graph_node_t *graph) {
    reset_graph_visited(graph);
    coal_graph_node_t *abs_vertex = NULL;
    _get_abs_vertex(&abs_vertex, graph);

    return abs_vertex;
}

/*
 * Note: Works only on DAGs and is slow, used for debugging
 */
static void ensure_graph_sanity(coal_graph_node_t *node) {
    if (node == NULL) {
        DIE_ERROR(1, "Node is NULL");
    }

    if (node->edges == NULL) {
        DIE_ERROR(1, "Node %zu has NULL edges",
                  node->data.vertex_index);
    }

    if (node->reverse_edges == NULL) {
        DIE_ERROR(1, "Node %zu has NULL reverse edges",
                  node->data.vertex_index);
    }

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); ++i) {
        for (size_t j = i+1; j < vector_length(node->edges); ++j) {
            if (values[i].node == values[j].node) {
                DIE_ERROR(1, "Graph has two edges. Node %zu to %zu",
                          node->data.vertex_index,
                          ((coal_graph_node_t*)values[i].node)->data.vertex_index);
            }
        }

        coal_graph_node_t *child = (coal_graph_node_t*)(values[i].node);

        ensure_graph_sanity(child);
    }
}

static void coal_graph_node_create(coal_graph_node_t **out,
        vec_entry_t *state_vec, void *state) {
    graph_node_create((graph_node_t**)out, sizeof(coal_graph_node_data_t));
    (*out)->data.state_vec = state_vec;
    (*out)->data.state = state;
    (*out)->data.vertex_index = -1;
    (*out)->data.reset_int = 0;
}

static int kingman_visit_vertex(coal_graph_node_t **out,
                        vec_entry_t *state,
                        avl_vec_node_t *bst,
                        size_t state_size,
                        size_t vector_length,
                        const size_t vec_nmemb) {
    avl_vec_node_t *bst_node = avl_vec_find(bst, state, vector_length);

    if (bst_node != NULL) {
        *out = bst_node->entry;
        return 0;
    } else {
        vec_entry_t *vertex_state = malloc(sizeof(vec_entry_t) * vector_length);
        memcpy(vertex_state, state, sizeof(vec_entry_t) * vector_length);

        coal_graph_node_create(out, vertex_state, vertex_state);

        avl_vec_insert(&bst, vertex_state, *out, vector_length);

        vec_entry_t *v = state;

        for (vec_entry_t i = 0; i < state_size; i++) {
            for (vec_entry_t j = i; j < state_size; j++) {
                if (((i == j && v[i] >= 2) || (i != j && v[i] > 0 && v[j] > 0))) {
                    mat_entry_t t = i == j ? v[i] * (v[i] - 1) / 2 : v[i] * v[j];

                    v[i]--;
                    v[j]--;
                    v[(i + j + 2) - 1]++;

                    coal_graph_node_t *new_vertex;
                    kingman_visit_vertex(&new_vertex,
                            v, bst,
                            state_size,
                            vector_length,
                            vec_nmemb);

                    v[i]++;
                    v[j]++;
                    v[(i + j + 2) - 1]--;

                    graph_add_edge((graph_node_t*)*out, (graph_node_t*)new_vertex, t);
                }
            }
        }

        return 0;
    }
}

static int kingman_visit_vertex_rw(coal_graph_node_t **out,
                                   size_t rw_index,
                                   vec_entry_t *state,
                                   avl_vec_node_t *bst,
                                   size_t state_size,
                                   size_t vector_length,
                                   const size_t vec_nmemb) {
    avl_vec_node_t *bst_node = avl_vec_find(bst, state, vector_length);

    if (bst_node != NULL) {
        *out = bst_node->entry;
        return 0;
    } else {
        /*fprintf(stderr, "ADDING: ");
        bool all_zero = true;
        for (size_t k = 0; k < state_size; ++k) {
            fprintf(stderr, "%zu, ", state[k]);

            if (state[k] > 0) {
                all_zero = false;
            }
        }

        if (all_zero) {
            DIE_ERROR(1, "ZERO");
        }
        fprintf(stderr, "\n");*/
        vec_entry_t *vertex_state = malloc(sizeof(vec_entry_t) * vector_length);
        memcpy(vertex_state, state, sizeof(vec_entry_t) * vector_length);

        coal_graph_node_create(out, vertex_state, vertex_state);

        avl_vec_insert(&bst, vertex_state, *out, vector_length);

        vec_entry_t *v = state;

        for (vec_entry_t i = 0; i < state_size; i++) {
            for (vec_entry_t j = i; j < state_size; j++) {
                if (((i == j && v[i] >= 2) || (i != j && v[i] > 0 && v[j] > 0))) {
                    mat_entry_t t = i == j ? v[i] * (v[i] - 1) / 2 : v[i] * v[j];

                    vec_entry_t *old = malloc(sizeof(vec_entry_t) * vector_length);
                    memcpy(old, v, sizeof(vec_entry_t) * vector_length);

                    v[i]--;
                    v[j]--;

                    if (i == state_size - 1 || j == state_size - 1) {
                        v[state_size - 1]++;
                        //fprintf(stderr, "PUSHING from %zu <> %zu, because %zu\n",
                        //        i, j, state_size - 1);
                    } else {
                        v[(i + j + 2) - 1]++;
                    }

                    bool *is_available = calloc(vector_length, sizeof(bool));

                    for (size_t k = 0; k < state_size; ++k) {
                        if (v[k] > 0) {
                            is_available[k] = true;
                            continue;
                        }

                        for (size_t k1 = 0; k1 < k; ++k1) {
                            size_t oppo = k - k1 - 1;

                            if (oppo == k1) {
                                //if (v[k1] >= 2) {
                                if (is_available[k1]) {
                                    is_available[k] = true;
                                }
                            } else {
                                //if (v[k1] > 0 && v[oppo] > 0) {
                                if (is_available[k1] && is_available[oppo]) {
                                    is_available[k] = true;
                                }
                            }
                        }
                    }

                    /*fprintf(stderr, "Vector: ");
                    for (size_t k = 0; k < state_size; ++k) {
                        fprintf(stderr, "%zu, ", v[k]);
                    }

                    fprintf(stderr, "\n");

                    fprintf(stderr, "isaval: ");
                    for (size_t k = 0; k < state_size; ++k) {
                        fprintf(stderr, "%i, ", is_available[k]);
                    }

                    fprintf(stderr, "\n");*/

                    for (ssize_t k = state_size - 2; k >= 0; --k) {
                        size_t oppo = rw_index - k - 1;

                        if (k > rw_index ||
                                (v[k] != 0 && k != rw_index && !is_available[oppo])) {
                            //fprintf(stderr, "should push all %zu\n", k);
                            v[state_size - 1] += v[k];
                            v[k] = 0;
                        }
                    }


                    coal_graph_node_t *new_vertex;
                    kingman_visit_vertex_rw(&new_vertex,
                                         rw_index,
                                         v, bst,
                                         state_size,
                                         vector_length,
                                         vec_nmemb);

                    v = old;

                    graph_add_edge((graph_node_t*)*out, (graph_node_t*)new_vertex, t);
                }
            }
        }

        return 0;
    }
}

static void print_graph_node(coal_graph_node_t *node, size_t vec_length, size_t indent) {
    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < indent; i++) {
        fprintf(stderr, "\t");
    }

    print_vector(stderr, node->data.state_vec, vec_length);
    fprintf(stderr, "\n");

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        //printf("Edge weight %f node data %f ", values[i].weight,
        //       ((struct node*)(&values[i].node))->data.alpha);

        //print_vector(((struct node*)(&values[i].node))->data.state_vec, vec_length);

        for (size_t k = 0; k < indent; k++) {
            fprintf(stderr, "\t");
        }

        fprintf(stderr, "Edge weight %Lf:\n ", values[i].weight);
        print_graph_node((coal_graph_node_t*)values[i].node,vec_length, indent+1);
        fprintf(stderr, "\n");
    }
}

static void _print_graph_list(FILE *stream, coal_graph_node_t *node,
        bool indexed,
        size_t vec_length, size_t vec_spacing) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;

    fprintf(stream, "Node: ");
    print_vector_spacing(stream, node->data.state_vec,
            vec_length, vec_spacing);
    if (indexed) {
        fprintf(stream, " (%zu)", node->data.vertex_index);
    }
    fprintf(stream, ":\n");

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        fprintf(stream, "\t");
        fprintf(stream, "(%Lf) ", values[i].weight);
        print_vector_spacing(stream,((coal_graph_node_t*)values[i].node)->data.state_vec,
                vec_length, vec_spacing);
        if (indexed) {
            fprintf(stream, " (%zu)", ((coal_graph_node_t*)values[i].node)->data.vertex_index);
        }
        if (vector_length(((coal_graph_node_t*)values[i].node)->edges) == 0) {
            fprintf(stream, " (abs)");
        }
        fprintf(stream, "\n");
    }

    fprintf(stream, "\n");
    for (size_t i = 0; i < vector_length(node->edges); i++) {
        _print_graph_list(stream, (coal_graph_node_t *) values[i].node,
                indexed,
                vec_length, vec_spacing);
    }
}

void coal_print_graph_list(FILE *stream, coal_graph_node_t *graph,
                           bool indexed,
                           size_t vec_length, size_t vec_spacing) {
    if (indexed) {
        size_t largest_index;
        coal_label_vertex_index(&largest_index, graph);
    }

    reset_graph_visited(graph);
    _print_graph_list(stream, graph, indexed, vec_length, vec_spacing);
    fflush(stream);
}

/*
 * Assumes visited state reset and indexed vertices.
 */
void insert_into_weight_mat(weight_t **weights, coal_graph_node_t *node) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;


    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

        weights[node->data.vertex_index][child->data.vertex_index] = values[i].weight;
        weights[node->data.vertex_index][node->data.vertex_index] -= values[i].weight;
    }

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

        insert_into_weight_mat(weights, child);
    }
}

int coal_graph_as_mat(weight_t ***weights, size_t *out_size, coal_graph_node_t *graph) {
    size_t largest_index;
    coal_label_vertex_index(&largest_index, graph);
    reset_graph_visited(graph);
    size_t size = largest_index+1;
    *out_size = size;
    *weights = malloc(sizeof(weight_t*)*size);

    for (size_t i = 0; i < size; ++i) {
        (*weights)[i] = calloc(size, sizeof(weight_t));
    }

    insert_into_weight_mat(*weights, graph);
    return 0;
}

void insert_into_gsl_weight_mat(gsl_matrix_long_double *weights,
        coal_graph_node_t *node,
        bool include_absorbing) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;

    weighted_edge_t *values = vector_get(node->edges);
    weight_t weight = 0;
    size_t node_index = (size_t) node->data.vertex_index;

    if (!include_absorbing) {
        if (node_index == 0) {
            return;
        }

        node_index--;
    }

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        weight += values[i].weight;
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

        size_t child_index = (size_t) child->data.vertex_index;

        if (!include_absorbing) {
            if (child_index == 0) {
                continue;
            }

            child_index--;
        }

        gsl_matrix_long_double_set(weights,
                                   (const size_t) node_index,
                                   (const size_t) child_index,
                                   values[i].weight);
    }

    gsl_matrix_long_double_set(weights,
                               (const size_t) node_index,
                               (const size_t) node_index,
                               -weight);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

        insert_into_gsl_weight_mat(weights, child, include_absorbing);
    }
}


int coal_graph_as_gsl_mat(gsl_matrix_long_double **weights, coal_graph_node_t *graph, bool include_absorbing) {
    size_t largest_index;
    coal_label_vertex_index(&largest_index, graph);
    reset_graph_visited(graph);
    size_t size = largest_index+1;

    if (!include_absorbing) {
        size--;
    }

    *weights = gsl_matrix_long_double_calloc(size, size);

    insert_into_gsl_weight_mat(*weights, graph, include_absorbing);

    return 0;
}

void insert_into_gsl_weight_mat_discrete(gsl_matrix_long_double *weights,
                                coal_graph_node_t *node,
                                bool include_absorbing) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;

    weighted_edge_t *values = vector_get(node->edges);
    weight_t sum = 0;
    size_t node_index = (size_t) node->data.vertex_index;

    if (!include_absorbing) {
        if (node_index == 0) {
            return;
        }

        node_index--;
    }

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;
        size_t child_index = (size_t) child->data.vertex_index;
        sum += values[i].weight;

        if (!include_absorbing) {
            if (child_index == 0) {
                continue;
            }

            child_index--;
        }

        gsl_matrix_long_double_set(weights,
                                   (const size_t) node_index,
                                   (const size_t) child_index,
                                   values[i].weight);
    }

    gsl_matrix_long_double_set(weights,
                               (const size_t) node_index,
                               (const size_t) node_index,
                               1 - sum);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

        insert_into_gsl_weight_mat_discrete(weights, child, include_absorbing);
    }
}

int coal_graph_as_gsl_mat_discrete(gsl_matrix_long_double **weights, coal_graph_node_t *graph, bool include_absorbing) {
    size_t largest_index;
    coal_label_vertex_index(&largest_index, graph);
    reset_graph_visited(graph);
    size_t size = largest_index+1;

    if (!include_absorbing) {
        size--;
    }

    *weights = gsl_matrix_long_double_calloc(size, size);

    insert_into_gsl_weight_mat_discrete(*weights, graph, include_absorbing);

    return 0;
}

void insert_into_gsl_weight_mat_double(gsl_matrix *weights,
                                coal_graph_node_t *node,
                                bool include_absorbing) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;

    weighted_edge_t *values = vector_get(node->edges);
    weight_t weight = 0;
    size_t node_index = (size_t) node->data.vertex_index;

    if (!include_absorbing) {
        if (node_index == 0) {
            return;
        }

        node_index--;
    }

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        weight += values[i].weight;
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

        size_t child_index = (size_t) child->data.vertex_index;

        if (!include_absorbing) {
            if (child_index == 0) {
                continue;
            }

            child_index--;
        }

        gsl_matrix_set(weights,
                                   (const size_t) node_index,
                                   (const size_t) child_index,
                       (const double) values[i].weight);
    }

    gsl_matrix_set(weights,
                               (const size_t) node_index,
                               (const size_t) node_index,
                   (const double) -weight);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

        insert_into_gsl_weight_mat_double(weights, child, include_absorbing);
    }
}

int coal_graph_as_gsl_mat_double(gsl_matrix **weights, coal_graph_node_t *graph, bool include_absorbing) {
    size_t largest_index;
    coal_label_vertex_index(&largest_index, graph);
    reset_graph_visited(graph);
    size_t size = largest_index+1;

    if (!include_absorbing) {
        size--;
    }

    *weights = gsl_matrix_calloc(size, size);

    insert_into_gsl_weight_mat_double(*weights, graph, include_absorbing);

    return 0;
}

int coal_gen_kingman_graph(coal_graph_node_t **graph, size_t n) {
    size_t state_size = n;

    vec_entry_t *initial = (vec_entry_t*)calloc(state_size, sizeof(vec_entry_t));
    initial[0] = state_size;

    vec_entry_t *mrca = (vec_entry_t*)calloc(state_size, sizeof(vec_entry_t));
    mrca[state_size-1] = 1;

    coal_graph_node_t *absorbing_vertex;
    coal_graph_node_create(&absorbing_vertex, mrca, mrca);

    avl_vec_node_t *BST;
    avl_vec_node_create(&BST, mrca, absorbing_vertex, NULL);

    coal_graph_node_t *state_graph;

    kingman_visit_vertex(&state_graph, initial, BST, state_size, state_size,
                         state_size);

    vec_entry_t *start_state = (vec_entry_t*)calloc(state_size, sizeof(vec_entry_t));

    coal_graph_node_t *start;
    coal_graph_node_create(&start, start_state, start_state);
    graph_add_edge((graph_node_t*) start,
                   (graph_node_t*) state_graph, 1);

    *graph = start;

    return 0;
}

int coal_gen_kingman_graph_rw(coal_graph_node_t **graph, size_t n, size_t rw_index) {
    size_t state_size = n+1;

    vec_entry_t *initial = (vec_entry_t*)calloc(state_size, sizeof(vec_entry_t));
    initial[0] = n;

    vec_entry_t *mrca = (vec_entry_t*)calloc(state_size, sizeof(vec_entry_t));
    mrca[state_size-1] = 1;

    coal_graph_node_t *absorbing_vertex;
    coal_graph_node_create(&absorbing_vertex, mrca, mrca);

    avl_vec_node_t *BST;
    avl_vec_node_create(&BST, mrca, absorbing_vertex, NULL);

    coal_graph_node_t *state_graph;

    kingman_visit_vertex_rw(&state_graph, rw_index, initial, BST, state_size, state_size,
                         state_size);

    vec_entry_t *start_state = (vec_entry_t*)calloc(state_size, sizeof(vec_entry_t));

    coal_graph_node_t *start;
    coal_graph_node_create(&start, start_state, start_state);
    graph_add_edge((graph_node_t*) start,
                   (graph_node_t*) state_graph, 1);

    *graph = start;

    return 0;
}

static int im_state_mat_init(vec_entry_t ***mat, size_t n1, size_t n2) {
    *mat = malloc(sizeof(vec_entry_t *) * (n1 + 1));

    for (size_t i = 0; i < n1 + 1; ++i) {
        (*mat)[i] = calloc(n2 + 1, sizeof(vec_entry_t));
    }

    return 0;
}

/*
 * Isolation with migration
 */
int im_state_init(im_state_t **out, size_t n1, size_t n2) {
    *out = malloc(sizeof(im_state_t));
    im_state_mat_init(&((*out)->mat1), n1, n2);
    im_state_mat_init(&((*out)->mat2), n1, n2);
    (*out)->n1 = n1;
    (*out)->n2 = n2;
    (*out)->in_iso = true;
    (*out)->flag_mig1to2 = true;
    (*out)->flag_mig2to1 = true;

    return 0;
}

int im_state_cpy(im_state_t **out,
                  const im_state_t *in,
                  size_t n1, size_t n2) {
    im_state_init(out, n1, n2);

    for (size_t i = 0; i < n1 + 1; ++i) {
        for (size_t j = 0; j < n2 + 1; ++j) {
            (*out)->mat1[i][j] = in->mat1[i][j];
            (*out)->mat2[i][j] = in->mat2[i][j];
        }
    }

    (*out)->in_iso = in->in_iso;
    (*out)->flag_mig1to2 = in->flag_mig1to2;
    (*out)->flag_mig2to1 = in->flag_mig2to1;

    return 0;
}

size_t im_state_length(size_t n1, size_t n2) {
    return (n1+1) * (n2+1) * 2 + 3;
}

void im_state_mat_as_vec(vec_entry_t *out, vec_entry_t **mat,
        size_t width, size_t height) {
    size_t k = 0;

    for (size_t j = 0; j < height; ++j) {
        for (size_t i = 0; i < width; ++i) {
            out[k] = mat[i][j];
            k++;
        }
    }
}

int im_state_as_vec(vec_entry_t **out, im_state_t *state,
                    size_t n1, size_t n2) {
    *out = malloc(sizeof(vec_entry_t)* im_state_length(n1, n2));

    size_t length = (n1+1)*(n2+1);
    im_state_mat_as_vec(*out, state->mat1, n1+1, n2+1);
    im_state_mat_as_vec(*out+length, state->mat2, n1+1, n2+1);
    *(*out+length*2) = (vec_entry_t)state->in_iso;
    *(*out+length*2+1) = (vec_entry_t)state->flag_mig1to2;
    *(*out+length*2+2) = (vec_entry_t)state->flag_mig2to1;

    return 0;
}

static int
im_visit_vertex(coal_graph_node_t **out, im_state_t *state,
                vec_entry_t *state_vec, avl_vec_node_t *bst,
                size_t num_coal_events,
                size_t vector_length,
                coal_gen_im_graph_args_t *args);

static size_t count_lineages_mat(vec_entry_t **mat, size_t n1, size_t n2);
static size_t count_lineages(im_state_t *state, size_t n1, size_t n2);

/*
 * Adds mat2 into mat
 */
static void combine_im_matrix(vec_entry_t ***mat_out,
        const vec_entry_t **mat1, const vec_entry_t **mat2,
        size_t n1, size_t n2) {
    im_state_mat_init(mat_out, n1, n2);

    for (vec_entry_t i = 0; i < n1+1; i++) {
        for (vec_entry_t j = 0; j < n2+1; j++) {
            (*mat_out)[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
}

static void empty_im_matrix(vec_entry_t ***mat_out,
                              size_t n1, size_t n2) {
    im_state_mat_init(mat_out, n1, n2);
}

static inline int im_visit_coal_loop(coal_graph_node_t **out,
        vec_entry_t **d,
        coal_param_real_t scale,
        im_state_t *state, avl_vec_node_t *bst,
        size_t num_coal_events,
        size_t vector_length,
        coal_gen_im_graph_args_t *args) {
    if (scale == 0) {
        return 0;
    }

    for (size_t i1 = 0; i1 < args->n1+1; i1++) {
        for (size_t j1 = 0; j1 < args->n2+1; j1++) {
            for (size_t i2 = 0; i2 < args->n1+1; i2++) {
                for (size_t j2 = 0; j2 < args->n2+1; j2++) {
                    if (i2 < i1 || (i2 == i1 && j2 < j1)) {
                        continue;
                    }

                    mat_entry_t t;

                    if (i1 == i2 && j1 == j2) {
                        if (d[i1][j1] >= 2) {
                            t = d[i1][j1] * (d[i1][j1] - 1) / 2.0f;
                        } else {
                            continue;
                        }
                    } else {
                        if (d[i1][j1] >= 1 && d[i2][j2] >= 1) {
                            t = d[i1][j1] * d[i2][j2];
                        } else {
                            continue;
                        }
                    }

                    t *= scale;

                    bool old_flag_mig1to2 = state->flag_mig1to2;
                    bool old_flag_mig2to1 = state->flag_mig2to1;
                    d[i1][j1]--;
                    d[i2][j2]--;
                    d[i1 + i2][j1 + j2]++;

                    if (state->in_iso) {
                        state->flag_mig1to2 = true;
                        state->flag_mig2to1 = true;
                    } else {
                        state->flag_mig1to2 = false;
                        state->flag_mig2to1 = false;
                    }

                    vec_entry_t *v;
                    coal_graph_node_t *new_vertex;

                    if (num_coal_events + 1 == args->num_iso_coal_events) {
                        vec_entry_t **combined_mat;
                        vec_entry_t **empty_mat;
                        vec_entry_t **old_mat1;
                        vec_entry_t **old_mat2;

                        old_mat1 = state->mat1;
                        old_mat2 = state->mat2;

                        combine_im_matrix(&combined_mat, (const vec_entry_t **) state->mat1,
                                          (const vec_entry_t **) state->mat2, args->n1, args->n2);
                        empty_im_matrix(&empty_mat, args->n1, args->n2);

                        state->flag_mig1to2 = false;
                        state->flag_mig2to1 = false;
                        state->mat1 = combined_mat;
                        state->mat2 = empty_mat;
                        state->in_iso = false;

                        im_state_as_vec(&v, state, args->n1, args->n2);
                        im_visit_vertex(&new_vertex,
                                        state, v,
                                        bst,
                                        num_coal_events + 1,
                                        vector_length,
                                        args);

                        state->in_iso = true;
                        state->mat2 = old_mat2;
                        state->mat1 = old_mat1;

                        // It is possible that two different events lead
                        // to the same state after the iso period. Therefore,
                        // we should possibly combine them into one edge.
                        graph_combine_edge((graph_node_t *) *out,
                                       (graph_node_t *) new_vertex, t);
                    } else {
                        im_state_as_vec(&v, state, args->n1, args->n2);
                        im_visit_vertex(&new_vertex,
                                        state, v,
                                        bst,
                                        num_coal_events + 1,
                                        vector_length,
                                        args);

                        graph_add_edge((graph_node_t *) *out,
                                       (graph_node_t *) new_vertex, t);
                    }

                    new_vertex->data.type = IM_VERTEX_TYPE_COAL;
                    state->flag_mig2to1 = old_flag_mig2to1;
                    state->flag_mig1to2 = old_flag_mig1to2;
                    d[i1 + i2][j1 + j2]--;
                    d[i2][j2]++;
                    d[i1][j1]++;
                }
            }
        }
    }

    return 0;
}

static inline int im_visit_mig_loop(coal_graph_node_t **out,
                                     vec_entry_t **d_from,
                                     vec_entry_t **d_to,
                                     bool *flag_mig_from,
                                     bool *flag_mig_to,
                                     coal_param_real_t scale,
                                     im_state_t *state, avl_vec_node_t *bst,
                                     size_t num_coal_events,
                                     size_t vector_length,
                                     coal_gen_im_graph_args_t *args) {
    if (scale == 0) {
        return 0;
    }

    for (size_t i = 0; i < args->n1+1; i++) {
        for (size_t j = 0; j < args->n2+1; j++) {
            // The MRCA cannot migrate
            if (i == args->n1 && j == args->n2) {
                continue;
            }

            // If flag is not set, we cannot migrate
            if (!*flag_mig_from) {
                continue;
            }

            mat_entry_t t;

            if (d_from[i][j] >= 1) {
                t = d_from[i][j];
            } else {
                continue;
            }

            t *= scale;

            bool old_flag_mig_to = *flag_mig_to;
            bool old_flag_mig_from = *flag_mig_from;
            d_from[i][j]--;
            d_to[i][j]++;

            switch (args->migration_type) {
                case MIG_ALL:
                    *flag_mig_to = true;
                    *flag_mig_from = true;
                    break;
                case MIG_DIR:
                    *flag_mig_to = false;
                    *flag_mig_from = true;
                    break;
                case MIG_ONCE:
                    *flag_mig_to = false;
                    *flag_mig_from = false;
                    break;
                default:
                    DIE_ERROR(1, "Illegal mig param, was %zu\n", (size_t)args->migration_type);
                    break;
            }

            vec_entry_t *v;
            im_state_as_vec(&v, state, args->n1, args->n2);

            coal_graph_node_t *new_vertex;
            im_visit_vertex(&new_vertex,
                            state, v,
                            bst,
                            num_coal_events,
                            vector_length,
                            args);

            new_vertex->data.type = IM_VERTEX_TYPE_MIG;
            *flag_mig_from = old_flag_mig_from;
            *flag_mig_to = old_flag_mig_to;
            d_to[i][j]--;
            d_from[i][j]++;

            graph_add_edge((graph_node_t *) *out,
                           (graph_node_t *) new_vertex, t);
        }
    }

    return 0;
}

static void reset_im_flag_matrix(vec_entry_t **flag_mat, size_t n1, size_t n2) {
    for (vec_entry_t i = 0; i < n1+1; i++) {
        for (vec_entry_t j = 0; j < n2 + 1; j++) {
            flag_mat[i][j] = 0;
        }
    }
}

static int
im_visit_vertex(coal_graph_node_t **out, im_state_t *state,
                vec_entry_t *state_vec, avl_vec_node_t *bst,
                size_t num_coal_events,
                size_t vector_length,
                coal_gen_im_graph_args_t *args) {
    avl_vec_node_t *bst_node = avl_vec_find(bst, state_vec, vector_length);

    if (bst_node != NULL) {
        *out = bst_node->entry;
        return 0;
    } else {
        //print_vector_spacing(stderr, state_vec, im_state_length(args->n1, args->n2), (args->n1+1));
        //fprintf(stderr, "\n");
        //fflush(stderr);
        // Never pass this on
        im_state_t *vertex_state;
        im_state_cpy(&vertex_state, state, args->n1, args->n2);

        vec_entry_t *vertex_state_vec;
        im_state_as_vec(&vertex_state_vec, vertex_state, args->n1, args->n2);

        coal_graph_node_create(out, state_vec, vertex_state);
        (*out)->data.coals = num_coal_events;

        avl_vec_insert(&bst, state_vec, *out, vector_length);

        if (num_coal_events < args->num_iso_coal_events) {
            // Special case. We do not want to run into dead ends
            bool may_mig;

            if (args->migration_type != MIG_ONCE) {
                may_mig = true;
            } else {
                size_t n1 = count_lineages_mat(state->mat1, args->n1, args->n2);
                size_t n2 = count_lineages_mat(state->mat2, args->n1, args->n2);

                if ((n1 == 2 && n2 == 0) || (n2 == 2 && n1 == 0)) {
                    may_mig = false;
                } else {
                    may_mig = true;
                }
            }

            if (may_mig) {
                im_visit_mig_loop(out, state->mat1, state->mat2,
                                  &(state->flag_mig1to2),
                                  &(state->flag_mig2to1),
                                  args->pop_scale1 * args->mig_scale1,
                                  state, bst, num_coal_events, vector_length, args);
                im_visit_mig_loop(out, state->mat2, state->mat1,
                                  &(state->flag_mig2to1),
                                  &(state->flag_mig1to2),
                                  args->pop_scale2 * args->mig_scale2,
                                  state, bst, num_coal_events, vector_length, args);
            }

            im_visit_coal_loop(out, state->mat1, args->pop_scale1, state,
                               bst, num_coal_events, vector_length, args);
            im_visit_coal_loop(out, state->mat2, args->pop_scale2, state,
                               bst, num_coal_events, vector_length, args);
        } else {
            coal_param_real_t scale = 1;
            im_visit_coal_loop(out, state->mat1, scale, state,
                               bst, num_coal_events, vector_length, args);
        }

        return 0;
    }
}


static int
im_ss_visit_vertex(coal_graph_node_t **out,
                avl_vec_node_t *bst,
                vec_entry_t npop1, vec_entry_t npop2,
                bool mig1to2, bool mig2to1,
                coal_gen_im_graph_args_t *args) {
    const size_t vector_length = 4;
    vec_entry_t *state = calloc(vector_length, sizeof(vec_entry_t));
    if (npop1 + npop2 == 1) {
        // If only one left, use left population as absorbing state
        state[0] = 1;
        state[1] = 0;
    } else {
        state[0] = npop1;
        state[1] = npop2;
    }
    state[2] = (vec_entry_t) mig1to2;
    state[3] = (vec_entry_t) mig2to1;

    size_t num_coal_events = (args->n1 + args->n2) - (npop1 + npop2);

    if (num_coal_events >= args->num_iso_coal_events) {
        npop1 += npop2;
        npop2 = 0;

        state[2] = (vec_entry_t) true;
        state[3] = (vec_entry_t) true;
    }

    avl_vec_node_t *bst_node = avl_vec_find(bst, state, vector_length);

    if (bst_node != NULL) {
        *out = bst_node->entry;
        return 0;
    } else {
        coal_graph_node_create(out, state, state);

        (*out)->data.coals = num_coal_events;

        avl_vec_insert(&bst, state, *out, vector_length);

        bool new_mig1to2;
        bool new_mig2to1;

        // Coal events
        if (npop1 > 1) {
            coal_param_real_t scale = args->pop_scale1;

            if (scale >= LDBL_EPSILON) {
                weight_t rate = scale * npop1 * (npop1 - 1) / 2.0f;

                new_mig1to2 = true;
                new_mig2to1 = true;

                coal_graph_node_t *to;
                im_ss_visit_vertex(&to, bst,
                                   npop1 - 1, npop2,
                                   new_mig1to2, new_mig2to1, args);
                graph_add_edge((graph_node_t *) *out, (graph_node_t *) to, rate);
            }
        }

        if (npop2 > 1) {
            coal_param_real_t scale = args->pop_scale2;

            if (scale >= LDBL_EPSILON) {
                weight_t rate = scale * npop2 * (npop2 - 1) / 2.0f;

                new_mig1to2 = true;
                new_mig2to1 = true;

                coal_graph_node_t *to;
                im_ss_visit_vertex(&to, bst,
                                   npop1, npop2 - 1,
                                   new_mig1to2, new_mig2to1, args);
                graph_add_edge((graph_node_t *) *out, (graph_node_t *) to, rate);
            }
        }

        if (num_coal_events < args->num_iso_coal_events) {
            if (npop1 + npop2 > 1) {
                // Mig events
                if (mig1to2 && npop1 > 0 &&
                    (args->migration_type != MIG_ONCE ||
                     !(npop1 == 2 && npop2 == 0))) {
                    coal_param_real_t scale = args->pop_scale1 * args->mig_scale1;

                    if (scale >= LDBL_EPSILON) {
                        weight_t rate = scale * npop1;

                        switch (args->migration_type) {
                            case MIG_DIR:
                                new_mig1to2 = true;
                                new_mig2to1 = false;
                                break;
                            case MIG_ALL:
                                new_mig1to2 = true;
                                new_mig2to1 = true;
                                break;
                            case MIG_ONCE:
                                new_mig1to2 = false;
                                new_mig2to1 = false;
                                break;
                            default:
                                DIE_ERROR(1, "Illegal mig param, was %zu\n", (size_t) args->migration_type);
                                break;
                        }

                        coal_graph_node_t *to;
                        im_ss_visit_vertex(&to, bst,
                                           npop1 - 1, npop2 + 1,
                                           new_mig1to2, new_mig2to1, args);
                        graph_add_edge((graph_node_t *) *out, (graph_node_t *) to, rate);
                    }
                }

                if (mig2to1 && npop2 > 0 &&
                    (args->migration_type != MIG_ONCE ||
                     !(npop2 == 2 && npop1 == 0))) {
                    coal_param_real_t scale = args->pop_scale2 * args->mig_scale2;

                    if (scale >= LDBL_EPSILON) {
                        weight_t rate = scale * npop2;

                        switch (args->migration_type) {
                            case MIG_DIR:
                                new_mig2to1 = true;
                                new_mig1to2 = false;
                                break;
                            case MIG_ALL:
                                new_mig2to1 = true;
                                new_mig1to2 = true;
                                break;
                            case MIG_ONCE:
                                new_mig2to1 = false;
                                new_mig1to2 = false;
                                break;
                            default:
                                DIE_ERROR(1, "Illegal mig param, was %zu\n", (size_t) args->migration_type);
                                break;
                        }

                        coal_graph_node_t *to;
                        im_ss_visit_vertex(&to, bst,
                                           npop1 + 1, npop2 - 1,
                                           new_mig1to2, new_mig2to1, args);
                        graph_add_edge((graph_node_t *) *out, (graph_node_t *) to, rate);
                    }
                }
            }
        }

        return 0;
    }
}

int coal_gen_im_graph(coal_graph_node_t **graph, avl_vec_node_t **bst, coal_gen_im_graph_args_t args) {
    size_t n1 = args.n1;
    size_t n2 = args.n2;

    if (args.num_iso_coal_events >= n1 + n2) {
        DIE_ERROR(1, "There can only be n1+n2-1 coalescence events. Got n1: %zu, n2: %zu, events: %zu\n",
                n1, n2, args.num_iso_coal_events);
    }

    im_state_t *initial;
    im_state_init(&initial, n1, n2);

    if (n1 > 0) {
        initial->mat1[1][0] = n1;
    }

    if (n2 > 0) {
        initial->mat2[0][1] = n2;
    }

    if (args.num_iso_coal_events == 0) {
        vec_entry_t **combined;
        vec_entry_t **empty;
        combine_im_matrix(&combined,
                          (const vec_entry_t **) initial->mat1,
                          (const vec_entry_t **) initial->mat2, n1, n2);
        empty_im_matrix(&empty, n1, n2);

        initial->mat1 = combined;
        initial->mat2 = empty;
        initial->in_iso = false;
        initial->flag_mig1to2 = false;
        initial->flag_mig2to1 = false;
    }

    vec_entry_t *initial_vec;
    im_state_as_vec(&initial_vec, initial, n1, n2);

    im_state_t *mrca;
    im_state_init(&mrca, n1, n2);
    mrca->mat1[n1][n2] = 1;
    mrca->in_iso = false;

    vec_entry_t *mrca_vectype;

    im_state_as_vec(&mrca_vectype, mrca, n1, n2);

    im_state_t *mrcatype1;
    im_state_init(&mrcatype1, n1, n2);
    im_state_t *mrcatype2;
    im_state_init(&mrcatype2, n1, n2);
    mrcatype1->mat1[n1][n2] = 1;
    mrcatype2->mat2[n1][n2] = 1;

    vec_entry_t *mrca_vectype1;
    vec_entry_t *mrca_vectype2;

    im_state_as_vec(&mrca_vectype1, mrcatype1, n1, n2);
    im_state_as_vec(&mrca_vectype2, mrcatype2, n1, n2);

    coal_graph_node_t *absorbing_vertex;
    coal_graph_node_create(&absorbing_vertex, mrca_vectype, mrca);


    // Either of the two matrices when n1,n2=1, is MRCA.
    avl_vec_node_t *BST;
    avl_vec_node_create(&BST, mrca_vectype, absorbing_vertex, NULL);
    avl_vec_insert(&BST, mrca_vectype1, absorbing_vertex, im_state_length(n1,n2));
    avl_vec_insert(&BST, mrca_vectype2, absorbing_vertex, im_state_length(n1,n2));

    coal_graph_node_t *state_graph;

    im_visit_vertex(&state_graph, initial, initial_vec, BST,
                    0, im_state_length(n1, n2), &args);

    state_graph->data.type = IM_VERTEX_TYPE_START;
    state_graph->data.coals = 0;

    *graph = state_graph;

    if (bst != NULL) {
        *bst = BST;
    }

    return 0;
}

int coal_gen_im_ss_graph(coal_graph_node_t **graph, coal_gen_im_graph_args_t args) {
    // Add an entry to the BST so we don't have a null pointer
    vec_entry_t *unused_vec = calloc(4, sizeof(vec_entry_t));
    unused_vec[0]=0;
    unused_vec[1]=0;
    unused_vec[2]=0;
    unused_vec[3]=0;
    avl_vec_node_t *BST;
    avl_vec_node_create(&BST, unused_vec, NULL, NULL);

    coal_graph_node_t *state_graph;

    im_ss_visit_vertex(&state_graph, BST, args.n1, args.n2,
            true, true, &args);

    *graph = state_graph;

    return 0;
}


static size_t count_lineages_mat(vec_entry_t **mat, size_t n1, size_t n2) {
    size_t sum = 0;

    for (size_t i = 0; i < n1+1; ++i) {
        for (size_t j = 0; j < n2+1; ++j) {
            sum += mat[i][j];
        }
    }

    return sum;
}

static size_t count_lineages(im_state_t *state, size_t n1, size_t n2) {
    return count_lineages_mat(state->mat1, n1, n2) +
            count_lineages_mat(state->mat2, n1, n2);
}

coal_graph_node_t *debug_graph;

int cutoff_remove_wrong_left(coal_graph_node_t *graph,
        coal_graph_node_t *absorbing_vertex,
        size_t n1, size_t n2,
        size_t left1, size_t left2) {
    if (graph->data.visited) {
        return 0;
    }

    graph->data.visited = true;

    size_t parent_lineages = graph->data.state_vec[0]+graph->data.state_vec[1];
    weighted_edge_t *values = vector_get(graph->edges);

    if (vector_length(graph->edges) > 0) {
        for (ssize_t i = vector_length(graph->edges) - 1; i >= 0; i--) {
            coal_graph_node_t *child = ((coal_graph_node_t *) values[i].node);
            size_t child_lineages = child->data.state_vec[0] + child->data.state_vec[1];

            bool has_coalesced;

            if (child_lineages < parent_lineages) {
                has_coalesced = true;
            } else {
                has_coalesced = false;
            }

            print_vector_spacing(stderr, graph->data.state_vec, 4, 4);
            fprintf(stderr, "->");
            print_vector_spacing(stderr, child->data.state_vec, 4, 4);
            fprintf(stderr, " coal: %d, clin: %zu", has_coalesced, child_lineages);
            fprintf(stderr, "\n");

            if (has_coalesced && child_lineages == left1 + left2) {
                // It has coalesced, and we have the right or a smaller total
                // number of lineages left

                if (child == absorbing_vertex) {
                    // This edge is fine
                    continue;
                }

                //if (child->data.state_vec[0] != left1 || child->data.state_vec[1] != left2) {
                if (false) {
                    fprintf(stderr, "Killing edge!\n");
                    graph_redistribute_edge((graph_node_t *) graph, (graph_node_t *) child);
                } else {
                    fprintf(stderr, "THIS SHOULD BE THE ABSORBING VERTEX:");
                    print_vector_spacing(stderr, child->data.state_vec, 4, 4);
                    fprintf(stderr, "\n");
                    fprintf(stderr, "Making edge to abs\n");
                    vector_clear(child->edges);

                    for (size_t j = 0; j < vector_length(graph->edges); j++) {
                        coal_graph_node_t *child2 = ((coal_graph_node_t *) values[j].node);

                        if (child2 == absorbing_vertex) {
                            DIE_ERROR(1, "Absorbing vertex already exists\n");
                        }
                    }

                    // Make edge go to the absorbing vertex instead of child
                    //values[i].node = (graph_node_t*)absorbing_vertex;
                    // TODO: for now we just add a very weighted edge
                    graph_add_edge((graph_node_t *) child, (graph_node_t *) absorbing_vertex, 99999);
                }
            }
        }
    }

    for (size_t i = 0; i < vector_length(graph->edges); i++) {
        coal_graph_node_t *child = ((coal_graph_node_t*)values[i].node);
        size_t child_lineages = child->data.state_vec[0] + child->data.state_vec[1];

        if (child_lineages > left1 + left2) {
            cutoff_remove_wrong_left(child, absorbing_vertex, n1, n2, left1, left2);
        }
    }
}

int remove_nonedge_state(coal_graph_node_t **removed,
        coal_graph_node_t *graph,
        coal_graph_node_t *absorbing_vertex) {
    if (graph->data.visited) {
        return 0;
    }

    graph->data.visited = true;

    weighted_edge_t *values = vector_get(graph->edges);

    for (size_t i = 0; i < vector_length(graph->edges); i++) {
        coal_graph_node_t *child = ((coal_graph_node_t*)values[i].node);
        size_t child_edges = vector_length(child->edges);


        print_vector_spacing(stderr, graph->data.state_vec, 4, 4);
        fprintf(stderr, "->");
        print_vector_spacing(stderr, child->data.state_vec, 4, 4);
        fprintf(stderr, " edges: %zu", child_edges);
        fprintf(stderr, "\n");

        if (child_edges == 0 && child != absorbing_vertex) {
            *removed = child;
            graph_redistribute_edge((graph_node_t *)graph, (graph_node_t *)child);
            fprintf(stderr, "NE: Killing edge");
            //print_vector_spacing(stderr, graph->data.state_vec, 4, 4);
            //fprintf(stderr, " to ");
            //print_vector_spacing(stderr, child->data.state_vec, 4, 4);
            //fprintf(stderr, "\n");
            return 0;
        }

        remove_nonedge_state(removed, child, absorbing_vertex);
    }

    return 0;
}

int cutoff_remove_after_left(coal_graph_node_t *graph,
                             coal_graph_node_t *absorbing_vertex,
                             size_t left) {
    if (graph->data.visited) {
        return 0;
    }

    graph->data.visited = true;

    size_t parent_lineages = graph->data.state_vec[0]+graph->data.state_vec[1];
    weighted_edge_t *values = vector_get(graph->edges);

    for (size_t i = 0; i < vector_length(graph->edges); i++) {
        coal_graph_node_t *child = ((coal_graph_node_t*)values[i].node);
        size_t child_lineages = child->data.state_vec[0]+child->data.state_vec[1];

        bool has_coalesced;

        if (child_lineages < parent_lineages) {
            has_coalesced = true;
        } else {
            has_coalesced = false;
        }

        if (has_coalesced && child_lineages <= left) {
            values[i].node = (graph_node_t*)absorbing_vertex;
        }

        cutoff_remove_after_left(child, absorbing_vertex, left);
    }
}

int cutoff_redirect_wrong_left(coal_graph_node_t **correct_vertex,
                               coal_graph_node_t *graph,
                               coal_graph_node_t *absorbing_vertex,
                               size_t n1, size_t n2,
                               size_t left1, size_t left2) {
    if (graph->data.visited) {
        return 0;
    }

    graph->data.visited = true;

    size_t parent_lineages = graph->data.state_vec[0]+graph->data.state_vec[1];
    weighted_edge_t *values = vector_get(graph->edges);

    for (size_t i = 0; i < vector_length(graph->edges); i++) {
        coal_graph_node_t *child = ((coal_graph_node_t*)values[i].node);
        size_t child_lineages = child->data.state_vec[0]+child->data.state_vec[1];

        bool has_coalesced;

        if (child_lineages < parent_lineages) {
            has_coalesced = true;
        } else {
            has_coalesced = false;
        }

        if (has_coalesced && child_lineages <= left1 + left2) {
            // It has coalesced, and we have the right or a smaller total
            // number of lineages left

            if (child->data.state_vec[0] != left1 ||
                child->data.state_vec[1] != left2) {
                // We do not kill the edge in this case, but make the child
                // vertex into the absorbing state
                values[i].node = (graph_node_t*)absorbing_vertex;
            } else {
                //fprintf(stderr, "We have the correct vertex\n");
                //fprintf(stderr, "Before:\n");
                //coal_print_graph_list(stderr, debug_graph, false, 4, 4);
                // This is the correct vertex
                *correct_vertex = child;
                // Remove all the edges from the child
                vector_clear(child->edges);
                //fprintf(stderr, "Post clear:\n");
                //coal_print_graph_list(stderr, debug_graph, false, 4, 4);
                if (absorbing_vertex != child) {
                    // Add edge to absorbing vertex
                    graph_add_edge((graph_node_t *) child,
                                   (graph_node_t *) absorbing_vertex, 1);
                }
                //fprintf(stderr, "Post abs edge:\n");
                //coal_print_graph_list(stderr, debug_graph, false, 4, 4);
            }
        }

        cutoff_redirect_wrong_left(correct_vertex, child, absorbing_vertex, n1, n2, left1, left2);
    }
}

int coal_gen_im_pure_cutoff_graph(coal_graph_node_t **graph, coal_gen_im_pure_cutoff_graph_args_t args) {
    size_t n1 = args.n1;
    size_t n2 = args.n2;
    size_t left = args.left;

    if (left > n1 + n2) {
        DIE_ERROR(1, "There has to be at most n1+n2 left. Got n1: %zu, n2: %zu, left: %zu\n",
                   n1, n2, left);
    }

    if (left == 0) {
        DIE_ERROR(1, "There has to be at least one lineage left. Got n1: %zu, n2: %zu, left: %zu\n",
                   n1, n2, left);
    }

    coal_gen_im_graph_args_t im_args = (coal_gen_im_graph_args_t) {
            .n1 = args.n1,
            .n2 = args.n2,
            .num_iso_coal_events = n1 + n2 - 1,
            .migration_type = args.migration_type,
            .pop_scale1 = args.pop_scale1,
            .pop_scale2 = args.pop_scale2,
            .mig_scale1 = args.mig_scale1,
            .mig_scale2 = args.mig_scale2
    };


    coal_graph_node_t *im_graph;
    coal_gen_im_ss_graph(&im_graph, im_args);
    coal_graph_node_t *absorbing_vertex = get_abs_vertex(im_graph);

    coal_graph_reset(im_graph);
    cutoff_remove_after_left(im_graph, absorbing_vertex, left);

    *graph = im_graph;

    return 0;
}


int coal_gen_im_cutoff_graph(coal_graph_node_t **graph, coal_gen_im_cutoff_graph_args_t args) {
    size_t n1 = args.n1;
    size_t n2 = args.n2;
    size_t left1 = args.left_n1;
    size_t left2 = args.left_n2;

    if (left1 + left2 > n1 + n2) {
        DIE_ERROR(1, "There has to be at most n1+n2 left. Got n1: %zu, n2: %zu, left1: %zu, left2: %zu\n",
                   n1, n2, left1, left2);
    }

    if (left1 + left2 == 0) {
        DIE_ERROR(1, "There has to be at least one lineage left. Got n1: %zu, n2: %zu, left1: %zu, left2: %zu\n",
                   n1, n2, left1, left2);
    }

    coal_gen_im_graph_args_t im_args = (coal_gen_im_graph_args_t) {
        .n1 = args.n1,
        .n2 = args.n2,
        .num_iso_coal_events = n1 + n2 - 1,
        .migration_type = args.migration_type,
        .pop_scale1 = args.pop_scale1,
        .pop_scale2 = args.pop_scale2,
        .mig_scale1 = args.mig_scale1,
        .mig_scale2 = args.mig_scale2
    };


    coal_graph_node_t *im_graph;
    coal_gen_im_ss_graph(&im_graph, im_args);
    coal_graph_node_t *absorbing_vertex = get_abs_vertex(im_graph);

    fprintf(stderr, "Args. n1: %zu, n2: %zu, left1: %zu, left2: %zu\n",
            args.n1, args.n2, args.left_n1, args.left_n2);
    fprintf(stderr, "First:\n");
    coal_print_graph_list(stderr, im_graph, false, 4, 4);

    coal_graph_reset(im_graph);
    cutoff_remove_wrong_left(im_graph, absorbing_vertex, n1, n2, left1, left2);
    fprintf(stderr, "Wrong left:\n");
    coal_print_graph_list(stderr, im_graph, false, 4, 4);

    coal_graph_node_t *removed;
    do {
        reset_graph_visited(im_graph);
        removed = NULL;
        remove_nonedge_state(&removed, im_graph, absorbing_vertex);
    } while (removed != NULL);

    fprintf(stderr, "Nowhere:\n");
    //coal_print_graph_list(stderr, im_graph, false, 4, 4);

    *graph = im_graph;

    return 0;
}

int coal_gen_im_prob_vertex_graph(coal_graph_node_t **graph,
        coal_graph_node_t **correct_vertex,
        coal_gen_im_cutoff_graph_args_t args) {
    size_t n1 = args.n1;
    size_t n2 = args.n2;
    size_t left1 = args.left_n1;
    size_t left2 = args.left_n2;

    if (left1 + left2 > n1 + n2) {
        DIE_ERROR(1, "There has to be at most n1+n2 left. Got n1: %zu, n2: %zu, left1: %zu, left2: %zu\n",
                   n1, n2, left1, left2);
    }

    if (left1 + left2 == 0) {
        DIE_ERROR(1, "There has to be at least one lineage left. Got n1: %zu, n2: %zu, left1: %zu, left2: %zu\n",
                   n1, n2, left1, left2);
    }

    coal_gen_im_graph_args_t im_args = (coal_gen_im_graph_args_t) {
            .n1 = args.n1,
            .n2 = args.n2,
            .num_iso_coal_events = n1 + n2 - 1,
            .migration_type = args.migration_type,
            .pop_scale1 = args.pop_scale1,
            .pop_scale2 = args.pop_scale2,
            .mig_scale1 = args.mig_scale1,
            .mig_scale2 = args.mig_scale2
    };


    coal_graph_node_t *im_graph;
    coal_gen_im_ss_graph(&im_graph, im_args);
    coal_graph_node_t *absorbing_vertex = get_abs_vertex(im_graph);

    coal_graph_reset(im_graph);
    debug_graph = im_graph;
    *correct_vertex = NULL;
    cutoff_redirect_wrong_left(correct_vertex, im_graph, absorbing_vertex, n1, n2, left1, left2);

    *graph = im_graph;

    return 0;
}

void _reset_graph_visited(coal_graph_node_t *node, size_t reset_int) {
    if (node->data.reset_int == reset_int) {
        return;
    }

    node->data.reset_int = reset_int;

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        _reset_graph_visited((coal_graph_node_t *) values[i].node, reset_int);
    }

    node->data.visited = false;
    node->data.visits = 0;
}

static void reset_graph_visited(coal_graph_node_t *node) {
    _reset_graph_visited(node, node->data.reset_int + 1);
}

void coal_graph_reset_visited(coal_graph_node_t *graph) {
    reset_graph_visited(graph);
}

void _reset_graph(coal_graph_node_t *node, size_t reset_int) {
    // TODO: This visits multiple times
    if (node->data.reset_int == reset_int) {
        return;
    }

    node->data.reset_int = reset_int;

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        _reset_graph((coal_graph_node_t *) values[i].node, reset_int);
    }

    node->data.full_path_value = -1;
    node->data.vertex_exp = -1;
    node->data.prob = -1;
    node->data.descendants_exp_sum = -1;
    node->data.visited = false;
    node->data.visits = 0;
    node->data.pointer = NULL;
    node->data.all_descendants_exp_sum = NULL;
    node->data.all_vertex_exp = NULL;
}

void reset_graph(coal_graph_node_t *node) {
    _reset_graph(node, node->data.reset_int + 1);
}

void coal_graph_reset(coal_graph_node_t *graph) {
    reset_graph(graph);
}

long double _coal_mph_expected(coal_graph_node_t *node, size_t reward_index) {
    if (node->data.full_path_value >= 0) {
        return node->data.full_path_value;
    }

    double rate = 0;
    double sum = 0;
    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        // Sum the rates
        rate += values[i].weight;
    }

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        weight_t prob = values[i].weight / rate;
        sum += prob * _coal_mph_expected((coal_graph_node_t*) values[i].node,
                reward_index);
    }

    if (rate != 0) {
        long double exp = 1 / rate;
        sum += exp * (long double)(node->data.state_vec)[reward_index];
    }

    node->data.full_path_value = sum;

    return sum;
}

long double coal_mph_expected(coal_graph_node_t *graph, size_t reward_index) {
    reset_graph(graph);
    return _coal_mph_expected(graph, reward_index);
}

gsl_matrix_long_double * _coal_mph_im_expected(coal_graph_node_t *node, size_t n1, size_t n2) {
    if (node->data.visits >= 1) {
        // TODO: This does not work if visits > 1...
        // We must pass on some probability parameter.
        return node->data.pointer;
    }

    weighted_edge_t *values = vector_get(node->edges);
    double rate = 0;

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        // Sum the rates
        rate += values[i].weight;
    }

    if (node->data.visits == 0) {
        node->data.pointer = gsl_matrix_long_double_calloc(n1+1, n2+1);

        for (size_t i = 0; i <= n1; ++i) {
            for (size_t j = 0; j <= n2; ++j) {
                long double reward;
                reward = ((im_state_t*)node->data.state)->mat1[i][j] +
                         ((im_state_t*)node->data.state)->mat2[i][j];

                if (rate != 0) {
                    long double exp = 1 / rate;
                    long double current = gsl_matrix_long_double_get(node->data.pointer, i, j);
                    gsl_matrix_long_double_set(node->data.pointer, i, j, current + exp * reward);
                }
            }
        }
    }

    node->data.visits++;

    gsl_matrix_long_double *scaled = gsl_matrix_long_double_alloc(n1+1, n2+1);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        weight_t prob = values[i].weight / rate;
        gsl_matrix_long_double *res = _coal_mph_im_expected((coal_graph_node_t*) values[i].node,
                                                            n1, n2);
        gsl_matrix_long_double_memcpy(scaled, res);
        gsl_matrix_long_double_scale(scaled, (double)prob);
        gsl_matrix_long_double_add(node->data.pointer, scaled);
    }

    gsl_matrix_long_double_free(scaled);

    return node->data.pointer;
}


gsl_matrix_long_double * coal_mph_im_expected(coal_graph_node_t *graph, size_t n1, size_t n2) {
    reset_graph(graph);

    return _coal_mph_im_expected(graph, n1, n2);
}

void coal_mph_cov_assign_vertex(coal_graph_node_t *node, size_t reward_index) {
    if (node->data.prob >= 0) {
        return;
    }

    long double rate = 0;

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        // Sum the rates
        rate += values[i].weight;
    }

    node->data.prob = 0;
    weighted_edge_t *reverse_values = vector_get(node->reverse_edges);

    for (size_t i = 0; i < vector_length(node->reverse_edges); i++) {
        weighted_edge_t *parent = &(reverse_values[i]);
        long double rate_parent = 0;

        weighted_edge_t *parent_values = vector_get(parent->node->edges);

        for (size_t j = 0; j < vector_length(parent->node->edges); j++) {
            // Sum the rates
            rate_parent += parent_values[j].weight;
        }

        coal_mph_cov_assign_vertex((coal_graph_node_t*) parent->node, reward_index);
        node->data.prob += reverse_values[i].weight/rate_parent * ((coal_graph_node_t*)(parent->node))->data.prob;
    }

    if (rate != 0) {
        node->data.vertex_exp = node->data.prob * ((long double)(node->data.state_vec)[reward_index] / rate);
    } else {
        node->data.vertex_exp = 0;
    }
}

void coal_mph_cov_assign_desc(coal_graph_node_t *node, size_t reward_index) {
    if (node->data.descendants_exp_sum >= 0) {
        return;
    }

    long double rate = 0;

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        // Sum the rates
        rate += values[i].weight;
    }

    node->data.descendants_exp_sum = 0;

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        coal_mph_cov_assign_desc((coal_graph_node_t*) values[i].node, reward_index);
        node->data.descendants_exp_sum += values[i].weight / rate * ((coal_graph_node_t*)(values[i].node))->data.descendants_exp_sum;
    }

    if (rate != 0) {
        node->data.descendants_exp_sum += ((long double)(node->data.state_vec)[reward_index] / rate);
    }
}

long double _coal_mph_cov(coal_graph_node_t *node) {
    if (node->data.visited) {
        return 0;
    }

    weighted_edge_t *values = vector_get(node->edges);
    long double sum = 0;
    sum += node->data.descendants_exp_sum * node->data.vertex_exp;
    node->data.visited = true;

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        sum += _coal_mph_cov((coal_graph_node_t*) values[i].node);
    }

    return sum;
}

long double coal_mph_cov(coal_graph_node_t *graph,
                     size_t reward_index_1,
                     size_t reward_index_2) {
    // TODO: Pre-calculate the entire vector of rewards
    long double sum = 0;

    reset_graph(graph);
    graph->data.prob = 1.0f;
    graph->data.vertex_exp = 0.0;
    coal_mph_cov_assign_vertex(get_abs_vertex(graph), reward_index_1);
    coal_mph_cov_assign_desc(graph, reward_index_2);
    sum += _coal_mph_cov(graph);

    reset_graph(graph);
    graph->data.prob = 1.0f;
    graph->data.vertex_exp = 0.0;
    coal_mph_cov_assign_vertex(get_abs_vertex(graph), reward_index_2);
    coal_mph_cov_assign_desc(graph, reward_index_1);
    sum += _coal_mph_cov(graph);

    sum -= coal_mph_expected(graph, reward_index_1) *
            coal_mph_expected(graph, reward_index_2);


    return sum;
}

void coal_mph_cov_assign_vertex_all(coal_graph_node_t *node, size_t m) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;

    long double rate = 0;

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        // Sum the rates
        rate += values[i].weight;
    }

    if (vector_length(node->reverse_edges) == 0) {
        // Starting vertex
        node->data.prob = 1.0f;
    } else {
        node->data.prob = 0.0f;
    }

    weighted_edge_t *reverse_values = vector_get(node->reverse_edges);

    for (size_t i = 0; i < vector_length(node->reverse_edges); i++) {
        weighted_edge_t *parent = &(reverse_values[i]);
        long double rate_parent = 0;

        coal_mph_cov_assign_vertex_all((coal_graph_node_t*) parent->node, m);
        weighted_edge_t *parent_values = vector_get(parent->node->edges);

        for (size_t j = 0; j < vector_length(parent->node->edges); j++) {
            // Sum the rates
            rate_parent += parent_values[j].weight;
        }

        node->data.prob += reverse_values[i].weight/rate_parent * ((coal_graph_node_t*)(parent->node))->data.prob;
    }


    node->data.all_vertex_exp = gsl_vector_alloc(m);

    for (size_t j = 0; j < m; ++j) {
        if (rate != 0) {
            gsl_vector_set(node->data.all_vertex_exp, j, (double) (node->data.prob * (((double)(node->data.state_vec)[j]) / rate)));
        } else {
            gsl_vector_set(node->data.all_vertex_exp, j, 0);
        }
    }
}

void coal_mph_cov_assign_desc_all(coal_graph_node_t *node, size_t m) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;

    long double rate = 0;

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        // Sum the rates
        rate += values[i].weight;
    }

    node->data.all_descendants_exp_sum = gsl_vector_alloc(m);
    gsl_vector_set_zero(node->data.all_descendants_exp_sum);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        coal_mph_cov_assign_desc_all((coal_graph_node_t *) values[i].node, m);
    }

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        gsl_vector_axpby((const double) (values[i].weight / rate),
                         ((coal_graph_node_t*)(values[i].node))->data.all_descendants_exp_sum,
                         1, node->data.all_descendants_exp_sum);
    }

    gsl_vector *exp = gsl_vector_alloc(m);

    for (size_t j = 0; j < m; ++j) {
        if (rate != 0) {
            gsl_vector_set(exp, j, (double) ((node->data.state_vec)[j] / rate));
        } else {
            gsl_vector_set(exp, j, 0);
        }
    }

    gsl_vector_add(node->data.all_descendants_exp_sum, exp);

    gsl_vector_free(exp);
}

double **cov;

void _coal_mph_cov_all(coal_graph_node_t *node, size_t m) {
    if (node->data.visited) {
        return;
    }

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            //print_vector(stderr, node->data.state_vec, m);
            cov[i][j] += gsl_vector_get(node->data.all_descendants_exp_sum, i) *
                    gsl_vector_get(node->data.all_vertex_exp, j);
            cov[i][j] += gsl_vector_get(node->data.all_descendants_exp_sum, j) *
                         gsl_vector_get(node->data.all_vertex_exp, i);
        }
    }

    node->data.visited = true;

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        _coal_mph_cov_all((coal_graph_node_t*) values[i].node, m);
    }
}

double*** coal_mph_cov_all(coal_graph_node_t *graph, size_t m) {
    coal_graph_node_t *abs = get_abs_vertex(graph);
    coal_graph_reset(graph);
    coal_mph_cov_assign_vertex_all(abs, m);
    reset_graph_visited(graph);
    coal_mph_cov_assign_desc_all(graph, m);

    cov = calloc(m, sizeof(double*));

    for (size_t i = 0; i < m; ++i) {
        cov[i] = calloc(m, sizeof(double));
    }

    reset_graph_visited(graph);
    _coal_mph_cov_all(graph, m);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            cov[i][j] -= gsl_vector_get(graph->data.all_descendants_exp_sum, i) *
                    gsl_vector_get(graph->data.all_descendants_exp_sum, j);
        }
    }

    return &cov;
}


/*
 * Also ensures that the absorbing vertex has index 0
 */
int coal_label_vertex_index(size_t *largest_index, coal_graph_node_t *graph) {
    coal_graph_node_t *abs_vertex;
    abs_vertex = get_abs_vertex(graph);
    reset_graph_visited(graph);
    queue_t *queue;
    queue_create(&queue, 8);
    size_t index = 0;

    if (abs_vertex == NULL) {
        DIE_ERROR(1, "Absorbing vertex was not found\n");
    }

    // The absorbing vertex should have index 0
    queue_enqueue(queue, abs_vertex);

    queue_enqueue(queue, graph);

    while(!queue_empty(queue)) {
        coal_graph_node_t *node = queue_dequeue(queue);

        if (node->data.visited) {
            continue;
        }

        node->data.visited = true;

        weighted_edge_t *values = vector_get(node->edges);

        node->data.vertex_index = index++;

        for (size_t i = 0; i < vector_length(node->edges); i++) {
            queue_enqueue(queue, values[i].node);
        }
    }

    if (largest_index != NULL) {
        *largest_index = index - 1;
    }

    return 0;
}

size_t coal_get_edges(coal_graph_node_t *graph) {
    reset_graph_visited(graph);
    queue_t *queue;
    queue_create(&queue, 8);
    size_t edges = 0;

    queue_enqueue(queue, graph);

    while(!queue_empty(queue)) {
        coal_graph_node_t *node = queue_dequeue(queue);

        if (node->data.visited) {
            continue;
        }

        node->data.visited = true;

        edges += vector_length(node->edges);

        weighted_edge_t *values = vector_get(node->edges);

        for (size_t i = 0; i < vector_length(node->edges); i++) {
            queue_enqueue(queue, values[i].node);
        }
    }

    return edges;
}


static void _print_graph_list_im(FILE *stream, coal_graph_node_t *node,
                                 bool indexed,
                                 size_t vec_length, size_t vec_spacing, size_t n1, size_t n2) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;

    fprintf(stream, "Node: ");
    print_vector_spacing(stream, node->data.state_vec,
                         vec_length, vec_spacing);
    if (indexed) {
        fprintf(stream, " (%zu)", node->data.vertex_index);
    }
    fprintf(stream, ":\n");

    weighted_edge_t *values = vector_get(node->edges);

    if (node->data.state != NULL) {
        for (size_t i = 0; i < vector_length(node->edges); i++) {
            coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;
            char *type;
            size_t ppop1, ppop2, cpop1, cpop2;
            bool iso;

            ppop1 = count_lineages_mat(((im_state_t *) node->data.state)->mat1, n1, n2);
            ppop2 = count_lineages_mat(((im_state_t *) node->data.state)->mat2, n1, n2);
            cpop1 = count_lineages_mat(((im_state_t *) child->data.state)->mat1, n1, n2);
            cpop2 = count_lineages_mat(((im_state_t *) child->data.state)->mat2, n1, n2);
            iso = ((im_state_t *) child->data.state)->in_iso;

            if (iso) {
                if (cpop1 < ppop1) {
                    if (cpop2 > ppop2) {
                        type = "M12";
                    } else {
                        type = "C1";
                    }
                } else if (cpop2 < ppop2) {
                    if (cpop1 > ppop1) {
                        type = "M21";
                    } else {
                        type = "C2";
                    }
                } else {
                    type = "NONE";
                }
            } else {
                type = "Cm";
            }

            fprintf(stream, "\t");
            fprintf(stream, "(%s) (%Lf) ", type, values[i].weight);
            print_vector_spacing(stream, child->data.state_vec,
                                 vec_length, vec_spacing);
            if (indexed) {
                fprintf(stream, " (%zu)", child->data.vertex_index);
            }

            if (vector_length(child->edges) == 0) {
                fprintf(stream, " (abs)");
            }
            fprintf(stream, "\n");
        }
    }

    fprintf(stream, "\n");
    for (size_t i = 0; i < vector_length(node->edges); i++) {
        _print_graph_list_im(stream, (coal_graph_node_t *) values[i].node,
                          indexed, vec_length, vec_spacing, n1, n2);
    }
}

void coal_print_graph_list_im(FILE *stream, coal_graph_node_t *graph,
                              bool indexed,
                              size_t vec_length, size_t vec_spacing,
                              size_t n1, size_t n2) {
    if (indexed) {
        size_t largest_index;
        coal_label_vertex_index(&largest_index, graph);
    }
    reset_graph_visited(graph);
    _print_graph_list_im(stream, graph, indexed, vec_length, vec_spacing, n1, n2);
    fflush(stream);
}



void print_mat(const gsl_matrix_long_double *M) {
    for (size_t i = 0; i < M->size1; i++) {
        for (size_t j = 0; j < M->size2; j++) {
            fprintf(stdout, "%Lf ", gsl_matrix_long_double_get(M, i, j));
        }

        fprintf(stdout, "\n");
    }
}

void convert_mat(gsl_matrix *out, const gsl_matrix_long_double *M) {
    for (size_t i = 0; i < M->size1; ++i) {
        for (size_t j = 0; j < M->size2; ++j) {
            gsl_matrix_set(out, i, j, (double) (gsl_matrix_long_double_get(M, i, j)));
        }
    }
}


void sum_prob(long double* out, const gsl_matrix *mat_exp, coal_graph_node_t *node) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;

    if (node->data.vertex_index == 0) {
        // Don't use the absorbing vertex
        return;
    }

    // Use row 0 as we always start at vertex 0
    out[node->data.coals] += gsl_matrix_get(mat_exp, 0, (size_t)node->data.vertex_index-1);

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        sum_prob(out, mat_exp, (coal_graph_node_t *) values[i].node);
    }
}

int coal_im_get_number_coals_probs(long double **out,
                                       double isolation_time,
                                       const coal_gen_im_graph_args_t *args) {
    size_t limit = args->n1 + args->n2 - 1;
    *out = calloc(limit+1, sizeof(long double));
    coal_graph_node_t *graph;
    coal_gen_im_graph_args_t *prob_args = malloc(sizeof(coal_gen_im_graph_args_t));
    memcpy(prob_args, args, sizeof(coal_gen_im_graph_args_t));
    prob_args->num_iso_coal_events = prob_args->n1 + prob_args->n2 - 1;
    coal_gen_im_ss_graph(&graph, *prob_args);

    gsl_matrix_long_double *mat;
    coal_graph_as_gsl_mat(&mat, graph, false);
    gsl_matrix_long_double_scale(mat, isolation_time);
    gsl_matrix *mat_exp = gsl_matrix_alloc(mat->size1, mat->size2);
    gsl_matrix *mat2 = gsl_matrix_alloc(mat->size1, mat->size2);
    convert_mat(mat2, mat);
    gsl_linalg_exponential_ss(mat2, mat_exp, GSL_PREC_DOUBLE);

    reset_graph_visited(graph);
    sum_prob(*out, mat_exp, graph);

    long double sum = 0;

    for (size_t coals = 0; coals < limit; ++coals) {
        sum += (*out)[coals];
    }

    // The last prob is the rest
    (*out)[limit] = 1 - sum;

    return 0;
}

int coal_get_as_mat(gsl_matrix **S, const coal_graph_node_t *graph) {
    coal_graph_as_gsl_mat_double(S, (coal_graph_node_t *) graph, false);

    return 0;
}

int coal_get_mat_cdf(long double *out,
        double t,
        gsl_matrix *S,
        coal_graph_node_t *start) {
    gsl_matrix *alpha = gsl_matrix_calloc(1, S->size2);
    gsl_matrix *e = gsl_matrix_alloc(S->size1, 1);
    gsl_matrix_set_all(e, 1);

    for (size_t i = 2; i < S->size1; ++i) {
        weighted_edge_t *values = vector_get(start->edges);
        weighted_edge_t *edge = NULL;

        for (size_t j = 0; j < vector_length(start->edges); ++j) {
            if (((coal_graph_node_t *) values[j].node)->data.vertex_index == i) {
                edge = &(values[j]);
                break;
            }
        }

        if (edge != NULL) {
            gsl_matrix_set(alpha, 0, (size_t)((coal_graph_node_t*)edge->node)->data.vertex_index, (double)edge->weight);
        }
    }

    gsl_matrix *Sscaled = gsl_matrix_alloc(S->size1, S->size2);
    gsl_matrix_memcpy(Sscaled, S);
    gsl_matrix_scale(Sscaled, t);

    gsl_matrix *mat_exp = gsl_matrix_calloc(S->size1, S->size2);
    gsl_linalg_exponential_ss(Sscaled, mat_exp, GSL_PREC_DOUBLE);

    gsl_matrix *mat_mul_alpha = gsl_matrix_alloc(1, mat_exp->size2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, alpha, mat_exp, 0, mat_mul_alpha);

    gsl_matrix *mat_mul_e = gsl_matrix_alloc(1, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mat_mul_alpha, e, 0, mat_mul_e);

    *out = 1-gsl_matrix_get(mat_mul_e, 0, 0);
    gsl_matrix_free(mat_mul_e);
    gsl_matrix_free(mat_mul_alpha);
    gsl_matrix_free(mat_exp);
    gsl_matrix_free(e);
    gsl_matrix_free(alpha);
    gsl_matrix_free(Sscaled);

    return 0;
}

static int _coal_rewards_set(coal_graph_node_t *node, double(*reward_function)(coal_graph_node_t *node)) {
    if (node->data.visited) {
        return 0;
    }

    node->data.visited = true;

    node->data.reward = reward_function(node);

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); ++i) {
        _coal_rewards_set((coal_graph_node_t *) values[i].node, reward_function);
    }

    return 0;
}

int coal_rewards_set(coal_graph_node_t *graph, double(*reward_function)(coal_graph_node_t *node)) {
    reset_graph_visited(graph);
    _coal_rewards_set(graph, reward_function);
    return 0;
}

int _coal_reward_transform(coal_graph_node_t *node) {
    if (node->data.visited) {
        return 0;
    }

    node->data.visited = true;

    if (node->data.vertex_index == 0) {
        // It is the absorbing vertex
        return 0;
    }

    weighted_edge_t *values;

    bool found_unvisited;

    do {
        found_unvisited = false;
        values = vector_get(node->edges);

        for (size_t i = 0; i < vector_length(node->edges); ++i) {
            coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

            if (!child->data.visited) {
                _coal_reward_transform(child);
                found_unvisited = true;
                break;
            }
        }
    } while (found_unvisited);

    values = vector_get(node->edges);

    if (node->data.reward == 0) {
        weighted_edge_t *parents = vector_get(node->reverse_edges);

        for (ssize_t j = vector_length(node->reverse_edges)-1; j >= 0; --j) {
            coal_graph_node_t *my_parent = (coal_graph_node_t *) parents[j].node;
            weight_t total_weight = 0;

            for (size_t i = 0; i < vector_length(node->edges); ++i) {
                total_weight += values[i].weight;
            }

            // Take all my edges and add to my parent instead.
            for (ssize_t i = vector_length(node->edges)-1; i >= 0; --i) {
                weighted_edge_t my_child_edge = values[i];
                weight_t prob = my_child_edge.weight / total_weight;

                weighted_edge_t *edge_to_me = NULL;
                weighted_edge_t *parent_edges = vector_get(my_parent->edges);

                for (size_t k = 0; k < vector_length(my_parent->edges); ++k) {
                    if ((coal_graph_node_t *) parent_edges[k].node == node) {
                        edge_to_me = &(parent_edges[k]);
                        break;
                    }
                }

                if (edge_to_me == NULL) {
                    DIE_ERROR(1, "No edge found from parent to me. This should not happen");
                }

                //fprintf(stderr, "\"Adding\" edge from %zd to %zd. I am node %zd\n", my_parent->data.vertex_index,
                //        ((coal_graph_node_t*)my_child_edge.node)->data.vertex_index, node->data.vertex_index);
                graph_combine_edge((graph_node_t *) my_parent,
                        my_child_edge.node,
                        edge_to_me->weight * prob);
            }

            // Remove us from the parent's edges
            //fprintf(stderr, "Parent: Removing edge from %zd to %zd\n", my_parent->data.vertex_index,
            //        node->data.vertex_index);
            graph_remove_edge((graph_node_t *) my_parent, (graph_node_t *) node);
            parents = vector_get(node->reverse_edges);
        }


        for (ssize_t i = vector_length(node->edges)-1; i >= 0; --i) {
            weighted_edge_t my_child_edge = values[i];
            // Remove our edges to our children
            //fprintf(stderr, "Child: Removing edge from %zd to %zd\n", node->data.vertex_index,
            //        ((coal_graph_node_t *) my_child_edge.node)->data.vertex_index);
            graph_remove_edge((graph_node_t *) node, my_child_edge.node);
            values = vector_get(node->edges);
        }
    } else {
        //fprintf(stderr, "Reward transforming %zd by reward %f\n", node->data.vertex_index,
        //        node->data.reward);

        for (size_t i = 0; i < vector_length(node->edges); ++i) {
            values[i].weight /= node->data.reward;
        }
    }

    return 0;
}

/*
 * Note: adds a new starting vertex.
 * This holds the initial probability vector values.
 * A pass must have been made that sets the node.reward values
 */
int coal_reward_transform(coal_graph_node_t *graph, coal_graph_node_t **start) {
    coal_graph_node_create(start, NULL, NULL);
    (*start)->data.reward = 1;
    (*start)->data.reset_int = graph->data.reset_int;

    graph_add_edge((graph_node_t*) *start, (graph_node_t*) graph, 1);
    size_t largest;

    coal_label_vertex_index(&largest, *start);
    reset_graph_visited(*start);
    _coal_reward_transform(graph);

    return 0;
}

coal_graph_node_t * find_similar_vertex(const vec_entry_t *state, const size_t state_length,
        const avl_vec_node_t *non_iso_bst) {
    return avl_vec_find(non_iso_bst, state, state_length)->entry;
}

int _coal_graph_im_redirect_at_coals(coal_graph_node_t *node, const size_t coals,
                                     const avl_vec_node_t *non_iso_bst) {
    if (node->data.visited) {
        return 0;
    }

    node->data.visited = true;

    bool changes = false;

    do {
        changes = false;
        weighted_edge_t *values = vector_get(node->edges);

        for (size_t i = 0; i < vector_length(node->edges); i++) {
            coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

            if (!((im_state_t*)child->data.state)->in_iso) {
                // This node is from the non-isolated part
                continue;
            }

            if (child->data.coals == coals && child->data.type == IM_VERTEX_TYPE_COAL) {
                // We have found a node that should no longer exist
                im_state_t *state = child->data.state;
                im_state_t *cpy_state;
                im_state_cpy(&cpy_state, state, state->n1, state->n2);
                cpy_state->flag_mig1to2 = false;
                cpy_state->flag_mig2to1 = false;
                cpy_state->in_iso = false;
                vec_entry_t **combined_mat;
                vec_entry_t **empty_mat;

                combine_im_matrix(&combined_mat, (const vec_entry_t **) state->mat1,
                                  (const vec_entry_t **) state->mat2, state->n1, state->n2);
                empty_im_matrix(&empty_mat, state->n1, state->n2);

                cpy_state->mat1 = combined_mat;
                cpy_state->mat2 = empty_mat;

                vec_entry_t *vec_state;
                im_state_as_vec(&vec_state, cpy_state, cpy_state->n1, cpy_state->n2);

                coal_graph_node_t *found = find_similar_vertex(vec_state, im_state_length(cpy_state->n1, cpy_state->n2),
                                                               non_iso_bst);

                if (found == NULL) {
                    DIE_ERROR(1, "Expected to find a similar vertex\n");
                }

                if (child == found) {
                    DIE_ERROR(1, "Found vertex was the same as the child\n");
                }

                //fprintf(stderr, "Combining/adding edge to %zu from %zu\n Removing to %zu\n",
                //        found->data.vertex_index, node->data.vertex_index, child->data.vertex_index);

                graph_combine_edge((graph_node_t *) node, (graph_node_t *) found, values[i].weight);
                graph_remove_edge((graph_node_t *) node, (graph_node_t *) child);
                changes = true;
                break;
            } else {
                _coal_graph_im_redirect_at_coals(child, coals, non_iso_bst);
            }
        }
    } while (changes);
}

int coal_graph_im_redirect_at_coals(coal_graph_node_t *graph, const size_t coals, const avl_vec_node_t *non_iso_bst) {
    coal_label_vertex_index(NULL, graph);
    coal_graph_reset_visited(graph);

    return _coal_graph_im_redirect_at_coals(graph, coals, non_iso_bst);
}

int _coal_graph_clone(coal_graph_node_t *node) {
    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); ++i) {
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

        if (!child->data.visited) {
            coal_graph_node_t *cloned;
            coal_graph_node_create(&cloned, child->data.state_vec, child->data.state);

            child->data.pointer = cloned;

            if (cloned == NULL) {
                DIE_ERROR(1, "Could not create node\n");
            }

            child->data.visited = true;
            _coal_graph_clone(child);
        }

        graph_add_edge(node->data.pointer, child->data.pointer, values[i].weight);
    }

    return 0;
}

int coal_graph_clone(coal_graph_node_t **out, coal_graph_node_t *graph) {
    coal_graph_node_create(out, graph->data.state_vec, graph->data.state);
    coal_label_vertex_index(NULL, graph);
    reset_graph_visited(graph);
    graph->data.visited = true;
    graph->data.pointer = *out;
    _coal_graph_clone(graph);

    return 0;
}

void set_indegree(coal_graph_node_t *node) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;
    node->data.value = vector_length(node->reverse_edges);

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); ++i) {
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

        set_indegree(child);
    }
}

int coal_label_topological_index(size_t *largest_index, coal_graph_node_t *graph) {
    reset_graph_visited(graph);
    set_indegree(graph);
    reset_graph_visited(graph);
    queue_t *queue;
    queue_create(&queue, 8);
    queue_enqueue(queue, graph);
    size_t index = 0;

    while(!queue_empty(queue)) {
        coal_graph_node_t *node = queue_dequeue(queue);

        node->data.vertex_index = index;
        index++;

        weighted_edge_t *values = vector_get(node->edges);

        for (size_t i = 0; i < vector_length(node->edges); ++i) {
            coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

            child->data.value -= 1;

            if (child->data.value == 0 && !child->data.visited) {
                child->data.visited = true;
                queue_enqueue(queue, child);
            }
        }
    }

    if (largest_index != NULL) {
        *largest_index = index - 1;
    }
}

int _coal_construct_unshifted_discrete(coal_graph_node_t *node,
        double theta) {
    if (node->data.visited) {
        return 0;
    }

    node->data.visited = true;

    weight_t total_reward = node->data.reward;

    weight_t rate = 0;
    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        rate += values[i].weight;
    }


    weight_t reward_rate = rate / total_reward;
    weight_t self_loop_prob = (theta/2)/(theta/2 + reward_rate);

    for (size_t i = 0; i < vector_length(node->edges); ++i) {
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

        _coal_construct_unshifted_discrete(child, theta);
    }

    size_t n_edges = vector_length(node->edges);

    for (size_t i = 0; i < n_edges; ++i) {
        weight_t probability = values[i].weight / rate;
        values[i].weight = probability * (1 - self_loop_prob);
    }

    return 0;
}

int coal_construct_unshifted_discrete(coal_graph_node_t *graph, double theta) {
    reset_graph_visited(graph);
    return _coal_construct_unshifted_discrete(graph, theta);
}

int _coal_unshifted_discrete_apply_weighting(coal_graph_node_t *node, size_t *weights, size_t weight_length) {
    if (node->data.visited) {
        return 0;
    }

    node->data.visited = true;

    if (vector_length(node->edges) == 0) {
        return 0;
    }

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); ++i) {
        coal_graph_node_t *child = (coal_graph_node_t *) values[i].node;

        _coal_unshifted_discrete_apply_weighting(child, weights, weight_length);
    }

    size_t total_reward = 0;

    for (size_t i = 0; i < weight_length; ++i) {
        total_reward += node->data.state_vec[i];
    }

    weight_t total_weight = 0;

    for (size_t i = 0; i < weight_length; ++i) {
        total_weight += weights[i];
    }

    weight_t self_loop_prob = 1;

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        self_loop_prob -= values[i].weight;
    }

    coal_graph_node_t **helpers;

    size_t max_weight = 0;

    for (size_t i = 0; i < weight_length; ++i) {
        if (node->data.state_vec[i] != 0) {
            if (weights[i] > max_weight) {
                max_weight = weights[i];
            }
        }
    }

    helpers = calloc(max_weight + 1, sizeof(coal_graph_node_t*));
    helpers[0] = node;

    for (size_t i = 1; i < max_weight + 1; ++i) {
        coal_graph_node_create(&(helpers[i]), node->data.state_vec, NULL);
        graph_add_edge((graph_node_t *) helpers[i], (graph_node_t *) helpers[i - 1], 1);
    }

    for (size_t i = 0; i < weight_length; ++i) {
        if (weights[i] == 0 || node->data.state_vec[i] == 0) {
            continue;
        }

        weight_t weight_prob = ((weight_t)node->data.state_vec[i]) / total_reward;
        weight_t new_weight = self_loop_prob * weight_prob;

        if (weights[i] != 0) {
            if (weights[i] != 1) {
                graph_combine_edge((graph_node_t *) node, (graph_node_t *) helpers[weights[i] - 1], new_weight);
            }
        }
    }

    return 0;
}

int coal_unshifted_discrete_apply_weighting(coal_graph_node_t *graph, size_t *weights, size_t weight_length) {
    coal_label_vertex_index(NULL, graph);
    reset_graph_visited(graph);
    return _coal_unshifted_discrete_apply_weighting(graph, weights, weight_length);
}
