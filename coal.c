#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_matrix_long_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "coal.h"
#include "bbst.h"
#include "utils.h"

static void reset_graph_visited(coal_graph_node_t *node);

typedef struct bst_node bst_node_t;
// TODO: make the entry a pointer not an offset
struct bst_node {
    bst_node_t *left;
    bst_node_t *right;
    vec_entry_t* value;
    size_t entry;
    coal_graph_node_t *graph_node;
};

static inline bool radix_cmp(const vec_entry_t* a, const vec_entry_t* b) {
    size_t i = 0;

    while (true) {
        if (a[i] < b[i]) {
            return false;
        } else if (a[i] > b[i]) {
            return true;
        }

        i++;
    }
}

static int bst_node_init(bst_node_t **bst_node, vec_entry_t *value, size_t entry) {
    if ((*bst_node = malloc(sizeof(bst_node_t))) == NULL) {
        return 1;
    }

    (*bst_node)->left = NULL;
    (*bst_node)->right = NULL;
    (*bst_node)->value = value;
    (*bst_node)->entry = entry;
    // TODO: Hack: We want to pass an actual pointer,
    // this is only backwards compatibility with previous
    // implementations
    (*bst_node)->graph_node = (coal_graph_node_t*) entry;
}

static bst_node_t * bst_init(vec_entry_t* new_value, size_t new_entry) {
    bst_node_t *n;
    bst_node_init(&n, new_value, new_entry);

    return n;
}

static bst_node_t * bst_add(bst_node_t * rootptr, vec_entry_t* new_value, size_t new_entry) {
    bst_node_t *n;
    bst_node_init(&n, new_value, new_entry);

    if (rootptr == NULL) {
        return n;
    }

    bst_node_t *node = rootptr;

    while (true) {
        if (radix_cmp(node->value, new_value)) {
            if (node->left == NULL) {
                node->left = n;

                break;
            } else {
                node = node->left;
            }
        } else {
            if (node->right == NULL) {
                node->right = n;

                break;
            } else {
                node = node->right;
            }
        }
    }

    return n;
}

static inline bool vector_eq(vec_entry_t* a, vec_entry_t* b,
        const size_t vec_nmemb) {
    return (memcmp(a, b, sizeof(vec_entry_t) * vec_nmemb)) == 0;
}

static bst_node_t * bst_find(bst_node_t *node, vec_entry_t *find_value,
                                const size_t vec_nmemb) {
    while (node != NULL) {
        if (vector_eq(node->value, find_value, vec_nmemb)) {
            return node;
        }

        if (radix_cmp(node->value, find_value)) {
            node = node->left;
        } else {
            node = node->right;
        }
    }

    return NULL;
}


// Returns true if found, returns false if potential new parent found
static bool bst_find_or_parent(bst_node_t **out,
        bst_node_t *node,
        vec_entry_t *find_value,
        const size_t vec_nmemb) {
    while (true) {
        if (vector_eq(node->value, find_value, vec_nmemb)) {
            *out = node;
            return true;
        }

        if (radix_cmp(node->value, find_value)) {
            if (node->left == NULL) {
                break;
            }

            node = node->left;
        } else {
            if (node->right == NULL) {
                break;
            }

            node = node->right;
        }
    }

    *out = node;
    return false;
}

struct stack {
    expanding_arr_t *stack;
    size_t head_index;
};

static int stack_init(struct stack *stack, size_t nmemb) {
    if ((stack->stack = malloc(sizeof(expanding_arr_t))) == NULL) {
        return 1;
    }

    stack->stack->length = 1;
    stack->stack->entry_size = sizeof(bst_node_t*);
    stack->stack->value = malloc(sizeof(bst_node_t**));
    *stack->stack->value = malloc(sizeof(bst_node_t*) * 1);
    stack->head_index = 0;

    return 0;
}

static int stack_destroy(struct stack *stack) {
    free(*stack->stack->value);
    free(stack->stack->value);
    free(stack->stack);
    free(stack);
}


static void stack_push(struct stack *stack, bst_node_t *entry) {
    expanding_arr_fit(stack->stack, stack->head_index+2);
    ((bst_node_t**)*stack->stack->value)[stack->head_index] = entry;
    stack->head_index++;
}

static bst_node_t* stack_pop(struct stack *stack) {
    bst_node_t *entry = ((bst_node_t**)*stack->stack->value)[stack->head_index-1];
    stack->head_index--;

    return entry;
}

static int stack_empty(struct stack *stack) {
    return (stack->head_index == 0);
}

struct _queue {
    bst_node_t **_queue;
    bst_node_t **head;
    bst_node_t **tail;
    size_t size;
    bst_node_t **bound;
};

static int _queue_init(struct _queue *_queue, size_t nmemb) {
    if ((_queue->_queue = malloc(sizeof(bst_node_t*)*nmemb)) == NULL) {
        return 1;
    }

    _queue->head = _queue->_queue;
    _queue->tail = _queue->_queue;
    _queue->size = nmemb;
    _queue->bound = _queue->_queue + _queue->size;
    memset(_queue->_queue, 0, _queue->size *sizeof(bst_node_t*));

    return 0;
}

static int _queue_destroy(struct _queue *_queue) {
    free(_queue->_queue);
}

static void _queue_print(struct _queue *_queue) {
    fprintf(stderr, "_queue: ");

    bst_node_t **i = _queue->_queue;
    size_t j = 0;

    while (i < _queue->_queue + _queue->size) {
        if (*i != NULL) {
            if (i == _queue->head) {
                fprintf(stderr, "h%zu, ", (*i)->entry);
            } else if (i == _queue->tail) {
                fprintf(stderr, "t%zu, ", (*i)->entry);
            } else {
                fprintf(stderr, "%zu, ", (*i)->entry);
            }
        } else {
            if (i == _queue->head) {
                fprintf(stderr, "h_, ");
            } else if (i == _queue->tail) {
                fprintf(stderr, "t_, ");
            } else {
                fprintf(stderr, "_, ");
            }
        }

        i++;
        j++;
    }

    fprintf(stderr, "\n");
}

static void _queue_en_queue(struct _queue *_queue, bst_node_t *entry) {
    *(_queue->tail) = entry;

    _queue->tail++;

    if (_queue->tail == _queue->bound) {
        _queue->tail = _queue->_queue;
    }

    if (_queue->tail == _queue->head) {
        size_t head_dist = _queue->size - (_queue->head - _queue->_queue);
        size_t head_offset = _queue->head - _queue->_queue;
        size_t tail_offset = _queue->tail - _queue->_queue;
        _queue->_queue = realloc(_queue->_queue, _queue->size * 2 *sizeof(bst_node_t**));

        memcpy(_queue->_queue + head_offset + _queue->size, _queue->_queue + head_offset, head_dist * sizeof(bst_node_t**));

        _queue->tail = _queue->_queue + tail_offset;

        _queue->head = _queue->_queue + _queue->size + head_offset;
        _queue->size *= 2;
        _queue->bound = _queue->_queue + _queue->size;
    }
}

static bst_node_t* _queue_de_queue(struct _queue *_queue) {
    bst_node_t *entry = *(_queue->head);
    *(_queue->head) = NULL;
    _queue->head++;

    if (_queue->head == _queue->bound) {
        _queue->head = _queue->_queue;
    }

    return entry;
}

static int _queue_empty(struct _queue *_queue) {
    return (_queue->head == _queue->tail);
}

#define MAX(a, b) (a>=b) ? a : b



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


int coal_gen_erlang_phdist(phdist_t **phdist, size_t samples) {
    size_t n_rows = samples-1;
    size_t n_cols = samples-1;

    *phdist = malloc(sizeof(phdist_t));

    (*phdist)->n_rw_rows = n_rows;
    (*phdist)->n_rw_cols = 1;

    if (((*phdist)->rw_arr = malloc(sizeof(size_t*)*n_rows)) == NULL) {
        return -1;
    }

    if (((*phdist)->si_mat = malloc(sizeof(mat_t))) == NULL) {
        return -1;
    }

    if (((*phdist)->si_mat->rows = malloc(sizeof(avl_flat_tuple_t*)*n_rows)) == NULL) {
        return -1;
    }

    if (((*phdist)->si_mat->cols = malloc(sizeof(avl_flat_tuple_t*)*n_cols)) == NULL) {
        return -1;
    }

    (*phdist)->si_mat->n_rows = n_rows;
    (*phdist)->si_mat->n_cols = n_cols;

    if (((*phdist)->si_mat->max_row_keys = malloc(sizeof(size_t)*n_rows)) == NULL) {
        return -1;
    }

    if (((*phdist)->si_mat->max_col_keys = malloc(sizeof(size_t)*n_cols)) == NULL) {
        return -1;
    }

    for (size_t i = 0; i < n_rows; i++) {
        mat_entry_t n = (mat_entry_t)(samples-i);

        if (((*phdist)->rw_arr[i] = malloc(sizeof(size_t)*1)) == NULL) {
            return -1;
        }

        (*phdist)->rw_arr[i][0] = samples - i;

        if (i == n_rows - 1) {
            if (((*phdist)->si_mat->rows[i] = malloc(sizeof(avl_flat_tuple_t) * 2)) == NULL) {
                return -1;
            }

            (*phdist)->si_mat->max_row_keys[i] = i;

            (*phdist)->si_mat->rows[i][0].entry =  -n * (n - 1) / 2;
            (*phdist)->si_mat->rows[i][0].key = i;
            (*phdist)->si_mat->rows[i][1].key = 0;
            (*phdist)->si_mat->rows[i][1].entry = 0;
        } else {
            if (((*phdist)->si_mat->rows[i] = malloc(sizeof(avl_flat_tuple_t) * 3)) == NULL) {
                return -1;
            }

            (*phdist)->si_mat->max_row_keys[i] = i+1;

            (*phdist)->si_mat->rows[i][0].entry = -n * (n - 1) / 2;
            (*phdist)->si_mat->rows[i][0].key = i;
            (*phdist)->si_mat->rows[i][1].entry = n * (n - 1) / 2;
            (*phdist)->si_mat->rows[i][1].key = i + 1;
            (*phdist)->si_mat->rows[i][2].key = 0;
            (*phdist)->si_mat->rows[i][2].entry = 0;
        }
    }

    for (size_t i = 0; i < n_cols; i++) {
        mat_entry_t n = (mat_entry_t)(samples-i+1);
        mat_entry_t n2 = (mat_entry_t)(samples-i);

        (*phdist)->si_mat->max_col_keys[i] = i;

        if (i == 0) {
            if (((*phdist)->si_mat->cols[i] = malloc(sizeof(avl_flat_tuple_t) * 2)) == NULL) {
                return -1;
            }

            (*phdist)->si_mat->cols[i][0].entry =  -n2 * (n2 - 1) / 2;
            (*phdist)->si_mat->cols[i][0].key = i;
            (*phdist)->si_mat->cols[i][1].key = i+1;
            (*phdist)->si_mat->cols[i][1].entry = 0;
        } else {
            if (((*phdist)->si_mat->cols[i] = malloc(sizeof(avl_flat_tuple_t) * 3)) == NULL) {
                return -1;
            }

            (*phdist)->si_mat->cols[i][0].entry =  n * (n - 1) / 2;
            (*phdist)->si_mat->cols[i][0].key = i-1;
            (*phdist)->si_mat->cols[i][1].entry =  -n2 * (n2 - 1) / 2;
            (*phdist)->si_mat->cols[i][1].key = i;
            (*phdist)->si_mat->cols[i][2].key = i+1;
            (*phdist)->si_mat->cols[i][2].entry = 0;
        }
    }

    return 0;
}

struct _queue_data {
    struct _queue *_queue;
    size_t vector_length;
    size_t vector_size;
    void *state;
};

static int coal_make_phdist(phdist_t **phdist,
                            size_t n_rows,
                            size_t n_cols,
                            avl_mat_node_t **rows,
                            avl_mat_node_t **cols,
                            size_t **StSpM,
                            size_t ri,
                            size_t state_size) {
    *phdist = malloc(sizeof(phdist_t));

    mat_malloc(&((*phdist)->si_mat), n_rows, n_cols);
    mat_flatten((*phdist)->si_mat, rows, cols);

    // Remove first row/column as it is the MRCA state_vec
    for (size_t r = 0; r < (*phdist)->si_mat->n_rows; r++) {
        (*phdist)->si_mat->max_row_keys[r]--;

        if ((*phdist)->si_mat->rows[r][0].key == 0) {
            (*phdist)->si_mat->rows[r]++;
        }

        avl_flat_tuple_t *p = (*phdist)->si_mat->rows[r];

        while (p->entry != 0) {
            p->key--;
            p++;
        }
    }

    for (size_t c = 0; c < (*phdist)->si_mat->n_cols; c++) {
        (*phdist)->si_mat->max_col_keys[c]--;

        avl_flat_tuple_t *p = (*phdist)->si_mat->cols[c];

        while (p->entry != 0) {
            p->key--;
            p++;
        }
    }

    // Remove the first row
    (*phdist)->si_mat->rows++;
    // Remove the first column
    (*phdist)->si_mat->cols++;
    (*phdist)->si_mat->n_rows--;
    (*phdist)->si_mat->n_cols--;

    // Remove the MRCA state_vec
    (*phdist)->rw_arr = StSpM+1;

    // Remove the MRCA state_vec
    (*phdist)->n_rw_rows = ri - 1;
    // Do not include the n-ton, therefore state_size - 1
    (*phdist)->n_rw_cols = state_size;

    return 0;
}

struct state_hobolth {
    bst_node_t *BST;
    size_t ri;
    size_t n_rows;
    size_t n_cols;
    expanding_arr_t *expanding_rows;
    expanding_arr_t *expanding_cols;
    size_t StSpM_n;
    size_t **StSpM;
};

static int add_edge(expanding_arr_t *expanding_rows, expanding_arr_t *expanding_cols,
        size_t from, size_t to,
        mat_entry_t weight) {
    expanding_arr_fit(expanding_rows, from);
    expanding_arr_fit(expanding_rows, to);
    expanding_arr_fit(expanding_cols, from);
    expanding_arr_fit(expanding_cols, to);

    avl_mat_node_t ***rows = (avl_mat_node_t***) expanding_rows->value;
    avl_mat_node_t ***cols = (avl_mat_node_t***) expanding_cols->value;

    fprintf(stderr, "Adding %zu to %zu with weight %f\n",
            from, to, weight);
    if (avl_insert_or_inc(&((*rows)[from]), to, weight)) {
        return 1;
    }

    fprintf(stderr, "Adding %zu to %zu with weight %f\n",
            from, from, -weight);

    if (avl_insert_or_inc(&((*rows)[from]), from, -weight)) {
        return 1;
    }

    if (avl_insert_or_inc(&((*cols)[to]), from, weight)) {
        return 1;
    }

    if (avl_insert_or_inc(&((*cols)[from]), from, -weight)) {
        return 1;
    }

    return 0;
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

    *weights = gsl_matrix_long_double_alloc(size, size);

    insert_into_gsl_weight_mat(*weights, graph, include_absorbing);

    return 0;
}


int coal_graph_as_phdist_rw(phdist_t **phdist, coal_graph_node_t *graph) {
    queue_t *queue;
    queue_t *revisit_queue;
    coal_graph_reset(graph);
    queue_create(&queue, 8);
    queue_create(&revisit_queue, 8);
    size_t index = 0;
    coal_graph_node_t *abs_vertex = get_abs_vertex(graph);

    if (abs_vertex == NULL) {
        DIE_ERROR(1, "Cannot find absorbing vertex. Not a phase-type");
    }

    queue_enqueue(queue, abs_vertex);
    queue_enqueue(queue, graph);

    while(!queue_empty(queue)) {
        coal_graph_node_t *node = queue_dequeue(queue);

        if (node->data.vertex_index != -1) {
            continue;
        }

        queue_enqueue(revisit_queue, node);
        weighted_edge_t *values = vector_get(node->edges);

        node->data.vertex_index = index++;

        for (size_t i = 0; i < vector_length(node->edges); i++) {
            queue_enqueue(queue, values[i].node);
        }
    }

    size_t n_rows = index;
    size_t n_cols = index;
    avl_mat_node_t **rows = malloc(sizeof(avl_mat_node_t*)*n_rows);
    avl_mat_node_t **cols = malloc(sizeof(avl_mat_node_t*)*n_cols);

    for (size_t i = 0; i < n_rows; ++i) {
        rows[i] = NULL;
    }

    for (size_t i = 0; i < n_cols; ++i) {
        cols[i] = NULL;
    }

    while(!queue_empty(revisit_queue)) {
        coal_graph_node_t *node = queue_dequeue(revisit_queue);

        weighted_edge_t *values = vector_get(node->edges);

        for (size_t i = 0; i < vector_length(node->edges); i++) {
            avl_insert_or_inc(&(rows[node->data.vertex_index]),
                              (size_t)((coal_graph_node_data_t*)&values[i].node->data)->vertex_index,
                              values[i].weight);

            avl_insert_or_inc(&(rows[node->data.vertex_index]),
                              (size_t)node->data.vertex_index,
                              -values[i].weight);

            avl_insert_or_inc(&(cols[(size_t)((coal_graph_node_data_t*)&values[i].node->data)->vertex_index]),
                              (size_t)node->data.vertex_index,
                              values[i].weight);

            avl_insert_or_inc(&(cols[node->data.vertex_index]),
                              (size_t)node->data.vertex_index,
                              -values[i].weight);

            // TODO: Remove last row
        }
    }

    *phdist = malloc(sizeof(phdist_t));

    mat_malloc(&((*phdist)->si_mat), n_rows, n_cols);
    mat_flatten((*phdist)->si_mat, rows, cols);

    // Remove first row/column as it is the MRCA state_vec
    for (size_t r = 0; r < (*phdist)->si_mat->n_rows; r++) {
        (*phdist)->si_mat->max_row_keys[r]--;

        if ((*phdist)->si_mat->rows[r][0].key == 0) {
            (*phdist)->si_mat->rows[r]++;
        }

        avl_flat_tuple_t *p = (*phdist)->si_mat->rows[r];

        while (p->entry != 0) {
            p->key--;
            p++;
        }
    }

    for (size_t c = 0; c < (*phdist)->si_mat->n_cols; c++) {
        (*phdist)->si_mat->max_col_keys[c]--;

        avl_flat_tuple_t *p = (*phdist)->si_mat->cols[c];

        while (p->entry != 0) {
            p->key--;
            p++;
        }
    }

    // Remove the first row
    (*phdist)->si_mat->rows++;
    // Remove the first column
    (*phdist)->si_mat->cols++;
    (*phdist)->si_mat->n_rows--;
    (*phdist)->si_mat->n_cols--;

    // Remove the MRCA state_vec
    //(*phdist)->rw_arr = StSpM+1;

    // Remove the MRCA state_vec
    //(*phdist)->n_rw_rows = ri - 1;
    // Do not include the n-ton, therefore state_size - 1
    //(*phdist)->n_rw_cols = state_size;
}

static int _queue_pop_ss_hobolth(phdist_t **phdist, struct _queue_data *_queue_data, void *args) {
    vec_entry_t *v;
    size_t myidx = 0;
    size_t idx;
    bst_node_t *idxN;
    bst_node_t *myidxN;
    struct _queue *_queue = _queue_data->_queue;
    struct state_hobolth *state = _queue_data->state;
    coal_args_hobolth_t *targs = args;

    size_t state_size = targs->n;

    size_t vector_size = _queue_data->vector_size;
    const size_t vec_nmemb = _queue_data->vector_length;

    bst_node_t *BST = state->BST;
    size_t ri = state->ri;
    size_t n_rows = state->n_rows;
    size_t n_cols = state->n_cols;
    expanding_arr_t *expanding_rows = state->expanding_rows;
    expanding_arr_t *expanding_cols = state->expanding_cols;
    avl_mat_node_t ***rows = (avl_mat_node_t***)expanding_rows->value;
    avl_mat_node_t ***cols = (avl_mat_node_t***)expanding_cols->value;
    size_t StSpM_n = state->StSpM_n;
    size_t **StSpM = state->StSpM;

    while (_queue_empty(_queue) == 0) {
        myidxN = _queue_de_queue(_queue);

        myidx = myidxN->entry;
        v = (vec_entry_t *) malloc(vector_size);
        memcpy(v, myidxN->value, vector_size);

        for (vec_entry_t i = 0; i < state_size; i++) {
            for (vec_entry_t j = i; j < state_size; j++) {
                if (((i == j && v[i] >= 2) || (i != j && v[i] > 0 && v[j] > 0))) {
                    ssize_t t = i == j ? v[i] * (v[i] - 1) / 2 : v[i] * v[j];

                    v[i]--;
                    v[j]--;
                    v[(i + j + 2) - 1]++;

                    if (!bst_find_or_parent(&idxN, BST, v, vec_nmemb)) {
                        vec_entry_t *nv = (vec_entry_t *) malloc(vector_size);
                        memcpy(nv, v, vector_size);

                        idx = ri;

                        while (ri >= StSpM_n) {
                            // TODO: Failure
                            StSpM = realloc(StSpM, sizeof(size_t *) * StSpM_n * 2);

                            for (size_t k = StSpM_n; k < StSpM_n * 2; k++) {
                                StSpM[k] = malloc(sizeof(size_t) * state_size);
                            }

                            StSpM_n *= 2;
                        }

                        memcpy(StSpM[ri], nv, vector_size);
                        idxN = bst_add(idxN, nv, ri);
                        _queue_en_queue(_queue, idxN);
                        ri = ri + 1;
                    } else {
                        idx = idxN->entry;
                    }

                    v[i]++;
                    v[j]++;
                    v[(i + j + 2) - 1]--;

                    add_edge(expanding_rows, expanding_cols,
                                myidx, idx,
                                t);

                    n_rows = MAX(n_rows, myidx + 1);
                    n_cols = MAX(n_cols, idx + 1);
                }
            }
        }

        free(v);
    }

    coal_make_phdist(phdist, n_rows, n_cols, *rows, *cols, StSpM, ri, state_size);
}

int _queue_init_standard_vector(struct _queue_data **out, size_t state_size, vec_entry_t n) {
    *out = malloc(sizeof(struct _queue_data));
    struct _queue *_queue = malloc(sizeof(struct _queue));
    struct state_hobolth *state = malloc(sizeof(struct state_hobolth));

    expanding_arr_t *expanding_rows = malloc(sizeof(expanding_arr_t));
    expanding_arr_t *expanding_cols = malloc(sizeof(expanding_arr_t));

    avl_mat_node_t** rows;
    avl_mat_node_t** cols;

    rows = malloc(sizeof(avl_mat_node_t*) * 1);
    cols = malloc(sizeof(avl_mat_node_t*) * 1);

    expanding_rows->value = malloc(sizeof(avl_mat_node_t***));
    *expanding_rows->value = rows;
    expanding_cols->value = malloc(sizeof(avl_mat_node_t***));
    *expanding_cols->value = cols;
    expanding_rows->length = 1;
    expanding_cols->length = 1;
    expanding_rows->entry_size = sizeof(avl_mat_node_t*);
    expanding_cols->entry_size = sizeof(avl_mat_node_t*);

    rows[0] = malloc(sizeof(avl_mat_node_t) * 1);
    avl_node_create(&(rows[0]), 0, 0, NULL);

    cols[0] = malloc(sizeof(avl_mat_node_t) * 1);
    avl_node_create(&(cols[0]), 0, 0, NULL);

    _queue_init(_queue, 2);

    vec_entry_t *initial = (vec_entry_t*)calloc(n, sizeof(vec_entry_t));
    initial[0] = n;

    vec_entry_t *mrca = (vec_entry_t*)calloc(n, sizeof(vec_entry_t));
    mrca[n-1] = 1;

    size_t StSpM_n = 2;
    size_t **StSpM;

    if ((StSpM = malloc(sizeof(size_t*)*StSpM_n)) == NULL) {
        return 1;
    }

    StSpM[0] = malloc(sizeof(size_t) * state_size);
    StSpM[1] = malloc(sizeof(size_t) * state_size);

    bst_node_t *BST = bst_init(mrca, 0);
    bst_node_t *BST_initial = bst_add(BST, initial, 1);

    _queue_en_queue(_queue, BST);
    _queue_en_queue(_queue, BST_initial);
    memcpy(StSpM[0], mrca, state_size);
    memcpy(StSpM[1], initial, state_size);

    // Set n_rows and n_cols to 1, as the MRCA state_vec does not
    // increase these.
    state->n_rows = 1;
    state->n_cols = 1;
    // Set ri to 2 as we add two nodes: Initial and MRCA
    state->ri = 2;
    state->BST = BST;

    state->expanding_rows = expanding_rows;
    state->expanding_cols = expanding_cols;

    state->StSpM_n = StSpM_n;
    state->StSpM = StSpM;

    (*out)->_queue = _queue;
    (*out)->vector_size = sizeof(vec_entry_t) * n;
    (*out)->vector_length = n;
    (*out)->state = state;
}

int _queue_init_ss_hobolth(struct _queue_data **out, size_t state_size) {
    _queue_init_standard_vector(out, state_size, (vec_entry_t)state_size);
}

int coal_gen_phdist(phdist_t **phdist, size_t state_size) {
    struct _queue_data *_queue_data;
    _queue_init_ss_hobolth(&_queue_data, state_size);
    struct coal_args_hobolth_t args;
    args.n = state_size;

    _queue_pop_ss_hobolth(phdist, _queue_data, &args);

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
    avl_vec_node_t *bst_node = avl_vec_find(bst, state, vector_length);

    if (bst_node != NULL) {
        *out = bst_node->entry;
        return 0;
    } else {
        coal_graph_node_create(out, state, state);

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
                            DIE_ERROR(1, "Illegal mig param, was %zu\n", (size_t)args->migration_type);
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
                coal_param_real_t scale = args->pop_scale2 *  args->mig_scale2;

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
                            DIE_ERROR(1, "Illegal mig param, was %zu\n", (size_t)args->migration_type);
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

        return 0;
    }
}

int coal_gen_im_graph(coal_graph_node_t **graph, coal_gen_im_graph_args_t args) {
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

    *graph = state_graph;

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

int d_ph_gen_fun(double **out, size_t from, size_t to, void *args) {
    d_phgen_args_t *arguments = (d_phgen_args_t*)args;

    mat_t *P;
    mat_t *sub;
    mat_t *inv_P;
    mat_t *id;
    mat_t *scaled;
    mat_t *p;

    mat_mul_scalar(&scaled, arguments->reward, 2/arguments->theta);
    mat_identity(&id, arguments->reward->n_rows);
    mat_sub(&inv_P, id, scaled);
    mat_inv(&P, inv_P);

    /*mat_inv(&inv, arguments->reward);
    mat_mul_scalar(&scaled, inv, 2/arguments->theta);
    mat_identity(&id, arguments->reward->n_rows);
    mat_sub(&sub, id, scaled);
    mat_sub(&sub2, id, sub);
    mat_row_sums(&p, sub2);*/

    mat_t *row_sums;
    mat_t *ones;

    mat_row_sums(&row_sums, P);
    mat_row_sums(&ones, id);
    mat_sub(&p, ones, row_sums);

    mat_t *pi;
    pi= malloc(sizeof(mat_t));
    pi->n_rows = 1;
    pi->n_cols = P->n_cols;
    pi->rows = malloc(sizeof(avl_flat_tuple_t*)*1);
    pi->cols = malloc(sizeof(avl_flat_tuple_t*)*pi->n_cols);

    pi->max_row_keys = malloc(sizeof(size_t)*1);
    pi->max_col_keys = malloc(sizeof(size_t)*pi->n_cols);

    pi->rows[0] = malloc(sizeof(avl_flat_tuple_t)*2);
    pi->rows[0][0].entry = 1;
    pi->rows[0][0].key = 0;
    pi->rows[0][1].entry = 0;
    pi->max_row_keys[0] = 0;

    for (size_t c = 0; c < pi->n_cols; c++) {
        pi->cols[c] = malloc(sizeof(avl_flat_tuple_t)*2);
        pi->cols[c][0].entry = 0;
        pi->cols[c][1].entry = 0;
        pi->max_col_keys[c] = 0;
    }

    pi->cols[0][0].key = 0;
    pi->cols[0][0].entry = 1;

    *out = malloc(sizeof(double)*(to-from+1));

    mat_t *power;
    mat_t *inv_power;
    mat_clone(&inv_power, inv_P);

    for (size_t i = 0; i<from;i++) {
        mat_mult(&inv_power, inv_power, inv_P);
    }

    mat_t *res1;
    mat_t *res2;

    for (size_t i = from; i<=to;i++) {
        mat_inv(&power, inv_power);
        mat_mult(&res1, pi, power);
        mat_mult(&res2, res1, p);


        (*out)[i-from] = res2->rows[0][0].entry;

        mat_mult(&inv_power, inv_power, inv_P);
    }

    return 0;
}

int coal_seg_sites(d_dist_t **dist, phdist_t *phdist) {
    DIE_ERROR(1, "Not implemented\n");

    /*phdist_t *reward;
    phdist_reward_transform(&reward, phdist);

    (*dist) = malloc(sizeof(d_dist_t));
    (*dist)->generator_fun = d_ph_gen_fun;
    d_phgen_args_t * args = malloc(sizeof(d_phgen_args_t));
    args->reward = reward->si_mat;
    (*dist)->args = args;*/




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
    node->data.pointer = NULL;
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

long double _coal_mph_im_expected(coal_graph_node_t *node, size_t rw_i, size_t rw_j) {
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
        sum += prob * _coal_mph_im_expected((coal_graph_node_t*) values[i].node,
                                         rw_i, rw_j);
    }

    long double reward;
    reward = ((im_state_t*)node->data.state)->mat1[rw_i][rw_j] +
            ((im_state_t*)node->data.state)->mat2[rw_i][rw_j];

    if (rate != 0) {
        long double exp = 1 / rate;
        sum += exp * reward;
    }

    node->data.full_path_value = sum;

    return sum;
}


long double coal_mph_im_expected(coal_graph_node_t *graph, size_t reward_index_i, size_t reward_index_j) {
    reset_graph(graph);

    return _coal_mph_im_expected(graph, reward_index_i, reward_index_j);
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

void convert_mat_o(gsl_matrix_long_double *out, const gsl_matrix *M) {
    for (size_t i = 0; i < M->size1; ++i) {
        for (size_t j = 0; j < M->size2; ++j) {
            gsl_matrix_long_double_set(out, i, j, ((long double)gsl_matrix_get(M, i, j)));
        }
    }
}

void mat_mul(gsl_matrix_long_double *C, const gsl_matrix_long_double *A, const gsl_matrix_long_double *B) {
    // New
    /*for (size_t i = 0; i < A->size1; ++i) {
        for (size_t j = 0; j < B->size2; ++j) {
            long double r = 0;

            for (size_t k = 0; k < A->size2; ++k) {
                r += gsl_matrix_long_double_get(A, i, k) *
                        gsl_matrix_long_double_get(B, k, j);
            }

            gsl_matrix_long_double_set(C, i, j, r);
        }
    }

    return;*/
    //old
    gsl_matrix *Ad = gsl_matrix_alloc(A->size1, A->size2);
    gsl_matrix *Bd = gsl_matrix_alloc(B->size1, B->size2);
    gsl_matrix *Cd = gsl_matrix_alloc(A->size1, B->size2);

    
    convert_mat(Ad, A);
    convert_mat(Bd, B);

    //fprintf(stderr, "Max %Lf, %Lf\n",gsl_matrix_long_double_max(A),
    //        gsl_matrix_long_double_max(B));

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, Ad, Bd, 0.0, Cd);
    
    convert_mat_o(C, Cd);
}

void mat_combine(gsl_matrix_long_double *C, const gsl_matrix_long_double *A, const gsl_matrix_long_double *B) {
    for (size_t i = 0; i < A->size1; ++i) {
        for (size_t j = 0; j < A->size2; ++j) {
            gsl_matrix_long_double_set(C, i, j, gsl_matrix_long_double_get(A, i, j));
        }
    }

    for (size_t i = 0; i < B->size1; ++i) {
        for (size_t j = 0; j < B->size2; ++j) {
            gsl_matrix_long_double_set(C, i+A->size1, j+A->size2, gsl_matrix_long_double_get(B, i, j));
        }
    }
}

#define ADJ_SUM_ITERS 100

long double adj_sum(double to,
                    gsl_matrix_long_double *A,
                    gsl_matrix_long_double *B,
                    gsl_matrix_long_double *C1,
                    gsl_matrix_long_double *C2) {
    gsl_matrix_long_double *sum;
    gsl_matrix_long_double *At;
    gsl_matrix_long_double *L;
    gsl_matrix_long_double *AtmL;
    gsl_matrix_long_double *LmAt;
    gsl_matrix_long_double *tmpC1;
    gsl_matrix_long_double *tmpC2;

    sum = gsl_matrix_long_double_alloc(B->size1, B->size2);
    gsl_matrix_long_double_set_zero(sum);

    At = gsl_matrix_long_double_alloc(A->size1, A->size2);
    gsl_matrix_long_double_memcpy(At, A);
    gsl_matrix_long_double_scale(At, to);

    L = gsl_matrix_long_double_alloc(B->size1, B->size2);
    gsl_matrix_long_double_memcpy(L, B);

    AtmL = gsl_matrix_long_double_alloc(At->size1, L->size2);
    LmAt = gsl_matrix_long_double_alloc(L->size1, At->size2);

    tmpC1 = gsl_matrix_long_double_alloc(C1->size1, sum->size2);
    tmpC2 = gsl_matrix_long_double_alloc(C1->size1, C2->size2);

    for (size_t k = 1; k <= ADJ_SUM_ITERS; ++k) {
        mat_mul(AtmL, At, L);
        mat_mul(LmAt, L, At);

        gsl_matrix_long_double_memcpy(L, AtmL);
        gsl_matrix_long_double_sub(L, LmAt);

        gsl_matrix_long_double_scale(L, 1.0/((double)k+1));

        gsl_matrix_long_double_add(sum, L);
    }

    mat_mul(tmpC1, C1, sum);
    mat_mul(tmpC2, tmpC1, C2);
    gsl_matrix_long_double_scale(tmpC2, to);

    return gsl_matrix_long_double_get(tmpC2, 0, 0);
}

void get_mat_rowsum(gsl_matrix_long_double *row, const gsl_matrix_long_double *M) {
    for (size_t i = 0; i < M->size1; ++i) {
        long double sum = 0;

        for (size_t j = 0; j < M->size2; ++j) {
            sum -= gsl_matrix_long_double_get(M, i, j);
        }

        gsl_matrix_long_double_set(row, i, 0, sum);
    }
}


void get_first(gsl_matrix_long_double **out,
               const size_t coals, const coal_gen_im_graph_args_t *args) {
    coal_gen_im_pure_cutoff_graph_args_t args_c = {
            .n1 = args->n1,
            .n2 = args->n2,
            .left = args->n1 + args->n2 - coals,
            .migration_type = args->migration_type,
            .pop_scale1 = args->pop_scale1,
            .pop_scale2 = args->pop_scale2,
            .mig_scale1 = args->mig_scale1,
            .mig_scale2 = args->mig_scale2
    };

    coal_graph_node_t *graph;

    coal_gen_im_pure_cutoff_graph(&graph, args_c);

    coal_graph_as_gsl_mat(out, graph, false);
}

void get_next(gsl_matrix_long_double **out,
        const size_t h1, const size_t h2,
        const coal_gen_im_graph_args_t *args) {
    coal_gen_im_pure_cutoff_graph_args_t args_c = {
            .n1 = h1,
            .n2 = h2,
            .left = h1 + h2 - 1,
            .migration_type = args->migration_type,
            .pop_scale1 = args->pop_scale1,
            .pop_scale2 = args->pop_scale2,
            .mig_scale1 = args->mig_scale1,
            .mig_scale2 = args->mig_scale2
    };

    coal_graph_node_t *graph;

    coal_gen_im_pure_cutoff_graph(&graph, args_c);

    coal_graph_as_gsl_mat(out, graph, false);
}

int get_left_prob(gsl_matrix_long_double **out, size_t *correct_index,
                   const size_t H1, const size_t H2,
                   const coal_gen_im_graph_args_t *args) {
    coal_graph_node_t *graph;

    coal_gen_im_cutoff_graph_args_t args_c = {
            .n1 = args->n1,
            .n2 = args->n2,
            .left_n1 = H1,
            .left_n2 = H2,
            .migration_type = args->migration_type,
            .pop_scale1 = args->pop_scale1,
            .pop_scale2 = args->pop_scale2,
            .mig_scale1 = args->mig_scale1,
            .mig_scale2 = args->mig_scale2
    };

    coal_graph_node_t *correct;
    coal_gen_im_prob_vertex_graph(&graph, &correct, args_c);

    // We may have been given an unreachable state, return NULL then
    if (correct == NULL) {
        *out = NULL;
        return 0;
    }

    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, graph);
    gsl_matrix_long_double *res = gsl_matrix_long_double_alloc(size, size);

    for (size_t i = 1; i < size; ++i) {
        weight_t rowsum = 0;
        weight_t taken = 0;

        for (size_t j = 1; j < size; ++j) {
            if (i == j) {
                rowsum = -mat[i][j];
            } else {
                taken += mat[i][j];
            }
        }

        for (size_t j = 1; j < size; ++j) {
            if (i == j) {
                gsl_matrix_long_double_set(res, i - 1, j - 1, 0);
            } else {
                gsl_matrix_long_double_set(res, i - 1, j - 1, mat[i][j]/rowsum);
            }
        }

        gsl_matrix_long_double_set(res, i - 1, size-1, (rowsum - taken)/rowsum);
    }

    for (size_t j = 0; j < size-1; ++j) {
        gsl_matrix_long_double_set(res, size-1, j, 0);
    }

    gsl_matrix_long_double_set(res, size-1, size-1, 1);
    *correct_index = (size_t) correct->data.vertex_index;
    *out = res;

    return 0;
}


int get_mat_prob(long double *prob, const size_t i, const size_t j, const gsl_matrix_long_double *M) {
    gsl_matrix_long_double *Mcpy;

    Mcpy = gsl_matrix_long_double_alloc(M->size1, M->size2);
    gsl_matrix_long_double_memcpy(Mcpy, M);
    long double sum = 0;

    for (size_t k = 0; k < 200; ++k) {
        long double inc = gsl_matrix_long_double_get(Mcpy, i, j);
        sum += inc;
        mat_mul(Mcpy, M, Mcpy);
    }

    *prob = sum;
    return 0;
}

int coal_im_get_number_coals_prob(long double *out,
        size_t coals, double isolation_time,
        const coal_gen_im_graph_args_t *args) {
        gsl_matrix_long_double *first;

    if (coals == 0) {
        get_first(&first, 1, args);
        gsl_matrix_long_double *alpha = gsl_matrix_long_double_alloc(1, first->size1);
        gsl_matrix_long_double_set_zero(alpha);
        gsl_matrix_long_double_set(alpha, 0, 0, 1.0f);
        gsl_matrix_long_double *e = gsl_matrix_long_double_alloc(first->size1, 1);
        gsl_matrix_long_double_set_all(e, 1);

        gsl_matrix *timeFirstD = gsl_matrix_alloc(first->size1, first->size2);
        convert_mat(timeFirstD, first);
        gsl_matrix_scale(timeFirstD, isolation_time);

        gsl_matrix *expmD = gsl_matrix_alloc(first->size1, first->size2);
        gsl_linalg_exponential_ss(timeFirstD, expmD, GSL_PREC_DOUBLE);

        gsl_matrix_long_double *expm = gsl_matrix_long_double_alloc(first->size1, first->size2);
        convert_mat_o(expm, expmD);
        gsl_matrix_long_double *C1 = gsl_matrix_long_double_alloc(alpha->size1, expm->size2);
        mat_mul(C1, alpha, expm);
        gsl_matrix_long_double *C2 = gsl_matrix_long_double_alloc(C1->size1, e->size2);
        mat_mul(C2, C1, e);

        *out = gsl_matrix_long_double_get(C2, 0, 0);
        return 0;
    }

    if (coals == args->n1 + args->n2 - 1) {
        get_first(&first, coals, args);
        gsl_matrix_long_double *alpha = gsl_matrix_long_double_alloc(1, first->size1);
        gsl_matrix_long_double_set_zero(alpha);
        gsl_matrix_long_double_set(alpha, 0, 0, 1.0f);
        gsl_matrix_long_double *e = gsl_matrix_long_double_alloc(first->size1, 1);
        gsl_matrix_long_double_set_all(e, 1);

        gsl_matrix *timeFirstD = gsl_matrix_alloc(first->size1, first->size2);
        convert_mat(timeFirstD, first);
        gsl_matrix_scale(timeFirstD, isolation_time);

        gsl_matrix *expmD = gsl_matrix_alloc(first->size1, first->size2);
        gsl_linalg_exponential_ss(timeFirstD, expmD, GSL_PREC_DOUBLE);

        gsl_matrix_long_double *expm = gsl_matrix_long_double_alloc(first->size1, first->size2);
        convert_mat_o(expm, expmD);
        gsl_matrix_long_double *C1 = gsl_matrix_long_double_alloc(alpha->size1, expm->size2);
        mat_mul(C1, alpha, expm);
        gsl_matrix_long_double *C2 = gsl_matrix_long_double_alloc(C1->size1, e->size2);
        mat_mul(C2, C1, e);


        *out = 1-gsl_matrix_long_double_get(C2, 0, 0);
        return 0;
    }

    get_first(&first, coals, args);
    long double sum = 0;

    for (size_t h1 = 0; h1 <= args->n1 + args->n2 - coals; ++h1) {
        size_t h2 = args->n1 + args->n2 - coals - h1;

        gsl_matrix_long_double *m;
        size_t correct_index;
        get_left_prob(&m, &correct_index, h1, h2, args);

        // Unreachable state
        if (m == NULL) {
            continue;
        }

        long double prob;
        get_mat_prob(&prob, 0, correct_index - 1, m);

        if (prob < 0.00001) {
            continue;
        }

        gsl_matrix_long_double *other;
        get_next(&other, h1, h2, args);

        gsl_matrix_long_double *Sc = gsl_matrix_long_double_alloc(first->size1 + other->size1,
                                                                  first->size2 + other->size2);
        mat_combine(Sc, first, other);

        gsl_matrix_long_double *sc = gsl_matrix_long_double_alloc(Sc->size1, 1);
        get_mat_rowsum(sc, Sc);

        gsl_matrix_long_double *alpha1 = gsl_matrix_long_double_alloc(1, sc->size1);
        gsl_matrix_long_double *alpha2 = gsl_matrix_long_double_alloc(1, sc->size1);
        gsl_matrix_long_double_set_zero(alpha1);
        gsl_matrix_long_double_set_zero(alpha2);
        gsl_matrix_long_double_set(alpha1, 0, 0, 1.0f);
        gsl_matrix_long_double_set(alpha2, 0, first->size1, 1.0f);

        gsl_matrix_long_double *e = gsl_matrix_long_double_alloc(sc->size1, 1);
        gsl_matrix_long_double_set_all(e, 1);

        gsl_matrix_long_double *A = gsl_matrix_long_double_alloc(Sc->size1, Sc->size2);
        gsl_matrix_long_double_memcpy(A, Sc);
        gsl_matrix_long_double_scale(A, -1);

        gsl_matrix_long_double *B = gsl_matrix_long_double_alloc(e->size1, alpha1->size2);
        mat_mul(B, e, alpha1);

        gsl_matrix_long_double *timeSc = gsl_matrix_long_double_alloc(Sc->size1, Sc->size2);
        gsl_matrix_long_double_memcpy(timeSc, Sc);
        gsl_matrix_long_double_scale(timeSc, isolation_time);

        gsl_matrix *timeScD = gsl_matrix_alloc(Sc->size1, Sc->size2);
        convert_mat(timeScD, timeSc);

        gsl_matrix *expmD = gsl_matrix_alloc(Sc->size1, Sc->size2);
        gsl_linalg_exponential_ss(timeScD, expmD, GSL_PREC_DOUBLE);

        gsl_matrix_long_double *expm = gsl_matrix_long_double_alloc(Sc->size1, Sc->size2);
        convert_mat_o(expm, expmD);

        gsl_matrix_long_double *C1 = gsl_matrix_long_double_alloc(alpha2->size1, expm->size2);
        mat_mul(C1, alpha2, expm);

        long double adj = adj_sum(isolation_time, A, B, C1, sc);
        sum += adj * prob;
    }

    *out = sum;

    return 0;
}

int coal_im_get_number_coals_probs(long double **out,
        double isolation_time,
        const coal_gen_im_graph_args_t *args) {
    size_t limit = args->n1 + args->n2 - 1;
    *out = calloc(limit, sizeof(double));
    for (size_t coals = 0; coals <= limit; ++coals) {
        coal_im_get_number_coals_prob(&((*out)[coals]), coals, isolation_time, args);
    }

    // We have a bug somewhere, TODO: Fix this
    (*out)[1] = 0;
    long double sum = 0;

    for (size_t coals = 0; coals <= limit; ++coals) {
        if (isnan(((*out)[coals]))) {
            ((*out)[coals]) = 0;
        }

        if ((*out)[coals] < -0.01) {
            DIE_ERROR(1, "Probability less than 1 occurred (was %Lf)\n",
                      (*out)[coals]);
        }

        sum += ((*out)[coals]);
    }

    if (sum > 1.01) {
        DIE_ERROR(1, "Probability sum larger than 1 occurred (was %Lf)\n",
                  sum);
    }

    (*out)[1] = 1-sum;

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