#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include "coal.h"
#include "bbst.h"
#include "utils.h"

// TODO: make the entry a pointer not an offset

typedef struct bst_node bst_node_t;

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
    (*out)->data.reset_flip = false;
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

        fprintf(stderr, "Edge weight %f:\n ", values[i].weight);
        print_graph_node((coal_graph_node_t*)values[i].node,vec_length, indent+1);
        fprintf(stderr, "\n");
    }
}

static void _print_graph_list(FILE *stream, coal_graph_node_t *node, size_t vec_length) {
    if (node->data.visited) {
        return;
    }

    node->data.visited = true;

    fprintf(stream, "Node: ");
    print_vector(stream, node->data.state_vec, vec_length);
    fprintf(stream, ":\n");

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        fprintf(stream, "\t");
        fprintf(stream, "(%f) ", values[i].weight);
        print_vector(stream,((coal_graph_node_t*)values[i].node)->data.state_vec, vec_length);
        fprintf(stream, "\n");
    }

    fprintf(stream, "\n");
    for (size_t i = 0; i < vector_length(node->edges); i++) {
        _print_graph_list(stream, (coal_graph_node_t *) values[i].node, vec_length);
    }
}

void coal_print_graph_list(FILE *stream, coal_graph_node_t *graph,
                           size_t vec_length) {
    coal_graph_reset(graph);
    _print_graph_list(stream, graph, vec_length);
    fflush(stream);
}

static coal_graph_node_t *get_abs_vertex(coal_graph_node_t *graph) {
    weighted_edge_t *values = vector_get(graph->edges);

    if (vector_length(graph->edges) == 0) {
        return graph;
    }

    for (size_t i = 0; i < vector_length(graph->edges); i++) {
        if (!((coal_graph_node_t*)values[i].node)->data.visited) {
            ((coal_graph_node_t*)values[i].node)->data.visited = true;
            coal_graph_node_t *res = get_abs_vertex((coal_graph_node_t*)values[i].node);

            if (res != NULL) {
                return res;
            }
        }
    }

    return NULL;
}

int coal_graph_as_phdist(phdist_t **phdist, coal_graph_node_t *graph) {
    queue_t *queue;
    queue_t *revisit_queue;
    coal_graph_reset(graph);
    queue_create(&queue, 8);
    queue_create(&revisit_queue, 8);
    size_t index = 0;
    coal_graph_node_t *abs_vertex = get_abs_vertex(graph);

    if (abs_vertex == NULL) {
        DIE_PERROR(1, "Cannot find absorbing vertex. Not a phase-type");
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

/*
 * Isolation with migration
 */
int im_state_init(im_state_t **out, size_t n1, size_t n2) {
    size_t state_size1 = n1+1;
    size_t state_size2 = n2+1;

    *out = malloc(sizeof(im_state_t*));
    (*out)->mat1 = malloc(sizeof(vec_entry_t*)*state_size1);
    (*out)->mat2 = malloc(sizeof(vec_entry_t*)*state_size1);

    for (size_t i = 0; i < state_size1; ++i) {
        (*out)->mat1[i] = calloc(state_size2, sizeof(vec_entry_t));
        (*out)->mat2[i] = calloc(state_size2, sizeof(vec_entry_t));
    }

    return 0;
}


int im_state_cpy(im_state_t **out,
                  im_state_t *in,
                  size_t n1, size_t n2) {
    im_state_init(out, n1, n2);

    for (size_t i = 0; i < n1 + 1; ++i) {
        for (size_t j = 0; j < n2 + 1; ++j) {
            (*out)->mat1[i][j] = in->mat1[i][j];
            (*out)->mat2[i][j] = in->mat2[i][j];
        }
    }

    return 0;
}

size_t im_state_length(size_t n1, size_t n2) {
    return (n1+1) * (n2+1) * 2;
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

    im_state_mat_as_vec(*out, state->mat1, n1+1, n2+1);
    im_state_mat_as_vec(*out+(n1+1)*(n2+1), state->mat2, n1+1, n2+1);

    return 0;
}

static int
im_visit_vertex(coal_graph_node_t **out, im_state_t *state,
                vec_entry_t *state_vec, avl_vec_node_t *bst,
                size_t vector_length,
                coal_gen_im_graph_args_t *args);

static inline int im_visit_coal_loop(coal_graph_node_t **out,
        vec_entry_t **d,
        coal_param_real_t scale,
        im_state_t *state, avl_vec_node_t *bst,
        size_t vector_length,
        coal_gen_im_graph_args_t *args) {
    if (scale == 0) {
        return 0;
    }

    for (vec_entry_t i1 = 0; i1 < args->n1+1; i1++) {
        for (vec_entry_t j1 = 0; j1 < args->n2+1; j1++) {
            for (vec_entry_t i2 = 0; i2 < args->n1+1; i2++) {
                for (vec_entry_t j2 = 0; j2 < args->n2+1; j2++) {
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

                    d[i1][j1]--;
                    d[i2][j2]--;
                    d[i1 + i2][j1 + j2]++;

                    vec_entry_t *v;
                    im_state_as_vec(&v, state, args->n1, args->n2);

                    coal_graph_node_t *new_vertex;
                    im_visit_vertex(&new_vertex,
                                    state, v,
                                    bst,
                                    vector_length,
                                    args);

                    d[i1][j1]++;
                    d[i2][j2]++;
                    d[i1 + i2][j1 + j2]--;

                    graph_add_edge((graph_node_t *) *out,
                            (graph_node_t *) new_vertex, t);
                }
            }
        }
    }

    return 0;
}

static inline int im_visit_mig_loop(coal_graph_node_t **out,
                                     vec_entry_t **d_from,
                                     vec_entry_t **d_to,
                                     coal_param_real_t scale,
                                     im_state_t *state, avl_vec_node_t *bst,
                                     size_t vector_length,
                                     coal_gen_im_graph_args_t *args) {
    if (scale == 0) {
        return 0;
    }

    for (vec_entry_t i = 0; i < args->n1+1; i++) {
        for (vec_entry_t j = 0; j < args->n2+1; j++) {
            // The MRCA cannot migrate
            if (i == args->n1 && j == args->n2) {
                continue;
            }

            mat_entry_t t;

            if (d_from[i][j] >= 1) {
                t = d_from[i][j];
            } else {
                continue;
            }

            t *= scale;

            d_from[i][j]--;
            d_to[i][j]++;

            vec_entry_t *v;
            im_state_as_vec(&v, state, args->n1, args->n2);

            coal_graph_node_t *new_vertex;
            im_visit_vertex(&new_vertex,
                            state, v,
                            bst,
                            vector_length,
                            args);

            d_from[i][j]++;
            d_to[i][j]--;

            graph_add_edge((graph_node_t *) *out,
                           (graph_node_t *) new_vertex, t);
        }
    }

    return 0;
}


static int
im_visit_vertex(coal_graph_node_t **out, im_state_t *state,
                vec_entry_t *state_vec, avl_vec_node_t *bst,
                size_t vector_length,
                coal_gen_im_graph_args_t *args) {
    avl_vec_node_t *bst_node = avl_vec_find(bst, state_vec, vector_length);

    if (bst_node != NULL) {
        *out = bst_node->entry;
        return 0;
    } else {
        im_state_t *vertex_state;
        im_state_cpy(&vertex_state, state, args->n1, args->n2);

        vec_entry_t *vertex_state_vec;
        im_state_as_vec(&vertex_state_vec, vertex_state, args->n1, args->n2);

        coal_graph_node_create(out, state_vec, vertex_state);

        avl_vec_insert(&bst, state_vec, *out, vector_length);

        im_visit_coal_loop(out, state->mat1, args->pop_scale1, state,
                bst, vector_length, args);
        im_visit_coal_loop(out, state->mat2, args->pop_scale2, state,
                           bst, vector_length, args);

        im_visit_mig_loop(out, state->mat1, state->mat2,
                args->pop_scale1 * args->mig_scale1 * args->migration_param,
                state, bst, vector_length, args);
        im_visit_mig_loop(out, state->mat2, state->mat1,
                          args->pop_scale2 * args->mig_scale2 * args->migration_param,
                          state, bst, vector_length, args);

        return 0;
    }
}

int coal_gen_im_graph(coal_graph_node_t **graph, coal_gen_im_graph_args_t args) {
    size_t n1 = args.n1;
    size_t n2 = args.n2;

    im_state_t *initial;
    im_state_init(&initial, n1, n2);
    initial->mat1[1][0] = n1;
    initial->mat2[0][1] = n2;

    vec_entry_t *initial_vec;
    im_state_as_vec(&initial_vec, initial, n1, n2);

    im_state_t *mrca;
    im_state_init(&mrca, n1, n2);
    mrca->mat1[n1][n2] = 1;
    mrca->mat2[n1][n2] = 1;

    vec_entry_t *mrca_vec;
    im_state_as_vec(&mrca_vec, mrca, n1, n2);

    coal_graph_node_t *absorbing_vertex;
    coal_graph_node_create(&absorbing_vertex, mrca_vec, mrca);

    avl_vec_node_t *BST;
    avl_vec_node_create(&BST, mrca_vec, absorbing_vertex, NULL);

    coal_graph_node_t *state_graph;

    im_visit_vertex(&state_graph, initial, initial_vec, BST,
                    im_state_length(n1, n2), &args);

    im_state_t *start_state;
    im_state_init(&start_state, n1, n2);


    vec_entry_t *start_state_vec;
    im_state_as_vec(&start_state_vec, start_state, n1, n2);

    coal_graph_node_t *start;
    coal_graph_node_create(&start, start_state_vec, start_state);
    graph_add_edge((graph_node_t*) start,
                   (graph_node_t*) state_graph, 1);

    *graph = start;

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
    phdist_t *reward;
    phdist_reward_transform(&reward, phdist);

    (*dist) = malloc(sizeof(d_dist_t));
    (*dist)->generator_fun = d_ph_gen_fun;
    d_phgen_args_t * args = malloc(sizeof(d_phgen_args_t));
    args->reward = reward->si_mat;
    (*dist)->args = args;




    return 0;
}

void _reset_graph(coal_graph_node_t *node, bool reset_flip) {
    // TODO: This visits multiple times
    if (node->data.reset_flip == reset_flip) {
        return;
    }

    node->data.reset_flip = !node->data.reset_flip;

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        _reset_graph((coal_graph_node_t *) values[i].node, reset_flip);
    }

    node->data.full_path_value = -1;
    node->data.vertex_exp = -1;
    node->data.prob = -1;
    node->data.descendants_exp_sum = -1;
    node->data.visited = false;
    node->data.pointer = NULL;
}

void reset_graph(coal_graph_node_t *node) {
    _reset_graph(node, !node->data.reset_flip);
}

void coal_graph_reset(coal_graph_node_t *graph) {
    reset_graph(graph);
}

double _coal_mph_expected(coal_graph_node_t *node, size_t reward_index) {
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
        double prob = values[i].weight / rate;
        sum += prob * _coal_mph_expected((coal_graph_node_t*) values[i].node,
                reward_index);
    }

    if (rate != 0) {
        double exp = 1 / rate;
        sum += exp * (double)(node->data.state_vec)[reward_index];
    }

    node->data.full_path_value = sum;

    return sum;
}

double coal_mph_expected(coal_graph_node_t *graph, size_t reward_index) {
    reset_graph(graph);
    return _coal_mph_expected(graph, reward_index);
}

void coal_mph_cov_assign_vertex(coal_graph_node_t *node, size_t reward_index) {
    if (node->data.prob >= 0) {
        return;
    }

    double rate = 0;

    weighted_edge_t *values = vector_get(node->edges);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        // Sum the rates
        rate += values[i].weight;
    }

    node->data.prob = 0;
    weighted_edge_t *reverse_values = vector_get(node->reverse_edges);

    for (size_t i = 0; i < vector_length(node->reverse_edges); i++) {
        weighted_edge_t *parent = &(reverse_values[i]);
        double rate_parent = 0;

        weighted_edge_t *parent_values = vector_get(parent->node->edges);

        for (size_t j = 0; j < vector_length(parent->node->edges); j++) {
            // Sum the rates
            rate_parent += parent_values[j].weight;
        }

        coal_mph_cov_assign_vertex((coal_graph_node_t*) parent->node, reward_index);
        node->data.prob += reverse_values[i].weight/rate_parent * ((coal_graph_node_t*)(parent->node))->data.prob;
    }

    if (rate != 0) {
        node->data.vertex_exp = node->data.prob * ((double)(node->data.state_vec)[reward_index] / rate);
    } else {
        node->data.vertex_exp = 0;
    }
}

void coal_mph_cov_assign_desc(coal_graph_node_t *node, size_t reward_index) {
    if (node->data.descendants_exp_sum >= 0) {
        return;
    }

    double rate = 0;

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
        node->data.descendants_exp_sum += ((double)(node->data.state_vec)[reward_index] / rate);
    }
}

double _coal_mph_cov(coal_graph_node_t *node) {
    if (node->data.visited) {
        return 0;
    }

    weighted_edge_t *values = vector_get(node->edges);
    double sum = 0;
    sum += node->data.descendants_exp_sum * node->data.vertex_exp;
    node->data.visited = true;

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        sum += _coal_mph_cov((coal_graph_node_t*) values[i].node);
    }

    return sum;
}

coal_graph_node_t *get_absorbing_vertex(coal_graph_node_t *graph) {
    if (vector_length(graph->edges) == 0) {
        return graph;
    }

    graph_node_t *first_child =
            ((weighted_edge_t*)(vector_get(graph->edges)))[0].node;

    return get_absorbing_vertex((coal_graph_node_t *) first_child);
}

double coal_mph_cov(coal_graph_node_t *graph,
                     size_t reward_index_1,
                     size_t reward_index_2) {
    // TODO: Pre-calculate the entire vector of rewards
    double sum = 0;

    reset_graph(graph);
    graph->data.prob = 1.0f;
    graph->data.vertex_exp = 0.0;
    coal_mph_cov_assign_vertex(get_absorbing_vertex(graph), reward_index_1);
    coal_mph_cov_assign_desc(graph, reward_index_2);
    sum += _coal_mph_cov(graph);

    reset_graph(graph);
    graph->data.prob = 1.0f;
    graph->data.vertex_exp = 0.0;
    coal_mph_cov_assign_vertex(get_absorbing_vertex(graph), reward_index_2);
    coal_mph_cov_assign_desc(graph, reward_index_1);
    sum += _coal_mph_cov(graph);

    sum -= coal_mph_expected(graph, reward_index_1) *
            coal_mph_expected(graph, reward_index_2);

    return sum;
}

int coal_label_vertex_index(size_t *largest_index, coal_graph_node_t *graph) {
    reset_graph(graph);
    queue_t *queue;
    queue_create(&queue, 8);
    size_t index = 0;
    queue_enqueue(queue, graph);

    while(!queue_empty(queue)) {
        coal_graph_node_t *node = queue_dequeue(queue);

        if (node->data.vertex_index != -1) {
            continue;
        }

        weighted_edge_t *values = vector_get(node->edges);

        node->data.vertex_index = index++;

        for (size_t i = 0; i < vector_length(node->edges); i++) {
            queue_enqueue(queue, values[i].node);
        }
    }

    *largest_index = index - 1;

    return 0;
}