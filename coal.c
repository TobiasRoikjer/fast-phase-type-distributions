#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include "coal.h"
#include "bbst.h"
#include "utils.h"

// TODO: make the entry a pointer not an offset

typedef size_t vec_entry_t;

typedef struct bst_node bst_node_t;

struct bst_node {
    bst_node_t *left;
    bst_node_t *right;
    vec_entry_t* value;
    size_t entry;
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

static size_t vec_nmemb;

static inline bool vector_eq(vec_entry_t* a, vec_entry_t* b) {
    return (memcmp(a, b, sizeof(vec_entry_t) * vec_nmemb)) == 0;
}

static bst_node_t * bst_find(bst_node_t *node, vec_entry_t *find_value) {
    while (node != NULL) {
        if (vector_eq(node->value, find_value)) {
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
static bool bst_find_or_parent(bst_node_t **out, bst_node_t *node, vec_entry_t *find_value) {
    while (true) {
        if (vector_eq(node->value, find_value)) {
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

struct queue {
    bst_node_t **queue;
    bst_node_t **head;
    bst_node_t **tail;
    size_t size;
    bst_node_t **bound;
};

static int queue_init(struct queue *queue, size_t nmemb) {
    if ((queue->queue = malloc(sizeof(bst_node_t*)*nmemb)) == NULL) {
        return 1;
    }

    queue->head = queue->queue;
    queue->tail = queue->queue;
    queue->size = nmemb;
    queue->bound = queue->queue + queue->size;
    memset(queue->queue, 0, queue->size *sizeof(bst_node_t*));

    return 0;
}

static int queue_destroy(struct queue *queue) {
    free(queue->queue);
}

static void queue_print(struct queue *queue) {
    fprintf(stderr, "Queue: ");

    bst_node_t **i = queue->queue;
    size_t j = 0;

    while (i < queue->queue + queue->size) {
        if (*i != NULL) {
            if (i == queue->head) {
                fprintf(stderr, "h%zu, ", (*i)->entry);
            } else if (i == queue->tail) {
                fprintf(stderr, "t%zu, ", (*i)->entry);
            } else {
                fprintf(stderr, "%zu, ", (*i)->entry);
            }
        } else {
            if (i == queue->head) {
                fprintf(stderr, "h_, ");
            } else if (i == queue->tail) {
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

static void queue_enqueue(struct queue *queue, bst_node_t *entry) {
    *(queue->tail) = entry;

    queue->tail++;

    if (queue->tail == queue->bound) {
        queue->tail = queue->queue;
    }

    if (queue->tail == queue->head) {
        size_t head_dist = queue->size - (queue->head - queue->queue);
        size_t head_offset = queue->head - queue->queue;
        size_t tail_offset = queue->tail - queue->queue;
        queue->queue = realloc(queue->queue, queue->size * 2 *sizeof(bst_node_t**));

        memcpy(queue->queue + head_offset + queue->size, queue->queue + head_offset, head_dist * sizeof(bst_node_t**));

        queue->tail = queue->queue + tail_offset;

        queue->head = queue->queue + queue->size + head_offset;
        queue->size *= 2;
        queue->bound = queue->queue + queue->size;
    }
}

static bst_node_t* queue_dequeue(struct queue *queue) {
    bst_node_t *entry = *(queue->head);
    *(queue->head) = NULL;
    queue->head++;

    if (queue->head == queue->bound) {
        queue->head = queue->queue;
    }

    return entry;
}

static int queue_empty(struct queue *queue) {
    return (queue->head == queue->tail);
}

#define MAX(a, b) (a>=b) ? a : b


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

struct queue_data {
    struct queue *queue;
    size_t vector_length;
    size_t vector_size;
    void *state;
};

static int coal_make_phdist(phdist_t **phdist,
                            size_t n_rows,
                            size_t n_cols,
                            avl_node_t **rows,
                            avl_node_t **cols,
                            size_t **StSpM,
                            size_t ri,
                            size_t state_size) {
    *phdist = malloc(sizeof(phdist_t));

    mat_malloc(&((*phdist)->si_mat), n_rows, n_cols);
    mat_flatten((*phdist)->si_mat, rows, cols);

    // Remove first row/column as it is the MRCA state
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

    // Remove the MRCA state
    (*phdist)->rw_arr = StSpM+1;

    // Remove the MRCA state
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

static int queue_pop_ss_hobolth_compress(phdist_t **phdist, struct queue_data *queue_data, void *args) {
    vec_entry_t *v;
    size_t myidx = 0;
    size_t idx;
    bst_node_t *idxN;
    bst_node_t *myidxN;
    struct queue *queue = queue_data->queue;

    // TODO: Quick hack, make a stack from the queue
    struct stack *stack = malloc(sizeof(struct stack));
    stack_init(stack, 1);
    bst_node_t *queue1 = queue_dequeue(queue);
    bst_node_t *queue2 = queue_dequeue(queue);
    stack_push(stack, queue2);
    stack_push(stack, queue1);

    struct state_hobolth *state = queue_data->state;
    coal_args_compress_t *targs = args;

    size_t state_size = targs->n;
    size_t reward_index = targs->reward_index;

    size_t vector_size = queue_data->vector_size;

    bst_node_t *BST = state->BST;
    size_t ri = state->ri;
    size_t n_rows = state->n_rows;
    size_t n_cols = state->n_cols;
    expanding_arr_t *expanding_rows = state->expanding_rows;
    expanding_arr_t *expanding_cols = state->expanding_cols;
    avl_node_t ***rows = (avl_node_t***)expanding_rows->value;
    avl_node_t ***cols = (avl_node_t***)expanding_cols->value;
    size_t StSpM_n = state->StSpM_n;
    size_t **StSpM = state->StSpM;

    while (!stack_empty(stack)) {
        myidxN = stack_pop(stack);

        myidx = myidxN->entry;
        v = (vec_entry_t *) malloc(vector_size);
        memcpy(v, myidxN->value, vector_size);

        for (vec_entry_t i = 0; i < state_size; i++) {
            for (vec_entry_t j = i; j < state_size; j++) {
                if (((i == j && v[i] >= 2) || (i != j && v[i] > 0 && v[j] > 0))) {
                    mat_entry_t t = i == j ? v[i] * (v[i] - 1) / 2 : v[i] * v[j];

                    if (v[reward_index] == 0) {
                        // TODO Push alpha
                        
                    } else {
                        //t /= v[reward_index];
                    }

                    v[i]--;
                    v[j]--;
                    v[(i + j + 2) - 1]++;

                    if (!bst_find_or_parent(&idxN, BST, v)) {
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
                        stack_push(stack, idxN);
                        ri = ri + 1;
                    } else {
                        idx = idxN->entry;
                    }

                    v[i]++;
                    v[j]++;
                    v[(i + j + 2) - 1]--;

                    expanding_arr_fit(expanding_rows, idx);
                    expanding_arr_fit(expanding_rows, myidx);
                    expanding_arr_fit(expanding_cols, idx);
                    expanding_arr_fit(expanding_cols, myidx);

                    if (avl_insert_or_inc(&((*rows)[myidx]), idx, t)) {
                        return 1;
                    }

                    if (avl_insert_or_inc(&((*rows)[myidx]), myidx, -t)) {
                        return 1;
                    }

                    if (avl_insert_or_inc(&((*cols)[idx]), myidx, t)) {
                        return 1;
                    }

                    if (avl_insert_or_inc(&((*cols)[myidx]), myidx, -t)) {
                        return 1;
                    }

                    n_rows = MAX(n_rows, myidx + 1);
                    n_cols = MAX(n_cols, idx + 1);
                }
            }
        }

        free(v);
    }

    coal_make_phdist(phdist, n_rows, n_cols, *rows, *cols, StSpM, ri, state_size);
}

static int queue_pop_ss_hobolth(phdist_t **phdist, struct queue_data *queue_data, void *args) {
    vec_entry_t *v;
    size_t myidx = 0;
    size_t idx;
    bst_node_t *idxN;
    bst_node_t *myidxN;
    struct queue *queue = queue_data->queue;
    struct state_hobolth *state = queue_data->state;
    coal_args_hobolth_t *targs = args;

    size_t state_size = targs->n;
    
    size_t vector_size = queue_data->vector_size;
    
    bst_node_t *BST = state->BST;
    size_t ri = state->ri;
    size_t n_rows = state->n_rows;
    size_t n_cols = state->n_cols;
    expanding_arr_t *expanding_rows = state->expanding_rows;
    expanding_arr_t *expanding_cols = state->expanding_cols;
    avl_node_t ***rows = (avl_node_t***)expanding_rows->value;
    avl_node_t ***cols = (avl_node_t***)expanding_cols->value;
    size_t StSpM_n = state->StSpM_n;
    size_t **StSpM = state->StSpM;

    while (queue_empty(queue) == 0) {
        myidxN = queue_dequeue(queue);

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

                    if (!bst_find_or_parent(&idxN, BST, v)) {
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
                        queue_enqueue(queue, idxN);
                        ri = ri + 1;
                    } else {
                        idx = idxN->entry;
                    }

                    v[i]++;
                    v[j]++;
                    v[(i + j + 2) - 1]--;

                    expanding_arr_fit(expanding_rows, idx);
                    expanding_arr_fit(expanding_rows, myidx);
                    expanding_arr_fit(expanding_cols, idx);
                    expanding_arr_fit(expanding_cols, myidx);

                    if (avl_insert_or_inc(&((*rows)[myidx]), idx, t)) {
                        return 1;
                    }

                    if (avl_insert_or_inc(&((*rows)[myidx]), myidx, -t)) {
                        return 1;
                    }

                    if (avl_insert_or_inc(&((*cols)[idx]), myidx, t)) {
                        return 1;
                    }

                    if (avl_insert_or_inc(&((*cols)[myidx]), myidx, -t)) {
                        return 1;
                    }

                    n_rows = MAX(n_rows, myidx + 1);
                    n_cols = MAX(n_cols, idx + 1);
                }
            }
        }

        free(v);
    }

    coal_make_phdist(phdist, n_rows, n_cols, *rows, *cols, StSpM, ri, state_size);
}

int queue_init_standard_vector(struct queue_data **out, size_t state_size, vec_entry_t n) {
    *out = malloc(sizeof(struct queue_data));
    struct queue *queue = malloc(sizeof(struct queue));
    struct state_hobolth *state = malloc(sizeof(struct state_hobolth));

    expanding_arr_t *expanding_rows = malloc(sizeof(expanding_arr_t));
    expanding_arr_t *expanding_cols = malloc(sizeof(expanding_arr_t));

    avl_node_t** rows;
    avl_node_t** cols;

    rows = malloc(sizeof(avl_node_t*) * 1);
    cols = malloc(sizeof(avl_node_t*) * 1);

    expanding_rows->value = malloc(sizeof(avl_node_t***));
    *expanding_rows->value = rows;
    expanding_cols->value = malloc(sizeof(avl_node_t***));
    *expanding_cols->value = cols;
    expanding_rows->length = 1;
    expanding_cols->length = 1;
    expanding_rows->entry_size = sizeof(avl_node_t*);
    expanding_cols->entry_size = sizeof(avl_node_t*);

    rows[0] = malloc(sizeof(avl_node_t) * 1);
    avl_node_create(&(rows[0]), 0, 0, NULL);

    cols[0] = malloc(sizeof(avl_node_t) * 1);
    avl_node_create(&(cols[0]), 0, 0, NULL);

    vec_nmemb = state_size;

    queue_init(queue, 2);

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

    queue_enqueue(queue, BST);
    queue_enqueue(queue, BST_initial);
    memcpy(StSpM[0], mrca, state_size);
    memcpy(StSpM[1], initial, state_size);

    // Set n_rows and n_cols to 1, as the MRCA state does not
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

    (*out)->queue = queue;
    (*out)->vector_size = sizeof(vec_entry_t) * n;
    (*out)->state = state;
}

int queue_init_ss_hobolth(struct queue_data **out, size_t state_size) {
    queue_init_standard_vector(out, state_size, (vec_entry_t)state_size);
}

int coal_gen_phdist(phdist_t **phdist, size_t state_size) {
    struct queue_data *queue_data;
    queue_init_ss_hobolth(&queue_data, state_size);
    struct coal_args_hobolth_t args;
    args.n = state_size;

    queue_pop_ss_hobolth(phdist, queue_data, &args);

    return 0;
}

int coal_gen_phdist_reward(phdist_t **phdist, size_t n, size_t reward_index) {
    struct queue_data *queue_data;
    queue_init_standard_vector(&queue_data, n, (vec_entry_t)n);
    coal_args_compress_t args;
    args.n = n;
    args.reward_index = reward_index;

    queue_pop_ss_hobolth_compress(phdist, queue_data, &args);

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
