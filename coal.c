#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "coal.h"
#include "bbst.h"

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

typedef bst_node_t* queue_entry_t;

struct queue {
    queue_entry_t *queue;
    queue_entry_t *head;
    queue_entry_t *tail;
};

static int queue_init(struct queue *queue, size_t nmemb) {
    if ((queue->queue = malloc(sizeof(queue_entry_t)*nmemb)) == NULL) {
        return 1;
    }

    queue->head = queue->queue;
    queue->tail = queue->queue - 1;

    return 0;
}

static int queue_destroy(struct queue *queue) {
    free(queue->queue);
}

static void queue_enqueue(struct queue *queue, queue_entry_t entry) {
    queue->tail++;
    *(queue->tail) = entry;
}

static queue_entry_t queue_dequeue(struct queue *queue) {
    queue_entry_t entry = *(queue->head);
    queue->head++;

    return entry;
}

static int queue_empty(struct queue *queue) {
    return (queue->head > queue->tail);
}

#define MAX(a, b) (a>=b) ? a : b

int coal_gen_phdist(phdist_t **phdist, size_t state_size) {
    avl_node_t** rows;
    avl_node_t** cols;

    size_t avl_node_arr_len = 1;

    rows = malloc(sizeof(avl_node_t*));
    cols = malloc(sizeof(avl_node_t*));
    rows[0] = NULL;
    cols[0] = NULL;

    vec_entry_t n = (vec_entry_t)state_size;
 //   size_t nSt = partitions(state_size);

    size_t SIZE = sizeof(vec_entry_t) * n;
    vec_nmemb = state_size;

    struct queue queue;

    // TODO
    queue_init(&queue, 999999);

    vec_entry_t *initial = (vec_entry_t*)calloc(n, sizeof(vec_entry_t));
    initial[0] = n;

    size_t StSpM_n = 1;
    size_t **StSpM;

    if ((StSpM = malloc(sizeof(size_t*)*StSpM_n)) == NULL) {
        return 1;
    }

    /*if ((rows = calloc(2, sizeof(avl_node_t*))) == NULL) {
        return 1;
    }

    if ((cols = calloc(nSt, sizeof(avl_node_t*))) == NULL) {
        return 1;
    }*/

    StSpM[0] = malloc(sizeof(size_t) * state_size);

    /*for (size_t r = 0; r < 100; r++) {
        if ((StSpM[r] = malloc(sizeof(size_t) * state_size)) == NULL) {
            return 1;
        }
    }*/


    bst_node_t *BST = bst_init(initial, 0);

    queue_enqueue(&queue, BST);
    memcpy(StSpM[0], initial, state_size);
    size_t ri = 1;

    vec_entry_t *v;
    size_t myidx = 0;
    size_t idx;
    bst_node_t *idxN;
    bst_node_t *myidxN;
    size_t n_rows = 0;
    size_t n_cols = 0;

    while (queue_empty(&queue) == 0) {
        myidxN = queue_dequeue(&queue);

        myidx = myidxN->entry;
        v = (vec_entry_t*)malloc(SIZE);
        memcpy(v, myidxN->value, SIZE);

        for (vec_entry_t i = 0; i < n; i++) {
            for (vec_entry_t j = i; j < n; j++) {
                if ((i == j && v[i] >= 2) || (i != j && v[i] > 0 && v[j] > 0)) {
                    ssize_t t = i == j ? v[i]*(v[i]-1)/2 : v[i]*v[j];

                    v[i]--;
                    v[j]--;
                    v[(i + j + 2) - 1]++;

                    if (!bst_find_or_parent(&idxN, BST, v)) {
                        vec_entry_t *nv = (vec_entry_t*)malloc(SIZE);
                        memcpy(nv, v, SIZE);

                        for (vec_entry_t k = 0; k < n; k++) {
                            nv[k] = v[k];
                        }

                        idx = ri;

                        while (ri >= StSpM_n) {
                            // TODO: Failure
                            StSpM = realloc(StSpM, sizeof(size_t*)*StSpM_n*2);

                            //memset(StSpM+StSpM_n, 0, StSpM_n * sizeof(size_t*));

                            for (size_t i = StSpM_n; i < StSpM_n * 2; i++) {
                              StSpM[i] = malloc(sizeof(size_t) * state_size);
                            }

                            StSpM_n *= 2;
                        }

                        memcpy(StSpM[ri], nv, SIZE);
                        idxN = bst_add(idxN, nv, ri);
                        queue_enqueue(&queue, idxN);
                        ri = ri + 1;
                    } else {
                        idx = idxN->entry;
                    }

                    v[i]++;
                    v[j]++;
                    v[(i + j + 2) - 1]--;

                    while (myidx >= avl_node_arr_len || idx >= avl_node_arr_len) {
                        // TODO: Failure
                        rows = realloc(rows, sizeof(avl_node_t*)*avl_node_arr_len*2);
                        cols = realloc(cols, sizeof(avl_node_t*)*avl_node_arr_len*2);

                        memset(rows+avl_node_arr_len, 0, avl_node_arr_len * sizeof(avl_node_t*));
                        memset(cols+avl_node_arr_len, 0, avl_node_arr_len * sizeof(avl_node_t*));

                        avl_node_arr_len *= 2;
                    }

                    if (avl_insert_or_inc(&(rows[myidx]), idx, t)) {
                        return 1;
                    }

                    if (avl_insert_or_inc(&(rows[myidx]), myidx, -t)) {
                        return 1;
                    }


                    if (avl_insert_or_inc(&(cols[idx]), myidx, t)) {
                        return 1;
                    }

                    if (avl_insert_or_inc(&(cols[myidx]), myidx, -t)) {
                        return 1;
                    }

                    n_rows = MAX(n_rows, myidx+1);
                    n_cols = MAX(n_cols, idx+1);
                }
            }
        }

        free(v);
    }

    // TODO free some more

    //TODO: Add final row

    ret:
    *phdist = malloc(sizeof(phdist_t));
    size_t rn;
    (*phdist)->si_mat = malloc(sizeof(mat_t));
    (*phdist)->si_mat->rows = malloc(sizeof(avl_flat_tuple_t*)*n_rows);
    (*phdist)->si_mat->cols = malloc(sizeof(avl_flat_tuple_t*)*avl_node_arr_len);
    (*phdist)->si_mat->n_rows = n_rows;
    (*phdist)->si_mat->n_cols= n_cols;

    for (size_t i = 0; i < n_rows; i++) {
        avl_flatten(&((*phdist)->si_mat->rows[i]), &rn, rows[i]);
    }

    for (size_t i = 0; i < n_cols; i++) {
        avl_flatten(&((*phdist)->si_mat->cols[i]), &rn, cols[i]);
    }

    // TODO: Use normal matrix for this
    (*phdist)->rw_mat = malloc(sizeof(mat_t));

    return 0;
}