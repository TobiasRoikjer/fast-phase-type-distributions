#ifndef AMAZEPHASE_BBST_H
#define AMAZEPHASE_BBST_H

#include <stdlib.h>

typedef size_t mat_key_t;
typedef double mat_entry_t;

typedef struct avl_node_t {
    mat_key_t key;
    struct avl_node_t *left;
    struct avl_node_t *right;
    struct avl_node_t *parent;
    mat_entry_t entry;
    signed short balance;
} avl_node_t;

typedef struct avl_flat_tuple {
    mat_key_t key;
    mat_entry_t entry;
} avl_flat_tuple_t;

int avl_node_create(avl_node_t **node, mat_key_t key, mat_entry_t entry, avl_node_t *parent);
void avl_node_destroy(avl_node_t *node);
int avl_insert_or_inc(avl_node_t **root, mat_key_t key, mat_entry_t entry);
int avl_flatten(avl_flat_tuple_t** arr, size_t *max_key, avl_node_t *root);
avl_node_t * avl_find(avl_node_t *rootptr, mat_key_t key);
void avl_print(avl_node_t *rootptr);

#endif // AMAZEPHASE_BSST_H