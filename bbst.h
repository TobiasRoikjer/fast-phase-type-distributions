#ifndef AMAZEPHASE_BBST_H
#define AMAZEPHASE_BBST_H

#include <stdlib.h>

typedef size_t mat_key_t;
typedef double mat_entry_t;

typedef struct avl_node {
    struct avl_node *left;
    struct avl_node *right;
    struct avl_node *parent;
    signed short balance;
} avl_node_t;

typedef struct avl_mat_node {
    struct avl_mat_node *left;
    struct avl_mat_node *right;
    struct avl_mat_node *parent;
    signed short balance;
    mat_key_t key;
    mat_entry_t entry;
} avl_mat_node_t;

typedef struct avl_flat_tuple {
    mat_key_t key;
    mat_entry_t entry;
} avl_flat_tuple_t;

int avl_node_create(avl_mat_node_t **node, mat_key_t key, mat_entry_t entry, avl_mat_node_t *parent);
void avl_node_destroy(avl_mat_node_t *node);
int avl_insert_or_inc(avl_mat_node_t **root, mat_key_t key, mat_entry_t entry);
int avl_flatten(avl_flat_tuple_t** arr, size_t *max_key, avl_mat_node_t *root);
avl_mat_node_t * avl_find(avl_mat_node_t *rootptr, mat_key_t key);
void avl_print(avl_mat_node_t *rootptr);

typedef size_t vec_entry_t;

typedef struct avl_vec_node {
    struct avl_vec_node *left;
    struct avl_vec_node *right;
    struct avl_vec_node *parent;
    signed short balance;
    vec_entry_t *key;
    void *entry;
} avl_vec_node_t;

int avl_vec_node_create(avl_vec_node_t **node, vec_entry_t *key, void *entry, avl_vec_node_t *parent);
void avl_vec_node_destroy(avl_vec_node_t *node);
int avl_vec_insert(avl_vec_node_t **root, vec_entry_t *key, void *entry, size_t vec_length);
avl_vec_node_t * avl_vec_find(avl_vec_node_t *rootptr, vec_entry_t *key, size_t vec_length);

#endif // AMAZEPHASE_BSST_H