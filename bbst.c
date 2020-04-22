#include "bbst.h"

#include <stdio.h>
#include <stdbool.h>
#include <string.h>

int avl_node_create(avl_mat_node_t **node, mat_key_t key, mat_entry_t entry, avl_mat_node_t *parent) {
    if ((*node = (avl_mat_node_t*) malloc(sizeof(avl_mat_node_t))) == NULL) {
        return 1;
    }

    (*node)->key = key;
    (*node)->entry = entry;
    (*node)->left = NULL;
    (*node)->right = NULL;
    (*node)->parent = parent;
    (*node)->balance = 0;

    return 0;
}

void avl_node_destroy(avl_mat_node_t *node) {
    if (node == NULL) {
        return;
    }

    avl_node_destroy((avl_mat_node_t *) node->left);
    avl_node_destroy((avl_mat_node_t *) node->right);

    free(node);
}

/* Example:
 *     A            A            A
 *   B   (left)    B  (right)   D
 * C       ->    D      ->    C   B
 *   D         C
 * In this case:
 *  C: child
 *  B: parent
 *  D: child_right
 */
avl_node_t *rotate_left_right(avl_node_t *parent, avl_node_t *child) {
    avl_node_t *child_right_left, *child_right_right;
    avl_node_t *child_right = child->right;
    child_right_left = child_right->left;
    child->right = child_right_left;

    if (child_right_left != NULL) {
        child_right_left->parent = child;
    }

    child_right->left = child;
    child->parent = child_right;
    child_right_right = child_right->right;
    parent->left = child_right_right;

    if (child_right_right != NULL) {
        child_right_right->parent = parent;
    }

    child_right->right = parent;
    parent->parent = child_right;

    if (child_right->balance > 0) {
        parent->balance = -1;
        child->balance = 0;
    } else if (child_right->balance == 0) {
        parent->balance = 0;
        child->balance = 0;
    } else {
        parent->balance = 0;
        child->balance = +1;
    }

    child_right->balance = 0;

    return child_right;
}

/* Example:
 *  A          A            A
 *   B  (right)  B   (left)   D
 *     C   ->      D    -> B    C
 *   D               C
 * In this case:
 *  C: child
 *  B: parent
 *  D: child_left
 */
avl_node_t *rotate_right_left(avl_node_t *parent, avl_node_t *child) {
    avl_node_t *child_left_right, *child_left_left;
    avl_node_t *child_left = child->left;

    child_left_right = child_left->right;

    child->left = child_left_right;

    if (child_left_right != NULL) {
        child_left_right->parent = child;
    }

    child_left->right = child;

    child->parent = child_left;
    child_left_left = child_left->left;
    parent->right = child_left_left;

    if (child_left_left != NULL) {
        child_left_left->parent = parent;
    }

    child_left->left = parent;
    parent->parent = child_left;

    if (child_left->balance > 0) {
        parent->balance = -1;
        child->balance = 0;
    } else if (child_left->balance == 0) {
        parent->balance = 0;
        child->balance = 0;
    } else {
        parent->balance = 0;
        child->balance = 1;
    }

    child_left->balance = 0;

    return child_left;
}

/*
 * Example:
 *  A              B
 *    B   (left) A   C
 *      C   ->
 */
avl_node_t *rotate_left(avl_node_t *parent, avl_node_t *child) {
    avl_node_t *child_left;

    child_left = child->left;
    parent->right = child_left;

    if (child_left != NULL) {
        child_left->parent = parent;
    }

    child->left = parent;
    parent->parent = child;

    if (child->balance == 0) {
        parent->balance = 1;
        child->balance = -1;
    } else {
        parent->balance = 0;
        child->balance = 0;
    }

    return child;
}

/*
 * Example:
 *      A            B
 *    B    (right) C   A
 *  C        ->
 */
avl_node_t *rotate_right(avl_node_t *parent, avl_node_t *child) {
    avl_node_t *child_right;

    child_right = child->right;
    parent->left = child_right;

    if (child_right != NULL) {
        child_right->parent = parent;
    }

    child->right = parent;
    parent->parent = child;

    if (child->balance == 0) {
        parent->balance = +1;
        child->balance = -1;
    } else {
        parent->balance = 0;
        child->balance = 0;
    }

    return child;
}

avl_mat_node_t * avl_find(avl_mat_node_t *rootptr, mat_key_t key) {
    if (rootptr == NULL) {
        return NULL;
    }

    avl_mat_node_t *node = rootptr;

    while (1) {
        if (key < node->key) {
            if (node->left == NULL) {
                return NULL;
            } else {
                node = node->left;
            }
        } else if (key > node->key) {
            if (node->right == NULL) {
                return NULL;
            } else {
                node = node->right;
            }
        } else {
            return node;
        }
    }
}

int find_or_insert(avl_mat_node_t **out, avl_mat_node_t *rootptr, mat_key_t key, mat_entry_t entry) {
    if (avl_node_create(out, key, entry, NULL)) {
        return -1;
    }

    if (rootptr == NULL) {
        return 1;
    }

    avl_mat_node_t *node = rootptr;

    while (1) {
        if (key < node->key) {
            if (node->left == NULL) {
                node->left = *out;
                break;
            } else {
                node = node->left;
            }
        } else if (key > node->key) {
            if (node->right == NULL) {
                node->right = *out;
                break;
            } else {
                node = node->right;
            }
        } else {
            *out = node;
            return 0;
        }
    }

    (*out)->parent = node;

    return 1;
}

int avl_insert_or_inc(avl_mat_node_t **root, mat_key_t key, mat_entry_t entry) {
    avl_mat_node_t *child;

    if (*root == NULL) {
        if (avl_node_create(root, key, entry, NULL)) {
            return 1;
        }

        return 0;
    }

    int res = find_or_insert(&child, *root, key, entry);

    if (res == -1) {
        return 1;
    }

    if (res == 0) {
        child->entry += entry;
        return 0;
    }

    avl_mat_node_t *pivot, *rotated_parent;

    for (avl_mat_node_t *parent = child->parent; parent != NULL; parent = child->parent) {
        if (child == parent->right) {
            if (parent->balance > 0) {
                pivot = parent->parent;

                if (child->balance < 0) {
                    rotated_parent = (avl_mat_node_t *) rotate_right_left((avl_node_t *) parent,
                                                                          (avl_node_t *) child);
                } else {
                    rotated_parent = (avl_mat_node_t *) rotate_left((avl_node_t *) parent,
                                                                    (avl_node_t *) child);
                }
            } else {
                if (parent->balance < 0) {
                    parent->balance = 0;

                    return 0;
                }

                parent->balance = 1;
                child = parent;

                continue;
            }
        } else {
            if (parent->balance < 0) {
                pivot = parent->parent;

                if (child->balance > 0) {
                    rotated_parent = (avl_mat_node_t *) rotate_left_right((avl_node_t *) parent, (avl_node_t *) child);
                } else {
                    rotated_parent = (avl_mat_node_t *) rotate_right((avl_node_t *) parent, (avl_node_t *) child);
                }
            } else {
                if (parent->balance > 0) {
                    parent->balance = 0;

                    return 0;
                }

                parent->balance = -1;
                child = parent;
                continue;
            }
        }

        rotated_parent->parent = pivot;

        if (pivot != NULL) {
            if (parent == pivot->left) {
                pivot->left = rotated_parent;
            } else {
                pivot->right = rotated_parent;
            }

            return 0;
        } else {
            *root = rotated_parent;
        }
    }

    return 0;
}

static size_t avl_get_size(avl_mat_node_t *node) {
    if (node == NULL) {
        return 0;
    }

    return 1 + avl_get_size(node->left) + avl_get_size(node->right);
}

avl_flat_tuple_t* avl_fill_arr(avl_mat_node_t *node, avl_flat_tuple_t *arr) {
    if (node == NULL) {
        return arr;
    }

    (*arr).key = node->key;
    (*arr).entry = node->entry;
    arr++;

    arr = avl_fill_arr(node->left, arr);
    arr = avl_fill_arr(node->right, arr);

    return arr;
}

static ssize_t compare(const avl_flat_tuple_t *a, const avl_flat_tuple_t *b)
{
    return (a->key - b->key);
}

int avl_flatten(avl_flat_tuple_t** arr, size_t *max_key, avl_mat_node_t *root) {
    size_t n = avl_get_size(root);

    if ((*arr = (avl_flat_tuple_t*) malloc(sizeof(avl_flat_tuple_t)*(n + 1))) == NULL) {
        return 1;
    }

    avl_fill_arr(root, *arr);
    (*arr)[n].entry = 0;
    qsort(*arr, n, sizeof(avl_flat_tuple_t),
            (int(*)(const void *, const void*)) compare);

    if (n > 0) {
        *max_key = (*arr)[n - 1].key;
    } else {
        *max_key = 0;
    }

    return 0;
}

void avl_print_impl(avl_mat_node_t *rootptr, size_t indent) {
    if (rootptr == NULL) {
        for (size_t i = 0; i < indent; i++) {
            fprintf(stderr, "\t");
        }
        fprintf(stderr, "NULL\n");
        return;
    }

    for (size_t i = 0; i < indent; i++) {
        fprintf(stderr, "\t");
    }

    fprintf(stderr, "%zu %f [\n", rootptr->key, rootptr->entry);

    for (size_t i = 0; i < indent; i++) {
        fprintf(stderr, "\t");
    }

    fprintf(stderr, "left: \n");
    avl_print_impl(rootptr->left, indent+1);

    for (size_t i = 0; i < indent; i++) {
        fprintf(stderr, "\t");
    }

    fprintf(stderr, "right: \n");
    avl_print_impl(rootptr->right, indent+1);

    for (size_t i = 0; i < indent; i++) {
        fprintf(stderr, "\t");
    }

    fprintf(stderr, "]");
}

void avl_print(avl_mat_node_t *rootptr) {
    avl_print_impl(rootptr, 0);
}

avl_flat_tuple_t *avl_bs_flat(avl_flat_tuple_t **values, mat_key_t key) {
    size_t low = 0;
    size_t high = 0;

}


static inline int radix_cmp(const vec_entry_t* a, const vec_entry_t* b,
        const size_t vec_length) {
    return (memcmp(a, b, sizeof(vec_entry_t) * vec_length));
}

int avl_vec_node_create(avl_vec_node_t **node, vec_entry_t *key, void *entry, avl_vec_node_t *parent) {
    if ((*node = (avl_vec_node_t*) malloc(sizeof(avl_vec_node_t))) == NULL) {
        return 1;
    }

    (*node)->key = key;
    (*node)->entry = entry;
    (*node)->left = NULL;
    (*node)->right = NULL;
    (*node)->parent = parent;
    (*node)->balance = 0;

    return 0;
}

void avl_vec_node_destroy(avl_vec_node_t *node) {
    if (node == NULL) {
        return;
    }

    avl_vec_node_destroy(node->left);
    avl_vec_node_destroy(node->right);

    free(node);
}

avl_vec_node_t * avl_vec_find(const avl_vec_node_t *rootptr, const vec_entry_t *key, const size_t vec_length) {
    if (rootptr == NULL) {
        return NULL;
    }

    avl_vec_node_t *node = (avl_vec_node_t *) rootptr;

    while (1) {
        int res = radix_cmp(key, node->key, vec_length);
        if (res < 0) {
            if (node->left == NULL) {
                return NULL;
            } else {
                node = node->left;
            }
        } else if (res > 0) {
            if (node->right == NULL) {
                return NULL;
            } else {
                node = node->right;
            }
        } else {
            return node;
        }
    }
}

int find_or_insert_vec(avl_vec_node_t **out, avl_vec_node_t *rootptr, vec_entry_t *key, void *entry, const size_t vec_length) {
    if (avl_vec_node_create(out, key, entry, NULL)) {
        return -1;
    }

    if (rootptr == NULL) {
        return 1;
    }

    avl_vec_node_t *node = rootptr;

    while(1) {
        int res = radix_cmp(key, node->key, vec_length);
        if (res < 0) {
            if (node->left == NULL) {
                node->left = *out;
                break;
            } else {
                node = node->left;
            }
        } else if (res > 0) {
            if (node->right == NULL) {
                node->right = *out;
                break;
            } else {
                node = node->right;
            }
        } else {
            *out = node;
        }
    }

    (*out)->parent = node;

    return 0;
}

int avl_vec_insert(avl_vec_node_t **root, vec_entry_t *key, void *entry, const size_t vec_length) {
    avl_vec_node_t *child;

    if (*root == NULL) {
        if (avl_vec_node_create(root, key, entry, NULL)) {
            return 1;
        }

        return 0;
    }

    int res = find_or_insert_vec(&child, *root, key, entry, vec_length);

    if (res == -1) {
        return 1;
    }

    if (res == 0) {
        return 0;
    }

    avl_vec_node_t *pivot, *rotated_parent;

    for (avl_vec_node_t *parent = child->parent; parent != NULL; parent = child->parent) {
        if (child == parent->right) {
            if (parent->balance > 0) {
                pivot = parent->parent;

                if (child->balance < 0) {
                    rotated_parent = (avl_vec_node_t *) rotate_right_left((avl_node_t *) parent, (avl_node_t *) child);
                } else {
                    rotated_parent = (avl_vec_node_t *) rotate_left((avl_node_t *) parent, (avl_node_t *) child);
                }
            } else {
                if (parent->balance < 0) {
                    parent->balance = 0;

                    return 0;
                }

                parent->balance = 1;
                child = parent;

                continue;
            }
        } else {
            if (parent->balance < 0) {
                pivot = parent->parent;

                if (child->balance > 0) {
                    rotated_parent = (avl_vec_node_t *) rotate_left_right((avl_node_t *) parent, (avl_node_t *) child);
                } else {
                    rotated_parent = (avl_vec_node_t *) rotate_right((avl_node_t *) parent, (avl_node_t *) child);
                }
            } else {
                if (parent->balance > 0) {
                    parent->balance = 0;

                    return 0;
                }

                parent->balance = -1;
                child = parent;
                continue;
            }
        }

        rotated_parent->parent = pivot;

        if (pivot != NULL) {
            if (parent == pivot->left) {
                pivot->left = rotated_parent;
            } else {
                pivot->right = rotated_parent;
            }

            return 0;
        } else {
            *root = rotated_parent;
        }
    }

    return 0;
}

static size_t avl_vec_get_size(avl_vec_node_t *node) {
    if (node == NULL) {
        return 0;
    }

    return 1 + avl_vec_get_size(node->left) + avl_vec_get_size(node->right);
}