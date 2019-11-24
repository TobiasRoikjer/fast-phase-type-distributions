#include "bbst.h"

int avl_node_create(avl_node_t **node, mat_key_t key, mat_entry_t entry, avl_node_t *parent) {
    if ((*node = (avl_node_t*) malloc(sizeof(avl_node_t))) == NULL) {
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

void avl_node_destroy(avl_node_t *node) {
    if (node == NULL) {
        return;
    }

    avl_node_destroy(node->left);
    avl_node_destroy(node->right);

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

/*
 * Returns 1 on insert, 0 on find, -1 on error.
 */
int find_or_insert(avl_node_t **out, avl_node_t *rootptr, mat_key_t key, mat_entry_t entry) {
    if (avl_node_create(out, key, entry, NULL)) {
        return -1;
    }

    if (rootptr == NULL) {
        return 1;
    }

    avl_node_t *node = rootptr;

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

int avl_insert_or_inc(avl_node_t **root, mat_key_t key, mat_entry_t entry) {
    avl_node_t *child;

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

    avl_node_t *pivot, *rotated_parent;

    for (avl_node_t *parent = child->parent; parent != NULL; parent = child->parent) {
        if (child == parent->right) {
            if (parent->balance > 0) {
                pivot = parent->parent;

                if (child->balance < 0) {
                    rotated_parent = rotate_right_left(parent, child);
                } else {
                    rotated_parent = rotate_left(parent, child);
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
                    rotated_parent = rotate_left_right(parent, child);
                } else {
                    rotated_parent = rotate_right(parent, child);
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

static size_t avl_get_size(avl_node_t *node) {
    if (node == NULL) {
        return 0;
    }

    return 1 + avl_get_size(node->left) + avl_get_size(node->right);
}

avl_flat_tuple_t* avl_fill_arr(avl_node_t *node, avl_flat_tuple_t *arr) {
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

static size_t compare(const avl_flat_tuple_t *a, const avl_flat_tuple_t *b)
{
    return (a->key - b->key);
}

int avl_flatten(avl_flat_tuple_t** arr, size_t *n, avl_node_t *root) {
    *n = avl_get_size(root);
    if ((*arr = (avl_flat_tuple_t*) malloc(sizeof(avl_flat_tuple_t)*(*n + 1))) == NULL) {
        return 1;
    }

    avl_fill_arr(root, *arr);
    (*arr)[*n].entry = 0;
    qsort(*arr, *n, sizeof(avl_flat_tuple_t),
            (int(*)(const void *, const void*)) compare);

    return 0;
}