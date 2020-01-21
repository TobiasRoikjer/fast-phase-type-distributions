#include <stdio.h>

#include "../bbst.h"

static void testBBSTavl_insert_or_inc1() {
    avl_node_t * root = NULL;

    avl_insert_or_inc(&root, 'M', 0);
    avl_insert_or_inc(&root, 'N', 0);
    avl_insert_or_inc(&root, 'O', 0);
    avl_insert_or_inc(&root, 'L', 0);
    avl_insert_or_inc(&root, 'K', 0);
    avl_insert_or_inc(&root, 'P', 0);
    avl_insert_or_inc(&root, 'Q', 0);
    avl_insert_or_inc(&root, 'H', 0);
    avl_insert_or_inc(&root, 'I', 0);
    avl_insert_or_inc(&root, 'A', 0);

    avl_print(root);
}

static void testBBSTavl_insert_or_inc2() {
    avl_node_t * root = NULL;

    avl_insert_or_inc(&root, 'J', 0);
    avl_insert_or_inc(&root, 'F', 0);
    avl_insert_or_inc(&root, 'D', 0);
    avl_insert_or_inc(&root, 'C', 0);
    avl_insert_or_inc(&root, 'G', 0);
    avl_insert_or_inc(&root, 'P', 0);
    avl_insert_or_inc(&root, 'L', 0);
    avl_insert_or_inc(&root, 'N', 0);
    avl_insert_or_inc(&root, 'V', 0);
    avl_insert_or_inc(&root, 'S', 0);
    avl_insert_or_inc(&root, 'Q', 0);
    avl_insert_or_inc(&root, 'U', 0);
    avl_insert_or_inc(&root, 'X', 0);

    avl_print(root);
}

static void testBBSTavl_flatten() {
    avl_node_t * root = NULL;

    avl_insert_or_inc(&root, 'J', 1);
    avl_insert_or_inc(&root, 'F', 2);
    avl_insert_or_inc(&root, 'D', 3);
    avl_insert_or_inc(&root, 'C', 4);
    avl_insert_or_inc(&root, 'G', 5);
    avl_insert_or_inc(&root, 'P', 6);
    avl_insert_or_inc(&root, 'L', 7);
    avl_insert_or_inc(&root, 'N', 8);
    avl_insert_or_inc(&root, 'V', 9);
    avl_insert_or_inc(&root, 'S', 11110);
    avl_insert_or_inc(&root, 'Q', 0);
    avl_insert_or_inc(&root, 'U', 0);
    avl_insert_or_inc(&root, 'X', 90);

    avl_flat_tuple_t *arr;
    size_t n;

    avl_flatten(&arr, &n, root);

    for (size_t i = 0; i <= n; i++) {
        printf("%c %u, ", arr[i].key, arr[i].entry);
    }
}

static void testBBSTinc() {
    avl_node_t * root = NULL;

    avl_insert_or_inc(&root, 'M', 1);
    avl_insert_or_inc(&root, 'M', 10);
    avl_insert_or_inc(&root, 'M', 100);
    avl_insert_or_inc(&root, 'A', 5);
    avl_insert_or_inc(&root, 'A', 100);
    avl_insert_or_inc(&root, 'M', 1000);

    avl_print(root);
}

int main(int argc, char **argv) {
    testBBSTavl_insert_or_inc1();
    testBBSTavl_insert_or_inc2();
    testBBSTinc();
    testBBSTavl_flatten();

    return 0;
}