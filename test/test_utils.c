#include <stdio.h>
#include <inttypes.h>

#include "../utils.h"

static void test_vector() {
    vector_t *vector;
    vector_init(&vector, sizeof(uint32_t), 4);

    uint32_t i = 7;

    *(uint32_t*)vector_add(vector) = i;
    i++;
    *(uint32_t*)vector_add(vector) = i;
    i++;
    *(uint32_t*)vector_add(vector) = i;
    i++;
    *(uint32_t*)vector_add(vector) = i;
    i++;
    *(uint32_t*)vector_add(vector) = i;
    i++;
    *(uint32_t*)vector_add(vector) = i;

    uint32_t *values = vector_get(vector);

    for (size_t i = 0; i < vector_length(vector); ++i) {
        printf("%"PRIu32",", values[i]);
    }

    printf("\n");
}

struct test {
    size_t a;
    size_t b;
};

static void test_vector2() {
    vector_t *vector;
    vector_init(&vector, sizeof(struct test), 4);

    struct test i = (struct test){.a = 1, .b = 2};

    *((struct test*)vector_add(vector)) = i;
    i = (struct test){.a = 9, .b = 4};
    *((struct test*)vector_add(vector)) = i;
    i.a = 44;
    *((struct test*)vector_add(vector)) = i;

    struct test *values = vector_get(vector);

    for (size_t i = 0; i < vector_length(vector); ++i) {
        printf("%zu %zu, ", values[i].a, values[i].b);
    }

    printf("\n");
}

static void test_vector3() {
    vector_t *vector;
    vector_init(&vector, sizeof(uint32_t), 4);

    uint32_t i = 7;

    *(uint32_t*)vector_add(vector) = i;
    i++;
    *(uint32_t*)vector_add(vector) = i;
    i++;
    *(uint32_t*)vector_add(vector) = i;
    vector_remove_entry(vector, 2);
    vector_remove_entry(vector, 1);
    // 7 8 9 -> 7
    i++;
    *(uint32_t*)vector_add(vector) = i;
    i++;
    *(uint32_t*)vector_add(vector) = i;
    i++;
    *(uint32_t*)vector_add(vector) = i;
    // 7 10 11 12 -> 7 10 11
    vector_remove_entry(vector, 3);
    // 7 10 11 12 -> 10 11
    vector_remove_entry(vector, 0);

    uint32_t *values = vector_get(vector);

    for (size_t i = 0; i < vector_length(vector); ++i) {
        printf("%"PRIu32",", values[i]);
    }

    printf("\n");
}

struct data {
    size_t foo;
    char t;
};

static void print_graph_node(graph_node_t *node) {
    weighted_edge_t *values = vector_get(node->edges);

    printf("Node data %zu %c\n", ((struct data*)(&node->data))->foo,
           ((struct data*)(&node->data))->t);

    for (size_t i = 0; i < vector_length(node->edges); i++) {
        printf("Edge weight %f node data %zu %c\n", values[i].weight,
                ((struct data*)(&values[i].node->data))->foo,
               ((struct data*)(&values[i].node->data))->t);
    }

    printf("\n");
}

static void test_graph() {
    graph_node_t *A;
    graph_node_t *B;
    graph_node_t *C;
    graph_node_t *D;


    graph_node_create(&A, sizeof(struct data));
    graph_node_create(&B, sizeof(struct data));
    graph_node_create(&C, sizeof(struct data));
    graph_node_create(&D, sizeof(struct data));
    ((struct data*)(&A->data))->foo = 1;
    ((struct data*)(&A->data))->t = 'A';
    ((struct data*)(&B->data))->foo = 2;
    ((struct data*)(&B->data))->t = 'B';
    ((struct data*)(&C->data))->foo = 3;
    ((struct data*)(&C->data))->t = 'C';
    ((struct data*)(&D->data))->foo = 4;
    ((struct data*)(&D->data))->t = 'D';

    graph_add_edge(A, B, 1);
    graph_add_edge(A, C, 2);
    graph_add_edge(C, B, 3);
    graph_add_edge(B, D, 4);

    print_graph_node(A);
    print_graph_node(B);
    print_graph_node(C);
    print_graph_node(D);
}

static void test_graph2() {
    graph_node_t *A;
    graph_node_t *B;
    graph_node_t *C;
    graph_node_t *D;


    graph_node_create(&A, sizeof(struct data));
    graph_node_create(&B, sizeof(struct data));
    graph_node_create(&C, sizeof(struct data));
    graph_node_create(&D, sizeof(struct data));
    ((struct data*)(&A->data))->foo = 1;
    ((struct data*)(&A->data))->t = 'A';
    ((struct data*)(&B->data))->foo = 2;
    ((struct data*)(&B->data))->t = 'B';
    ((struct data*)(&C->data))->foo = 3;
    ((struct data*)(&C->data))->t = 'C';
    ((struct data*)(&D->data))->foo = 4;
    ((struct data*)(&D->data))->t = 'D';

    graph_add_edge(A, B, 1);
    graph_add_edge(A, C, 2);
    graph_add_edge(C, B, 3);
    graph_add_edge(B, D, 4);

    graph_redistribute_edge(A, B);
    graph_redistribute_edge(C, B);

    print_graph_node(A);
    print_graph_node(B);
    print_graph_node(C);
    print_graph_node(D);
}

static void test_queue() {
    queue_t *queue;
    queue_create(&queue, 1);
    char *A = malloc(sizeof(char));
    *A = 'A';
    char *B = malloc(sizeof(char));;
    *B = 'B';
    char *C = malloc(sizeof(char));;
    *C = 'C';
    char *D = malloc(sizeof(char));;
    *D = 'D';
    
    queue_enqueue(queue, A);
    queue_enqueue(queue, B);
    queue_enqueue(queue, C);
    printf("%c" , *(char*)queue_dequeue(queue));
    printf("%c" , *(char*)queue_dequeue(queue));
    queue_enqueue(queue, D);
    printf("%c" , *(char*)queue_dequeue(queue));
    printf("%c" , *(char*)queue_dequeue(queue));
}

int main(int argc, char **argv) {
    //test_vector();
    //test_vector2();
    //test_vector3();
    //test_graph();
    test_graph2();
    //test_queue();

    return 0;
}