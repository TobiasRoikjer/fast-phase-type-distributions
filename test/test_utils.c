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

    for (size_t i = 0; i < vector_size(vector); ++i) {
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

    for (size_t i = 0; i < vector_size(vector); ++i) {
        printf("%zu %zu, ", values[i].a, values[i].b);
    }

    printf("\n");
}


int main(int argc, char **argv) {
    test_vector();
    test_vector2();

    return 0;
}