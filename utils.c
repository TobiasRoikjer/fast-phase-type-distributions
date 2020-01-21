#include "utils.h"
#include <string.h>
#include <stdint.h>

int expanding_arr_init(expanding_arr_t **arr, size_t entry_size, size_t initial_size) {
    *arr = malloc(sizeof(expanding_arr_t));
    (*arr)->length = initial_size;
    (*arr)->entry_size = entry_size;
    (*arr)->value = malloc(sizeof(void**));
    *(*arr)->value = calloc(initial_size, entry_size);

    return 0;
}

void expanding_arr_destroy(expanding_arr_t *arr) {
    free(*arr->value);
    free(arr->value);
    free(arr);
}

int expanding_arr_fit(expanding_arr_t *arr, size_t min_length) {
    if (min_length > 1000) {
        exit(1);
    }

    while (min_length >= arr->length - 1) {
        if (((*arr->value) = realloc((*arr->value), arr->entry_size * arr->length * 2)) == NULL) {
            return 1;
        }

        // Note: cast to char pointer in order to increment pointer
        memset((char *)(*arr->value) + (arr->entry_size * arr->length),
               0, arr->length * arr->entry_size);

        arr->length *= 2;
    }

    return 0;
}

int vector_init(vector_t **vector, size_t entry_size, size_t initial_size) {
    *vector = malloc(sizeof(vector_t));
    (*vector)->head_index = 0;
    expanding_arr_init(&(*vector)->expanding_arr, entry_size, initial_size);

    return 0;
}

void *vector_add(vector_t *vector) {
    expanding_arr_fit(vector->expanding_arr, vector->head_index+1);
    char **p = (char**)(vector->expanding_arr->value);

    return *p + vector->head_index++ * vector->expanding_arr->entry_size;
}

void *vector_get(vector_t *vector) {
    return *vector->expanding_arr->value;
}

size_t vector_size(vector_t *vector) {
    return vector->head_index;
}
