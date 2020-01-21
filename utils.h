#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define DIE_PERROR(error_code, error, ...) do {     \
    char error_formatted[1024];                     \
    char error_formatted_line[1024];                \
                                                    \
    snprintf(error_formatted,                       \
            sizeof(error_formatted),                \
            error, ##__VA_ARGS__);                  \
    snprintf(error_formatted_line,                  \
            sizeof(error_formatted_line),           \
            "%s @ %s (%d)", error_formatted,        \
            __FILE__, __LINE__);                    \
                                                    \
    perror(error_formatted_line);                   \
    exit(error_code);                              \
} while(0)


typedef struct expanding_arr {
    void **value;
    size_t length;
    size_t entry_size;
} expanding_arr_t;

int expanding_arr_init(expanding_arr_t **arr, size_t entry_size, size_t initial_size);
void expanding_arr_destroy(expanding_arr_t *arr);
int expanding_arr_fit(expanding_arr_t *arr, size_t min_length);

typedef struct vector {
    size_t head_index;
    expanding_arr_t *expanding_arr;
} vector_t;

int vector_init(vector_t **vector, size_t entry_size, size_t initial_size);
void *vector_add(vector_t *vector);
void *vector_get(vector_t *vector);
size_t vector_size(vector_t *vector);

#endif
