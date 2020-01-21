#include "utils.h"
#include <string.h>

int expanding_arr_fit(expanding_arr_t *arr, size_t min_length) {
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