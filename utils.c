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

/*
 * Note, this is a very expensive method of doing this.
 */
int vector_remove_entry(vector_t *vector, size_t index) {
    char **p = (char**)(vector->expanding_arr->value);

    for (size_t i = 0; i < vector->expanding_arr->entry_size; ++i) {
        *(*p + i + index * (vector->expanding_arr->entry_size)) =
                *(*p + i + (vector->head_index - 1) * (vector->expanding_arr->entry_size));
    }

    vector_remove_head(vector);

    return 0;
}

void vector_remove_head(vector_t *vector) {
    vector->head_index--;
}

inline void *vector_get(vector_t *vector) {
    return *vector->expanding_arr->value;
}

size_t vector_length(vector_t *vector) {
    return vector->head_index;
}

int graph_node_create(graph_node_t **node, size_t data_size) {
    *node = malloc(sizeof(graph_node_t) + data_size);
    vector_init(&(*node)->edges, sizeof(weighted_edge_t), 1);
    vector_init(&(*node)->reverse_edges, sizeof(weighted_edge_t), 1);
    return 0;
}

int graph_node_destroy(graph_node_t *node) {
    // TODO
}

int graph_add_edge(graph_node_t *from, graph_node_t *to, weight_t weight) {
    *((weighted_edge_t*)(vector_add(from->edges))) =
            (weighted_edge_t) {.node = to, .weight = weight};
    *((weighted_edge_t*)(vector_add(to->reverse_edges))) =
            (weighted_edge_t) {.node = from, .weight = weight};

    return 0;
}

int graph_redistribute_edge(graph_node_t *from, graph_node_t *to) {
    weight_t weight = 0;
    weight_t weight_increase = 0;

    for (size_t i = 0; i < vector_length(from->edges); i++) {
        weighted_edge_t *edge = &(((weighted_edge_t*)vector_get(from->edges))[i]);

        if (edge->node == to) {
            weight = edge->weight;
            vector_remove_entry(from->edges, i);
        }
    }

    if (weight == 0) {
        DIE_PERROR(1, "The weight was zero. Likely the edge was not found");
    }

    weight_increase = weight / vector_length(from->edges);

    for (size_t i = 0; i < vector_length(from->edges); i++) {
        weighted_edge_t *edge = &(((weighted_edge_t*)vector_get(from->edges))[i]);

        edge->weight += weight_increase;
    }

    return 0;
}

inline size_t graph_get_no_children(graph_node_t *node) {
    return vector_length(node->edges);
}

int queue_create(queue_t **queue, size_t nmemb) {
    *queue = malloc(sizeof(queue_t));
    
    if (((*queue)->queue = malloc(sizeof(void*)*nmemb)) == NULL) {
        return 1;
    }

    (*queue)->head = (*queue)->queue;
    (*queue)->tail = (*queue)->queue;
    (*queue)->size = nmemb;
    (*queue)->bound = (*queue)->queue + (*queue)->size;
    memset((*queue)->queue, 0, (*queue)->size *sizeof(void*));

    return 0;
}

int queue_destroy(queue_t *queue) {
    free(queue->queue);
    return 0;
}

/*void queue_print(queue_t *queue) {
    fprintf(stderr, "Queue: ");

    void **i = queue->queue;
    size_t j = 0;

    while (i < queue->queue + queue->size) {
        if (*i != NULL) {
            if (i == queue->head) {
                fprintf(stderr, "h%zu, ", (*i)->entry);
            } else if (i == queue->tail) {
                fprintf(stderr, "t%zu, ", (*i)->entry);
            } else {
                fprintf(stderr, "%zu, ", (*i)->entry);
            }
        } else {
            if (i == queue->head) {
                fprintf(stderr, "h_, ");
            } else if (i == queue->tail) {
                fprintf(stderr, "t_, ");
            } else {
                fprintf(stderr, "_, ");
            }
        }

        i++;
        j++;
    }

    fprintf(stderr, "\n");
}*/

int queue_enqueue(queue_t *queue, void *entry) {
    *(queue->tail) = entry;

    queue->tail++;

    if (queue->tail == queue->bound) {
        queue->tail = queue->queue;
    }

    if (queue->tail == queue->head) {
        size_t head_dist = queue->size - (queue->head - queue->queue);
        size_t head_offset = queue->head - queue->queue;
        size_t tail_offset = queue->tail - queue->queue;
        queue->queue = realloc(queue->queue, queue->size * 2 *sizeof(void**));

        memcpy(queue->queue + head_offset + queue->size,
                queue->queue + head_offset,
                head_dist * sizeof(void**));

        queue->tail = queue->queue + tail_offset;

        queue->head = queue->queue + queue->size + head_offset;
        queue->size *= 2;
        queue->bound = queue->queue + queue->size;
    }
    
    return 0;
}

void* queue_dequeue(queue_t *queue) {
    void *entry = *(queue->head);
    *(queue->head) = NULL;
    queue->head++;

    if (queue->head == queue->bound) {
        queue->head = queue->queue;
    }

    return entry;
}

int queue_empty(queue_t *queue) {
    return (queue->head == queue->tail);
}