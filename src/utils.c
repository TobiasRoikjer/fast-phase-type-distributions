#include "utils.h"
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

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

void vector_clear(vector_t *vector) {
    size_t entries = vector_length(vector);
    
    for (size_t i = 0; i < entries; ++i) {
        vector_remove_head(vector);
    }
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

/*
 * Adds the edge if it does not exists. Otherwise increments the weight
 */
int graph_combine_edge(graph_node_t *from, graph_node_t *to, weight_t weight) {
    weighted_edge_t *values = vector_get(from->edges);

    for (size_t i = 0; i < vector_length(from->edges); ++i) {
        if (values[i].node == to) {
            values[i].weight += weight;
            return 0;
        }
    }

    *((weighted_edge_t*)(vector_add(from->edges))) =
            (weighted_edge_t) {.node = to, .weight = weight};
    *((weighted_edge_t*)(vector_add(to->reverse_edges))) =
            (weighted_edge_t) {.node = from, .weight = 0};

    return 0;
}



int graph_remove_edge_forwards(graph_node_t *from, graph_node_t *to) {
    size_t length = vector_length(from->edges);
    bool found = false;
    weighted_edge_t *edges = (weighted_edge_t*)vector_get(from->edges);

    for (ssize_t i = length-1; i >= 0; i--) {
        weighted_edge_t edge = edges[i];

        if (edge.node == to) {
            if (found) {
                DIE_ERROR(1, "Edge already found\n");
            }

            vector_remove_entry(from->edges, (size_t) i);
            found = true;
        }
    }

    if (!found) {
        DIE_ERROR(1, "The edge was not found");
    }

    return 0;
}
int graph_remove_edge_backwards(graph_node_t *from, graph_node_t *to) {
    size_t length = vector_length(to->reverse_edges);
    bool found = false;

    for (ssize_t i = length-1; i >= 0; i--) {
        weighted_edge_t *edge = &(((weighted_edge_t*)vector_get(to->reverse_edges))[i]);

        if (edge->node == from) {
            if (found) {
                DIE_ERROR(1, "Edge already found\n");
            }

            vector_remove_entry(to->reverse_edges, (size_t) i);
            found = true;
        }
    }

    if (!found) {
        DIE_ERROR(1, "The edge was not found");
    }

    return 0;
}
int graph_remove_edge(graph_node_t *from, graph_node_t *to) {
   graph_remove_edge_forwards(from, to);
   graph_remove_edge_backwards(from, to);
   return 0;
}

int graph_redistribute_edge(graph_node_t *from, graph_node_t *to) {
    weight_t total_new_weight = 0;
    weight_t total_weight = 0;
    size_t length = vector_length(from->edges);
    bool found = false;

    for (size_t i = 0; i < length; i++) {
        weighted_edge_t *edge = &(((weighted_edge_t*)vector_get(from->edges))[i]);
        total_weight += edge->weight;

        if (edge->node == to) {
            if (found) {
                DIE_ERROR(1, "Edge already found\n");
            }

            vector_remove_entry(from->edges, i);
            found = true;
        } else {
            total_new_weight += edge->weight;
        }
    }

    if (!found) {
        DIE_ERROR(1, "The edge was not found");
    }

    for (size_t i = 0; i < vector_length(from->edges); i++) {
        weighted_edge_t *edge = &(((weighted_edge_t*)vector_get(from->edges))[i]);
        weight_t proportion = edge->weight / total_new_weight;
        edge->weight = proportion * total_weight;
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