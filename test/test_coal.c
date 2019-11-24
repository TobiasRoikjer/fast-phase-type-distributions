#include <stdio.h>

#include "../phdist.h"
#include "../coal.h"

int main(int argc, char **argv) {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 5);
    phdist_print_as_matrix(phdist);
    return 0;
}