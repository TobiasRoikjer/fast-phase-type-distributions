#include <stdio.h>

#include "../phdist.h"
#include "../coal.h"

void test_zeroes() {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 55);
    printf("\n%zu", phdist_count_non_zeros(phdist));
    mat_t *mat;

    mat_mult(&mat, phdist->si_mat, phdist->si_mat);
    phdist->si_mat = mat;

    printf(" %zu\n", phdist_count_non_zeros(phdist));
}

int main(int argc, char **argv) {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 5);
    phdist_print_as_matrix(phdist);
    printf("\n%zu\n", phdist_count_non_zeros(phdist));
    mat_t *mat;

    mat_mult(&mat, phdist->si_mat, phdist->si_mat);
    phdist->si_mat = mat;

    phdist_print_as_matrix(phdist);
    printf("\n%zu\n", phdist_count_non_zeros(phdist));

    test_zeroes();
    return 0;
}