#include <stdio.h>
#include <time.h>
#include "../phdist.h"
#include "../coal.h"


void test_zeroes2() {
    phdist_t *phdist;
    time_t t;
    time(&t);
    coal_gen_phdist(&phdist, 60);
    printf("\n%zu", phdist_count_non_zeros(phdist));
    printf("\n time: %zu\n", time(NULL) - t);
}

void test_zeroes() {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 20);
    printf("\n%zu", phdist_count_non_zeros(phdist));
    mat_t *mat;

    time_t t;
    time(&t);

    mat_mult(&mat, phdist->si_mat, phdist->si_mat);
    phdist->si_mat = mat;

    printf(" %zu\n", phdist_count_non_zeros(phdist));
    printf("\n time: %zu\n", time(NULL) - t);
}

void test_gen() {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 5);
    phdist_print_as_matrix(phdist);
    printf("\n%zu\n", phdist_count_non_zeros(phdist));
    mat_t *mat;

    mat_mult(&mat, phdist->si_mat, phdist->si_mat);
    phdist->si_mat = mat;

    phdist_print_as_matrix(phdist);
    printf("\n%zu\n", phdist_count_non_zeros(phdist));
}

int main(int argc, char **argv) {
    test_gen();
    test_zeroes();
    test_zeroes2();
    return 0;
}