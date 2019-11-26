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


void test_scalar() {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 5);
    mat_t *unscaled = phdist->si_mat;
    double scalars[] = {0, 1, 2, 3, 4, 5, 6};
    mat_t *scaled;
    mat_scale_rows(&scaled, phdist->si_mat, scalars);
    phdist->si_mat = scaled;
    phdist_print_as_matrix(phdist);
    phdist_print_as_matrix_col(phdist);
    phdist->si_mat = unscaled;
    phdist_print_as_matrix(phdist);
}


void test_clone() {
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

void test_gen() {
    phdist_t *phdist2;
    coal_gen_phdist(&phdist2, 5);
    phdist_t *phdist;
    phdist_clone(&phdist, phdist2);
    phdist_print_as_matrix(phdist);
    phdist_print_as_matrix_col(phdist);
    printf("\n%zu\n", phdist_count_non_zeros(phdist));
}

void test_reward_sites() {
    phdist_t *phdist2;
    coal_gen_phdist(&phdist2, 5);
    phdist_t *phdist;
    phdist_print_as_matrix(phdist2);
    phdist_reward_transform(&phdist, phdist2);
    phdist_print_as_matrix(phdist);
    phdist_print_as_matrix_col(phdist);
}

int main(int argc, char **argv) {
    //test_gen();
    //test_clone();
    //test_zeroes();
    //test_zeroes2();
    //test_scalar();
    test_reward_sites();
    return 0;
}