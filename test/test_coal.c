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
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 5);
    phdist_print_as_matrix(phdist);
    phdist_print_as_matrix_col(phdist);
    printf("\n%zu\n", phdist_count_non_zeros(phdist));
    phdist_t *phdist2;
    phdist_reward_transform(&phdist2, phdist);
    phdist_print_as_matrix(phdist2);
    phdist_print_as_matrix_col(phdist2);
}

void test_gen_erlang() {
    phdist_t *phdist;
    coal_gen_erlang_phdist(&phdist, 5);
    phdist_print_as_matrix(phdist);
    phdist_print_as_matrix_col(phdist);
    printf("\n%zu\n", phdist_count_non_zeros(phdist));
    phdist_t *phdist2;
    phdist_reward_transform(&phdist2, phdist);
    phdist_print_as_matrix(phdist2);
    phdist_print_as_matrix_col(phdist2);
}


void test_gen_reward() {
    coal_graph_node_t *graph;
    phdist_t *phdist;
    coal_gen_kingman_graph(&graph, 5);
    coal_graph_as_phdist_rw(&phdist, graph);
    phdist_print_as_matrix(phdist);
    phdist_print_as_matrix_col(phdist);
    printf("\n%zu\n", phdist_count_non_zeros(phdist));
}

void test_exp_reward() {
    coal_graph_node_t *graph;
    phdist_t *phdist;
    coal_gen_kingman_graph(&graph, 10);

    for (size_t i = 0; i < 10; i++) {
        printf("Expected %zu: %f\n", i, coal_mph_expected(graph, i));
    }
}


void test_cov_reward() {
    coal_graph_node_t *graph;
    phdist_t *phdist;
    coal_gen_kingman_graph(&graph, 30);

    for (size_t j = 0; j < 4; j++) {
        printf("Cov %zu: ", j);
        for (size_t i = 0; i < 4; i++) {
            printf("%f\t", coal_mph_cov(graph, i, j));
        }
        printf("\n");
    }
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

void test_mat_sub() {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 5);
    phdist->n_rw_cols = 0;
    phdist->n_rw_rows = 0;
    mat_t *id;
    mat_identity(&id, 7);
    mat_t *sub;
    mat_t *orig = phdist->si_mat;
    mat_sub(&sub, phdist->si_mat, id);
    phdist->si_mat = sub;
    phdist_print_as_matrix(phdist);
    phdist_print_as_matrix_col(phdist);

    mat_sub(&sub, id, orig);
    phdist->si_mat = sub;
    phdist_print_as_matrix(phdist);
    phdist_print_as_matrix_col(phdist);
}


void test_mat_scalar() {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 5);
    phdist->n_rw_cols = 0;
    phdist->n_rw_rows = 0;
    mat_t *scaled;
    mat_mul_scalar(&scaled, phdist->si_mat, 0.5);
    phdist->si_mat = scaled;
    phdist_print_as_matrix(phdist);
    phdist_print_as_matrix_col(phdist);
}




void test_mat_rowsum() {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 5);
    phdist->n_rw_cols = 0;
    phdist->n_rw_rows = 0;
    mat_t *sum;

    for (size_t r = 0; r < phdist->si_mat->n_rows; r++) {
        phdist->si_mat->rows[r]->entry = -3;
    }
    mat_row_sums(&sum, phdist->si_mat);
    phdist_print_as_matrix(phdist);
    phdist->si_mat = sum;
    phdist_print_as_matrix(phdist);
    phdist_print_as_matrix_col(phdist);
}



void test_inverse() {
    phdist_t phdist3;
    phdist3.n_rw_rows = 0;
    phdist3.n_rw_cols = 0;
    mat_t *id;
    mat_identity(&id, 5);
    phdist3.si_mat = id;
    phdist_print_as_matrix(&phdist3);
    phdist_print_as_matrix_col(&phdist3);

    phdist_t *phdist2;
    coal_gen_phdist(&phdist2, 5);

    phdist2->n_rw_rows = 0;
    phdist_print_as_matrix(phdist2);

    mat_t *m;
    mat_clone(&m, phdist2->si_mat);
    mat_t *inv;
    mat_inv_slow(&inv, m);
    phdist2->si_mat = inv;
    phdist_print_as_matrix(phdist2);
    phdist_print_as_matrix_col(phdist2);
    mat_t *s;
    mat_clone(&s, inv);
    mat_inv(&inv, m);
    phdist2->si_mat = inv;
    phdist_print_as_matrix(phdist2);
    phdist_print_as_matrix_col(phdist2);

    mat_t *res;
    mat_mult(&res, m, inv);
    phdist2->si_mat = res;
    phdist_print_as_matrix(phdist2);
}

void test_time_seg() {
    phdist_t *phdist;
    time_t t = time(NULL);
    coal_gen_phdist(&phdist, 10);
    d_dist_t *dist;
    coal_seg_sites(&dist, phdist);
    d_phgen_args_t *args = dist->args;
    args->theta = 4;
    double *out;
    dist->generator_fun(&out, 0, 20, args);

    for (size_t i = 0; i <= 20; i++) {
        printf("%f ", out[i]);
    }

    printf("\n Time: %zu\n", time(NULL) - t);
}


void test_time_seg_erlang() {
    phdist_t *phdist;
    time_t t = time(NULL);
    coal_gen_erlang_phdist(&phdist, 10);
    d_dist_t *dist;
    coal_seg_sites(&dist, phdist);
    d_phgen_args_t *args = dist->args;
    args->theta = 4;
    double *out;
    dist->generator_fun(&out, 0, 20, args);

    for (size_t i = 0; i <= 20; i++) {
        printf("%f ", out[i]);
    }

    printf("\n Time: %zu\n", time(NULL) - t);
}


void test_seg() {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 6);
    d_dist_t *dist;
    coal_seg_sites(&dist, phdist);
    d_phgen_args_t *args = dist->args;
    args->theta = 4;
    double *out;
    dist->generator_fun(&out, 0, 10, args);

    for (size_t i = 0; i <= 10; i++) {
        printf("%f ", out[i]);
    }

    printf("\n");
}


void test_time_inverse() {
    phdist_t *phdist2;
    coal_gen_phdist(&phdist2, 25);

    phdist2->n_rw_rows = 0;
    mat_t *inv;
    mat_mult(&(phdist2->si_mat), phdist2->si_mat, phdist2->si_mat);
    mat_mult(&(phdist2->si_mat), phdist2->si_mat, phdist2->si_mat);
    time_t t = time(NULL);
    mat_inv(&inv, phdist2->si_mat);
    phdist2->si_mat = inv;
    printf("Fast Time %zu\n", time(NULL) - t);


    t = time(NULL);
    mat_inv_slow(&inv, phdist2->si_mat);
    phdist2->si_mat = inv;
    printf("Slow Time %zu\n", time(NULL) - t);
}


void test_time_mul() {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 25);
    clock_t t = clock();
    mat_t *mat = phdist->si_mat;
    mat_inv(&mat, mat);
    //mat_print_as_matrix(mat);
    mat_mult(&mat, mat, mat);
    //mat_print_as_matrix(mat);
    printf("Time %zu\n", (clock() - t)/100000);
}


void test_mat_mul() {
    phdist_t *phdist;
    coal_gen_phdist(&phdist, 8);
    mat_t *mat = phdist->si_mat;
    mat_t *inv;
    mat_t *id;
    mat_identity(&id, mat->n_rows);
    mat_t *res;
    mat_mult(&res, id, mat);
    mat_print_as_matrix(res);
    mat_mult(&res, mat, id);
    mat_print_as_matrix(res);
    mat_inv(&inv, mat);
    mat_mult(&res, inv, mat);
    mat_print_as_matrix(res);

    mat_mult(&mat, mat, mat);
    mat_mult(&inv, inv, inv);
    mat_mult(&res, inv, mat);
    mat_print_as_matrix(res);

    mat_mult(&res, mat, mat);
    mat_mult(&res, inv, res);
    mat_mult(&res, res, inv);
    mat_print_as_matrix(res);
}

void test_im_mat_utils() {
    im_state_t *state;
    im_state_init(&state, 3, 4);
    state->mat1[1][0]= 3;
    state->mat1[0][1]= 4;
    state->mat1[3][0]= 1;
    state->mat1[0][4]= 2;
    state->mat1[3][4]= 9;

    state->mat2[1][0]= 3;
    state->mat2[0][1]= 4;
    state->mat2[3][0]= 1;
    state->mat2[0][4]= 2;
    state->mat2[0][2]= 3;
    state->mat2[1][2]= 4;
    state->mat2[2][2]= 5;
    state->mat2[3][4]= 8;

    vec_entry_t *vec;
    im_state_as_vec(&vec, state, 3, 4);
    fflush(stdout);

    for (size_t i = 0; i < (3+1)*(4+1)*2; ++i) {
        fprintf(stderr, "%zu,", vec[i]);
        if ((i+1) % 4 == 0) {
            fprintf(stderr, "\n");
        }
    }

    fprintf(stderr, "\n");
    fflush(stderr);
}

void test_gen_im() {
    coal_graph_node_t *graph;
    phdist_t *phdist;
    coal_gen_im_graph_args_t args = {
            .n1 = 1,
            .n2 = 2,
            .allow_back_migrations = false,
            .num_iso_coal_events = 1,
            .migration_param = 100.0f,
            .pop_scale1 = 10.0f,
            .pop_scale2 = 1000.0f,
            .mig_scale1 = 1.0f,
            .mig_scale2 = 1.0f
    };

    coal_gen_im_graph(&graph, args);
    coal_print_graph_list(stdout, graph, (args.n1 + 1)*(args.n2 + 1)*2+3,
                          (args.n1+1));
    coal_graph_as_phdist_rw(&phdist, graph);
    phdist_print_as_matrix(phdist);
    mat_print_as_matrix_with_abs(phdist->si_mat);
}

void test_gen_im_time() {
    for (size_t i = 0; i < 10; ++i) {
        coal_graph_node_t *graph;
        coal_gen_im_graph_args_t args = {
                .n1 = 10,
                .n2 = 10,
                .allow_back_migrations = false,
                .num_iso_coal_events = i,
                .migration_param = 100.0f,
                .pop_scale1 = 10.0f,
                .pop_scale2 = 1000.0f,
                .mig_scale1 = 1.0f,
                .mig_scale2 = 1.0f
        };

        coal_gen_im_graph(&graph, args);
    }
}

void test_gen_im_cutoff() {
    coal_graph_node_t *graph;
    phdist_t *phdist;
    coal_gen_im_cutoff_graph_args_t args = {
            .n1 = 5,
            .n2 = 0,
            .left_n1 = 2,
            .left_n2 = 0,
            .allow_back_migrations = true,
            .migration_param = 0.0f,
            .pop_scale1 = 1.0f,
            .pop_scale2 = 1.0f,
            .mig_scale1 = 1.0f,
            .mig_scale2 = 1.0f
    };

    coal_gen_im_cutoff_graph(&graph, args);
    coal_print_graph_list_im(stdout, graph, (args.n1 + 1)*(args.n2 + 1)*2+3,
                          (args.n1+1), args.n1, args.n2);
    weight_t **mat;
    size_t size;
    coal_graph_as_mat(&mat, &size, graph);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            fprintf(stdout, "%f ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }
}



int main(int argc, char **argv) {
    //test_gen();
    //printf("\n..\n");
    //test_gen_erlang();
    //printf("\n..\n");
    /*//test_clone();
    printf("\n..\n");
    //test_zeroes();
    printf("\n..\n");
    //test_zeroes2();
    printf("\n..\n");
    test_scalar();
    printf("\n..\n");
    test_reward_sites();
    printf("\n..\n");*/
    //test_inverse();
    //printf("\n..\n");
    //test_time_inverse();
    //printf("\n..\n");
    /*test_mat_sub();
    printf("\n..\n");
    test_mat_scalar();
    printf("\n..\n");
    test_mat_rowsum();
    printf("\n..\n");
    test_seg();
    printf("\n..\n");*/
    //test_time_seg();
    //printf("\n..\n");
    //test_time_seg_erlang();
    //printf("\n..\n");
    //test_time_mul();
    //printf("\n..\n");
    //test_mat_mul();
    //printf("\n..\n");
    //test_gen_reward();
    //printf("\n..\n");
    //test_exp_reward();
    //printf("\n..\n");
    //test_cov_reward();
    //printf("\n..\n");
    //test_im_mat_utils();
    //printf("\n..\n");
  //  test_gen_im();
    //printf("\n..\n");
    //test_gen_im_time();
    //printf("\n..\n");
    test_gen_im_cutoff();
    printf("\n..\n");

    return 0;
}