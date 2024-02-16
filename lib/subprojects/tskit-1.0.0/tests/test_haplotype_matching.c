/*
 * MIT License
 *
 * Copyright (c) 2019-2022 Tskit Developers
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "testlib.h"
#include <tskit/haplotype_matching.h>

#include <unistd.h>
#include <stdlib.h>

/****************************************************************
 * TestHMM
 ****************************************************************/

static double
tsk_ls_hmm_compute_normalisation_factor_site_test(tsk_ls_hmm_t *TSK_UNUSED(self))
{
    return 1.0;
}

static int
tsk_ls_hmm_next_probability_test(tsk_ls_hmm_t *TSK_UNUSED(self),
    tsk_id_t TSK_UNUSED(site_id), double TSK_UNUSED(p_last), bool TSK_UNUSED(is_match),
    tsk_id_t TSK_UNUSED(node), double *result)
{
    *result = rand();
    /* printf("next proba = %f\n", *result); */
    return 0;
}

static int
run_test_hmm(tsk_ls_hmm_t *hmm, int32_t *haplotype, tsk_compressed_matrix_t *output)
{
    int ret = 0;

    srand(1);

    ret = tsk_ls_hmm_run(hmm, haplotype, tsk_ls_hmm_next_probability_test,
        tsk_ls_hmm_compute_normalisation_factor_site_test, output);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

/****************************************************************
 * TestHMM
 ****************************************************************/

static void
test_single_tree_missing_alleles(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_ls_hmm_t ls_hmm;
    tsk_compressed_matrix_t forward;
    tsk_viterbi_matrix_t viterbi;

    double rho[] = { 0, 0.25, 0.25 };
    double mu[] = { 0.125, 0.125, 0.125 };
    int32_t h[] = { 0, 0, 0, 0 };

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);

    ret = tsk_ls_hmm_init(&ls_hmm, &ts, rho, mu, TSK_ALLELES_ACGT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_ls_hmm_forward(&ls_hmm, h, &forward, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ALLELE_NOT_FOUND);
    ret = tsk_ls_hmm_viterbi(&ls_hmm, h, &viterbi, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_ALLELE_NOT_FOUND);

    tsk_ls_hmm_free(&ls_hmm);
    tsk_compressed_matrix_free(&forward);
    tsk_viterbi_matrix_free(&viterbi);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_exact_match(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_ls_hmm_t ls_hmm;
    tsk_compressed_matrix_t forward;
    tsk_viterbi_matrix_t viterbi;

    double rho[] = { 0.0, 0.25, 0.25 };
    double mu[] = { 0, 0, 0 };
    int32_t h[] = { 1, 1, 1 };
    tsk_id_t path[3];
    double decoded_compressed_matrix[12];
    unsigned int precision;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);

    ret = tsk_ls_hmm_init(&ls_hmm, &ts, rho, mu, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_ls_hmm_forward(&ls_hmm, h, &forward, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_compressed_matrix_print_state(&forward, _devnull);
    tsk_ls_hmm_print_state(&ls_hmm, _devnull);
    ret = tsk_compressed_matrix_decode(&forward, decoded_compressed_matrix);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_ls_hmm_viterbi(&ls_hmm, h, &viterbi, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_viterbi_matrix_print_state(&viterbi, _devnull);
    tsk_ls_hmm_print_state(&ls_hmm, _devnull);
    ret = tsk_viterbi_matrix_traceback(&viterbi, path, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(path[0], 2);
    CU_ASSERT_EQUAL(path[1], 1);
    CU_ASSERT_EQUAL(path[2], 1);

    /* Should get the same answer at lower precision */
    for (precision = 1; precision < 24; precision++) {
        ret = tsk_ls_hmm_set_precision(&ls_hmm, precision);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_ls_hmm_viterbi(&ls_hmm, h, &viterbi, TSK_NO_INIT);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tsk_viterbi_matrix_print_state(&viterbi, _devnull);
        tsk_ls_hmm_print_state(&ls_hmm, _devnull);
        ret = tsk_viterbi_matrix_traceback(&viterbi, path, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(path[0], 2);
        CU_ASSERT_EQUAL(path[1], 1);
        CU_ASSERT_EQUAL(path[2], 1);
    }

    tsk_ls_hmm_free(&ls_hmm);
    tsk_compressed_matrix_free(&forward);
    tsk_viterbi_matrix_free(&viterbi);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_missing_haplotype_data(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_ls_hmm_t ls_hmm;
    tsk_compressed_matrix_t forward;
    tsk_viterbi_matrix_t viterbi;

    double rho[] = { 0.0, 0.25, 0.25 };
    double mu[] = { 0, 0, 0 };
    int32_t h[] = { 1, TSK_MISSING_DATA, 1 };
    tsk_id_t path[3];
    double decoded_compressed_matrix[12];

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);

    ret = tsk_ls_hmm_init(&ls_hmm, &ts, rho, mu, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_ls_hmm_forward(&ls_hmm, h, &forward, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_compressed_matrix_print_state(&forward, _devnull);
    tsk_ls_hmm_print_state(&ls_hmm, _devnull);
    ret = tsk_compressed_matrix_decode(&forward, decoded_compressed_matrix);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_ls_hmm_viterbi(&ls_hmm, h, &viterbi, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_viterbi_matrix_print_state(&viterbi, _devnull);
    tsk_ls_hmm_print_state(&ls_hmm, _devnull);
    ret = tsk_viterbi_matrix_traceback(&viterbi, path, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(path[0], 2);
    CU_ASSERT_EQUAL(path[1], 2);
    CU_ASSERT_EQUAL(path[2], 2);

    tsk_ls_hmm_free(&ls_hmm);
    tsk_compressed_matrix_free(&forward);
    tsk_viterbi_matrix_free(&viterbi);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_match_impossible(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_ls_hmm_t ls_hmm;
    tsk_compressed_matrix_t forward;
    tsk_viterbi_matrix_t viterbi;

    double rho[] = { 0.0, 0.25, 0.25 };
    double mu[] = { 0, 0, 0 };
    /* This haplotype can't happen with a mutation rate of 0 */
    int32_t h[] = { 0, 0, 0 };

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);

    ret = tsk_ls_hmm_init(&ls_hmm, &ts, rho, mu, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_ls_hmm_forward(&ls_hmm, h, &forward, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MATCH_IMPOSSIBLE);
    tsk_compressed_matrix_print_state(&forward, _devnull);
    tsk_ls_hmm_print_state(&ls_hmm, _devnull);

    ret = tsk_ls_hmm_viterbi(&ls_hmm, h, &viterbi, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MATCH_IMPOSSIBLE);
    tsk_viterbi_matrix_print_state(&viterbi, _devnull);
    tsk_ls_hmm_print_state(&ls_hmm, _devnull);

    tsk_ls_hmm_free(&ls_hmm);
    tsk_compressed_matrix_free(&forward);
    tsk_viterbi_matrix_free(&viterbi);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_errors(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_ls_hmm_t ls_hmm;
    tsk_compressed_matrix_t forward;
    tsk_viterbi_matrix_t viterbi;
    tsk_value_transition_t T[1];
    double decoded[3][4];

    double rho[] = { 0.0, 0.25, 0.25 };
    double mu[] = { 0, 0, 0 };
    int32_t h[] = { 0, 0, 0 };

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);

    ret = tsk_viterbi_matrix_init(&viterbi, &ts, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_compressed_matrix_init(&forward, &ts, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_ls_hmm_init(&ls_hmm, &ts, rho, mu, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    viterbi.matrix.tree_sequence = NULL;
    ret = tsk_ls_hmm_viterbi(&ls_hmm, h, &viterbi, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    viterbi.matrix.tree_sequence = &ts;

    forward.tree_sequence = NULL;
    ret = tsk_ls_hmm_forward(&ls_hmm, h, &forward, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    forward.tree_sequence = &ts;

    ret = tsk_compressed_matrix_store_site(&forward, 3, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    ret = tsk_compressed_matrix_store_site(&forward, 4, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);

    T[0].tree_node = -1;
    T[0].value = 0;
    ret = tsk_compressed_matrix_store_site(&forward, 0, 1, 1, T);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_compressed_matrix_decode(&forward, (double *) decoded);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    T[0].tree_node = 7;
    T[0].value = 0;
    ret = tsk_compressed_matrix_store_site(&forward, 0, 1, 1, T);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_compressed_matrix_decode(&forward, (double *) decoded);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    tsk_ls_hmm_free(&ls_hmm);
    tsk_compressed_matrix_free(&forward);
    tsk_viterbi_matrix_free(&viterbi);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_compressed_matrix(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_compressed_matrix_t matrix;
    tsk_ls_hmm_t ls_hmm;
    tsk_size_t max_transitions = 1024;
    tsk_value_transition_t T[max_transitions];
    double decoded[3][4];
    int j;

    double rho[] = { 0.0, 0.25, 0.25 };
    double mu[] = { 0.1, 0.1, 0.1 };
    int32_t h[] = { 0, 0, 0 };

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);

    ret = tsk_compressed_matrix_init(&matrix, &ts, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_compressed_matrix_print_state(&matrix, _devnull);

    T[0].tree_node = 6;
    T[0].value = 0;
    for (j = 0; j < 3; j++) {
        T[1].tree_node = j;
        T[1].value = 1;
        ret = tsk_compressed_matrix_store_site(&matrix, j, 1.0, 2, T);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
    tsk_compressed_matrix_print_state(&matrix, _devnull);

    ret = tsk_compressed_matrix_decode(&matrix, (double *) decoded);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(decoded[0][0], 1.0);
    CU_ASSERT_EQUAL(decoded[0][1], 0.0);
    CU_ASSERT_EQUAL(decoded[0][2], 0.0);
    CU_ASSERT_EQUAL(decoded[1][0], 0.0);
    CU_ASSERT_EQUAL(decoded[1][1], 1.0);
    CU_ASSERT_EQUAL(decoded[1][2], 0.0);
    CU_ASSERT_EQUAL(decoded[2][0], 0.0);
    CU_ASSERT_EQUAL(decoded[2][1], 0.0);
    CU_ASSERT_EQUAL(decoded[2][2], 1.0);

    /* Cleared matrix should be zero everywhere */
    tsk_compressed_matrix_clear(&matrix);
    ret = tsk_compressed_matrix_decode(&matrix, (double *) decoded);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < 3; j++) {
        CU_ASSERT_EQUAL(decoded[j][0], 0.0);
        CU_ASSERT_EQUAL(decoded[j][1], 0.0);
        CU_ASSERT_EQUAL(decoded[j][2], 0.0);
    }

    ret = tsk_ls_hmm_init(&ls_hmm, &ts, rho, mu, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_ls_hmm_forward(&ls_hmm, h, &matrix, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_compressed_matrix_print_state(&matrix, _devnull);
    ret = tsk_compressed_matrix_decode(&matrix, (double *) decoded);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_compressed_matrix_free(&matrix);
    tsk_ls_hmm_free(&ls_hmm);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_viterbi_matrix(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_viterbi_matrix_t viterbi;
    tsk_ls_hmm_t ls_hmm;
    double rho[] = { 0.0, 0.25, 0.25 };
    double mu[] = { 0, 0, 0 };
    int32_t h[] = { 1, 1, 1 };
    tsk_id_t path[3];
    tsk_value_transition_t T[2];
    int j;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);

    ret = tsk_viterbi_matrix_init(&viterbi, &ts, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_viterbi_matrix_print_state(&viterbi, _devnull);
    ret = tsk_viterbi_matrix_traceback(&viterbi, path, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NULL_VITERBI_MATRIX);

    T[0].tree_node = 6;
    T[0].value = 0;
    T[1].tree_node = 1;
    T[1].value = 1;
    for (j = 0; j < 3; j++) {
        ret = tsk_compressed_matrix_store_site(&viterbi.matrix, j, 1.0, 2, T);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        /* We need to have one record per site, so we put in a record
         * at the root saying we don't need to recombine */
        ret = tsk_viterbi_matrix_add_recombination_required(&viterbi, j, 6, false);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
    ret = tsk_viterbi_matrix_traceback(&viterbi, path, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(path[0], 1);
    CU_ASSERT_EQUAL_FATAL(path[1], 1);
    CU_ASSERT_EQUAL_FATAL(path[2], 1);

    ret = tsk_ls_hmm_init(&ls_hmm, &ts, rho, mu, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_ls_hmm_viterbi(&ls_hmm, h, &viterbi, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_viterbi_matrix_print_state(&viterbi, _devnull);
    ret = tsk_viterbi_matrix_traceback(&viterbi, path, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_viterbi_matrix_clear(&viterbi);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_viterbi_matrix_traceback(&viterbi, path, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NULL_VITERBI_MATRIX);

    tsk_viterbi_matrix_free(&viterbi);

    ret = tsk_viterbi_matrix_init(&viterbi, &ts, 1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Make sure we hit the realloc case for recombination records */
    for (j = 0; j < 100; j++) {
        ret = tsk_viterbi_matrix_add_recombination_required(&viterbi, 0, 6, false);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
    tsk_viterbi_matrix_print_state(&viterbi, _devnull);

    tsk_viterbi_matrix_free(&viterbi);
    tsk_ls_hmm_free(&ls_hmm);
    tsk_treeseq_free(&ts);
}

static void
test_multi_tree_exact_match(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_ls_hmm_t ls_hmm;
    tsk_compressed_matrix_t forward;
    tsk_viterbi_matrix_t viterbi;

    double rho[] = { 0.0, 0.25, 0.25 };
    double mu[] = { 0, 0, 0 };
    int32_t h[] = { 1, 1, 1 };
    tsk_id_t path[3];
    double decoded_compressed_matrix[12];
    unsigned int precision;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, paper_ex_sites,
        paper_ex_mutations, paper_ex_individuals, NULL, 0);

    ret = tsk_ls_hmm_init(&ls_hmm, &ts, rho, mu, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_ls_hmm_forward(&ls_hmm, h, &forward, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_ls_hmm_print_state(&ls_hmm, _devnull);
    tsk_compressed_matrix_print_state(&forward, _devnull);
    ret = tsk_compressed_matrix_decode(&forward, decoded_compressed_matrix);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_ls_hmm_viterbi(&ls_hmm, h, &viterbi, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_viterbi_matrix_print_state(&viterbi, _devnull);
    tsk_ls_hmm_print_state(&ls_hmm, _devnull);
    ret = tsk_viterbi_matrix_traceback(&viterbi, path, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(path[0], 2);
    CU_ASSERT_EQUAL(path[1], 0);
    CU_ASSERT_EQUAL(path[2], 1);

    /* Should get the same answer at lower precision */
    for (precision = 4; precision < 24; precision++) {
        ret = tsk_ls_hmm_set_precision(&ls_hmm, precision);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_ls_hmm_viterbi(&ls_hmm, h, &viterbi, TSK_NO_INIT);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tsk_viterbi_matrix_print_state(&viterbi, _devnull);
        tsk_ls_hmm_print_state(&ls_hmm, _devnull);
        ret = tsk_viterbi_matrix_traceback(&viterbi, path, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(path[0], 2);
        CU_ASSERT_EQUAL(path[1], 0);
        CU_ASSERT_EQUAL(path[2], 1);
    }

    tsk_ls_hmm_free(&ls_hmm);
    tsk_compressed_matrix_free(&forward);
    tsk_viterbi_matrix_free(&viterbi);
    tsk_treeseq_free(&ts);
}

static void
test_multi_tree_errors(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_compressed_matrix_t forward;
    tsk_value_transition_t T[1];
    double decoded[3][4];

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, paper_ex_sites,
        paper_ex_mutations, paper_ex_individuals, NULL, 0);

    ret = tsk_compressed_matrix_init(&forward, &ts, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* We want a tree node that is not in the first tree */
    T[0].tree_node = 7;
    T[0].value = 0;
    ret = tsk_compressed_matrix_store_site(&forward, 0, 1, 1, T);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_compressed_matrix_decode(&forward, (double *) decoded);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_COMPRESSED_MATRIX_NODE);

    tsk_compressed_matrix_free(&forward);
    tsk_treeseq_free(&ts);
}

static void
test_caterpillar_tree_many_values(void)
{
    int ret = 0;
    tsk_ls_hmm_t ls_hmm;
    tsk_compressed_matrix_t matrix;
    double unused[] = { 0, 0, 0, 0, 0 };
    int32_t h[] = { 0, 0, 0, 0, 0 };
    tsk_size_t n[] = {
        8,
        16,
        32,
        64,
    };
    tsk_treeseq_t *ts;
    tsk_size_t j;

    for (j = 0; j < sizeof(n) / sizeof(*n); j++) {
        ts = caterpillar_tree(n[j], 5, n[j] - 2);
        ret = tsk_ls_hmm_init(&ls_hmm, ts, unused, unused, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_compressed_matrix_init(&matrix, ts, 1 << 10, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = run_test_hmm(&ls_hmm, h, &matrix);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tsk_compressed_matrix_print_state(&matrix, _devnull);
        tsk_ls_hmm_print_state(&ls_hmm, _devnull);

        tsk_ls_hmm_free(&ls_hmm);
        tsk_compressed_matrix_free(&matrix);
        tsk_treeseq_free(ts);
        free(ts);
    }

    j = 40;
    ts = caterpillar_tree(j, 5, j - 2);
    ret = tsk_ls_hmm_init(&ls_hmm, ts, unused, unused, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_compressed_matrix_init(&matrix, ts, 1 << 20, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Short circuit this value so we can run the test in reasonable time */
    ls_hmm.max_parsimony_words = 1;
    ret = run_test_hmm(&ls_hmm, h, &matrix);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TOO_MANY_VALUES);

    tsk_ls_hmm_free(&ls_hmm);
    tsk_compressed_matrix_free(&matrix);
    tsk_treeseq_free(ts);
    free(ts);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_single_tree_missing_alleles", test_single_tree_missing_alleles },
        { "test_single_tree_exact_match", test_single_tree_exact_match },
        { "test_single_tree_missing_haplotype_data",
            test_single_tree_missing_haplotype_data },
        { "test_single_tree_match_impossible", test_single_tree_match_impossible },
        { "test_single_tree_errors", test_single_tree_errors },
        { "test_single_tree_compressed_matrix", test_single_tree_compressed_matrix },
        { "test_single_tree_viterbi_matrix", test_single_tree_viterbi_matrix },

        { "test_multi_tree_exact_match", test_multi_tree_exact_match },
        { "test_multi_tree_errors", test_multi_tree_errors },

        { "test_caterpillar_tree_many_values", test_caterpillar_tree_many_values },
        { NULL, NULL },
    };

    return test_main(tests, argc, argv);
}
