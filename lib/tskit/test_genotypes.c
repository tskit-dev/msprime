#include "testlib.h"
#include "tsk_genotypes.h"

#include <unistd.h>
#include <stdlib.h>

static void
test_single_tree_hapgen_char_alphabet(void)
{
    int ret = 0;
    const char *sites =
        "0.0    A\n"
        "0.1    A\n"
        "0.2    C\n"
        "0.4    A\n";
    const char *mutations =
        "0    0     T\n"
        "1    1     T\n"
        "2    0     G\n"
        "2    1     A\n"
        "2    2     T\n"  // A bunch of different sample mutations
        "3    4     T\n"
        "3    0     A\n"; // A back mutation from T -> A
    uint32_t num_samples = 4;
    tsk_treeseq_t ts;
    char *haplotype;
    size_t j;
    tsk_hapgen_t hapgen;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
            NULL, NULL, NULL, NULL);
    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < num_samples; j++) {
        ret = tsk_hapgen_get_haplotype(&hapgen, (tsk_id_t) j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, "");
    }
    ret = tsk_hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
            sites, mutations, NULL, NULL);
    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);

    ret = tsk_hapgen_get_haplotype(&hapgen, 0, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "TAGA");
    ret = tsk_hapgen_get_haplotype(&hapgen, 1, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "ATAT");
    ret = tsk_hapgen_get_haplotype(&hapgen, 2, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "AATA");
    ret = tsk_hapgen_get_haplotype(&hapgen, 3, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "AACA");

    ret = tsk_hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_treeseq_free(&ts);
}

static void
test_single_tree_vargen_char_alphabet(void)
{
    int ret = 0;
    const char *sites =
        "0.0    A\n"
        "0.1    A\n"
        "0.2    C\n"
        "0.4    A\n";
    const char *mutations =
        "0    0     T   -1\n"
        "1    1     TTTAAGGG   -1\n"
        "2    0     G   -1\n"
        "2    1     AT  -1\n"
        "2    2     T   -1\n"  // A bunch of different sample mutations
        "3    4     T   -1\n"
        "3    0     A   5\n"; // A back mutation from T -> A
    tsk_treeseq_t ts;
    tsk_vargen_t vargen;
    tsk_variant_t *var;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
            sites, mutations, NULL, NULL);
    ret = tsk_vargen_alloc(&vargen, &ts, NULL, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->site->position, 0.0);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_EQUAL(var->allele_lengths[0], 1);
    CU_ASSERT_EQUAL(var->allele_lengths[1], 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "A", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "T", 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[2], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[3], 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->site->position, 0.1);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_EQUAL(var->allele_lengths[0], 1);
    CU_ASSERT_EQUAL(var->allele_lengths[1], 8);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "A", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "TTTAAGGG", 8);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[2], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[3], 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->site->position, 0.2);
    CU_ASSERT_EQUAL(var->num_alleles, 4);
    CU_ASSERT_EQUAL(var->allele_lengths[0], 1);
    CU_ASSERT_EQUAL(var->allele_lengths[1], 1);
    CU_ASSERT_EQUAL(var->allele_lengths[2], 2);
    CU_ASSERT_EQUAL(var->allele_lengths[3], 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "C", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "G", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[2], "AT", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[3], "T", 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 2);
    CU_ASSERT_EQUAL(var->genotypes.u8[2], 3);
    CU_ASSERT_EQUAL(var->genotypes.u8[3], 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->site->position, 0.4);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_EQUAL(var->allele_lengths[0], 1);
    CU_ASSERT_EQUAL(var->allele_lengths[1], 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "A", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "T", 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[2], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[3], 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_vargen_free(&vargen);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_hapgen_binary_alphabet(void)
{
    int ret = 0;
    uint32_t num_samples = 4;
    tsk_treeseq_t ts;
    char *haplotype;
    size_t j;
    tsk_hapgen_t hapgen;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
            NULL, NULL, NULL, NULL);
    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < num_samples; j++) {
        ret = tsk_hapgen_get_haplotype(&hapgen, (tsk_id_t) j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, "");
    }
    ret = tsk_hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
            single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL);
    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);

    ret = tsk_hapgen_get_haplotype(&hapgen, 0, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "001");
    ret = tsk_hapgen_get_haplotype(&hapgen, 1, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "011");
    ret = tsk_hapgen_get_haplotype(&hapgen, 2, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "101");
    ret = tsk_hapgen_get_haplotype(&hapgen, 3, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "001");

    ret = tsk_hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_treeseq_free(&ts);
}

static void
test_single_tree_vargen_binary_alphabet(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_vargen_t vargen;
    tsk_variant_t *var;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
            single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL);
    ret = tsk_vargen_alloc(&vargen, &ts, NULL, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_vargen_print_state(&vargen, _devnull);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[2], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[3], 0);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->site->id, 0);
    CU_ASSERT_EQUAL(var->site->mutations_length, 1);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[2], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[3], 0);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->site->id, 1);
    CU_ASSERT_EQUAL(var->site->mutations_length, 2);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[2], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[3], 1);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->site->id, 2);
    CU_ASSERT_EQUAL(var->site->mutations_length, 4);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_vargen_free(&vargen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_vargen_errors(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_vargen_t vargen;
    tsk_id_t samples[] = {0, 3};

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
            single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL);
    ret = tsk_vargen_alloc(&vargen, &ts, samples, 2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_vargen_free(&vargen);

    samples[0] = -1;
    ret = tsk_vargen_alloc(&vargen, &ts, samples, 2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_OUT_OF_BOUNDS);
    tsk_vargen_free(&vargen);

    samples[0] = 7;
    ret = tsk_vargen_alloc(&vargen, &ts, samples, 2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_OUT_OF_BOUNDS);
    tsk_vargen_free(&vargen);

    samples[0] = 3;
    ret = tsk_vargen_alloc(&vargen, &ts, samples, 2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);
    tsk_vargen_free(&vargen);

    tsk_treeseq_free(&ts);
}

static void
test_single_tree_vargen_subsample(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_vargen_t vargen;
    tsk_variant_t *var;
    tsk_id_t samples[] = {0, 3};

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
            single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL);
    ret = tsk_vargen_alloc(&vargen, &ts, samples, 2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_vargen_print_state(&vargen, _devnull);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 0);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->site->id, 0);
    CU_ASSERT_EQUAL(var->site->mutations_length, 1);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 0);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->site->id, 1);
    CU_ASSERT_EQUAL(var->site->mutations_length, 2);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 1);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->site->id, 2);
    CU_ASSERT_EQUAL(var->site->mutations_length, 4);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_vargen_free(&vargen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Zero samples */
    ret = tsk_vargen_alloc(&vargen, &ts, samples, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_vargen_print_state(&vargen, _devnull);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->site->id, 0);
    CU_ASSERT_EQUAL(var->site->mutations_length, 1);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->site->id, 1);
    CU_ASSERT_EQUAL(var->site->mutations_length, 2);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->site->id, 2);
    CU_ASSERT_EQUAL(var->site->mutations_length, 4);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_vargen_free(&vargen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_treeseq_free(&ts);
}

static void
test_single_tree_vargen_many_alleles(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_vargen_t vargen;
    tsk_variant_t *var;
    tsk_tbl_size_t num_alleles = 257;
    tsk_id_t j, k, l;
    int flags;
    char alleles[num_alleles];
    tsk_tbl_collection_t tables;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
            NULL, NULL, NULL, NULL);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_FATAL(ret == 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_FATAL(ret == 0);
    tsk_treeseq_free(&ts);
    memset(alleles, 'X', (size_t) num_alleles);
    ret = tsk_site_tbl_add_row(tables.sites, 0, "Y", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);

    /* Add j mutations over a single node. */
    for (j = 0; j < (tsk_id_t) num_alleles; j++) {
        /* When j = 0 we get a parent of -1, which is the NULL_NODE */
        ret = tsk_mutation_tbl_add_row(tables.mutations, 0, 0, j - 1, alleles,
                (tsk_tbl_size_t) j, NULL, 0);
        CU_ASSERT_FATAL(ret >= 0);
        ret = tsk_treeseq_alloc(&ts, &tables, TSK_BUILD_INDEXES);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (l = 0; l < 2; l++) {
            flags = 0;
            if (l == 1) {
                flags = TSK_16_BIT_GENOTYPES;
            }
            ret = tsk_vargen_alloc(&vargen, &ts, NULL, 0, flags);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            tsk_vargen_print_state(&vargen, _devnull);
            ret = tsk_vargen_next(&vargen, &var);
            /* We have j + 2 alleles. So, if j >= 254, we should fail with 8bit
             * genotypes */
            if (l == 0 && j >= 254) {
                CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TOO_MANY_ALLELES);
            } else {
                CU_ASSERT_EQUAL_FATAL(ret, 1);
                CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "Y", 1);
                for (k = 1; k < (tsk_id_t) var->num_alleles; k++) {
                    CU_ASSERT_EQUAL(k - 1, var->allele_lengths[k]);
                    CU_ASSERT_NSTRING_EQUAL(var->alleles[k], alleles, var->allele_lengths[k]);
                }
                CU_ASSERT_EQUAL(var->num_alleles, j + 2);
            }
            ret = tsk_vargen_free(&vargen);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
        }
        tsk_treeseq_free(&ts);
    }
    tsk_tbl_collection_free(&tables);
}

static void
test_single_unary_tree_hapgen(void)
{
    int ret = 0;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  1   0\n"
        "0  2   0\n"
        "0  3   0\n"
        "0  4   0\n";
    const char *edges =
        "0 1 2 0\n"
        "0 1 3 1\n"
        "0 1 4 2,3\n"
        "0 1 5 4\n"
        "0 1 6 5\n";
    const char *sites =
        "0     0\n"
        "0.1   0\n"
        "0.2   0\n"
        "0.3   0\n";
    const char *mutations =
        "0    0     1\n"
        "1    1     1\n"
        "2    4     1\n"
        "3    5     1\n";
    tsk_treeseq_t ts;
    size_t num_samples = 2;
    size_t j;
    tsk_hapgen_t hapgen;
    char *haplotype;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL);

    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < num_samples; j++) {
        ret = tsk_hapgen_get_haplotype(&hapgen, (tsk_id_t) j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, "");
    }
    ret = tsk_hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, sites, mutations, NULL, NULL);
    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);

    ret = tsk_hapgen_get_haplotype(&hapgen, 0, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "1011");
    ret = tsk_hapgen_get_haplotype(&hapgen, 1, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "0111");

    ret = tsk_hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_treeseq_free(&ts);
}

static void
test_single_tree_inconsistent_mutations(void)
{
    const char *sites =
        "0.0     0\n"
        "0.1     0\n"
        "0.2     0\n";
    const char *mutations =
        "0    0     1\n"
        "1    1     1\n"
        "2    4     1\n"
        "2    0     1\n";
    tsk_treeseq_t ts;
    tsk_variant_t *var;
    tsk_vargen_t vargen;
    tsk_hapgen_t hapgen;
    int flags[] = {0, TSK_16_BIT_GENOTYPES};
    tsk_id_t all_samples[] = {0, 1, 2, 3};
    tsk_id_t *samples[] = {NULL, all_samples};
    size_t num_samples = 4;
    size_t s, f;
    int ret;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
            sites, mutations, NULL, NULL);

    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INCONSISTENT_MUTATIONS);
    ret = tsk_hapgen_free(&hapgen);

    for (s = 0; s < 2; s++) {
        for (f = 0; f < sizeof(flags) / sizeof(*flags); f++) {
            ret = tsk_vargen_alloc(&vargen, &ts, samples[s], num_samples, flags[f]);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            ret = tsk_vargen_next(&vargen, &var);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
            ret = tsk_vargen_next(&vargen, &var);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
            ret = tsk_vargen_next(&vargen, &var);
            CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INCONSISTENT_MUTATIONS);
            tsk_vargen_free(&vargen);
        }
    }

    tsk_treeseq_free(&ts);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        {"test_single_tree_hapgen_char_alphabet", test_single_tree_hapgen_char_alphabet},
        {"test_single_tree_hapgen_binary_alphabet", test_single_tree_hapgen_binary_alphabet},
        {"test_single_unary_tree_hapgen", test_single_unary_tree_hapgen},
        {"test_single_tree_vargen_char_alphabet", test_single_tree_vargen_char_alphabet},
        {"test_single_tree_vargen_binary_alphabet", test_single_tree_vargen_binary_alphabet},
        {"test_single_tree_vargen_errors", test_single_tree_vargen_errors},
        {"test_single_tree_vargen_subsample", test_single_tree_vargen_subsample},
        {"test_single_tree_vargen_many_alleles", test_single_tree_vargen_many_alleles},
        {"test_single_tree_inconsistent_mutations", test_single_tree_inconsistent_mutations},
        {NULL},
    };

    return test_main(tests, argc, argv);
}
