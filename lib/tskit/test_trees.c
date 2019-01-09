#include "testlib.h"
#include "tsk_trees.h"
#include "tsk_genotypes.h"

#include <unistd.h>
#include <stdlib.h>


/*=======================================================
 * Verification utilities.
 *======================================================*/

static void
verify_compute_mutation_parents(tsk_treeseq_t *ts)
{
    int ret;
    size_t size = tsk_treeseq_get_num_mutations(ts) * sizeof(tsk_id_t);
    tsk_id_t *parent = malloc(size);
    tsk_tbl_collection_t tables;

    CU_ASSERT_FATAL(parent != NULL);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_copy_tables(ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    memcpy(parent, tables.mutations->parent, size);
    /* tsk_tbl_collection_print_state(&tables, stdout); */
    /* Make sure the tables are actually updated */
    memset(tables.mutations->parent, 0xff, size);

    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(memcmp(parent, tables.mutations->parent, size), 0);
    /* printf("after\n"); */
    /* tsk_tbl_collection_print_state(&tables, stdout); */

    free(parent);
    tsk_tbl_collection_free(&tables);
}

static void
verify_individual_nodes(tsk_treeseq_t *ts)
{
    int ret;
    tsk_individual_t individual;
    tsk_id_t k;
    size_t num_nodes = tsk_treeseq_get_num_nodes(ts);
    size_t num_individuals = tsk_treeseq_get_num_individuals(ts);
    size_t j;

    for (k = 0; k < (tsk_id_t) num_individuals; k++) {
        ret = tsk_treeseq_get_individual(ts, (size_t) k, &individual);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (j = 0; j < individual.nodes_length; j++) {
            CU_ASSERT_FATAL(individual.nodes[j] < (tsk_id_t) num_nodes);
            CU_ASSERT_EQUAL_FATAL(k,
                    ts->tables->nodes->individual[individual.nodes[j]]);
        }
    }
}

static void
verify_trees(tsk_treeseq_t *ts, uint32_t num_trees, tsk_id_t* parents)
{
    int ret;
    tsk_id_t u, j, v;
    uint32_t mutation_index, site_index;
    tsk_tbl_size_t k, l, tree_sites_length;
    tsk_site_t *sites = NULL;
    tsk_tree_t tree;
    size_t num_nodes = tsk_treeseq_get_num_nodes(ts);
    size_t num_sites = tsk_treeseq_get_num_sites(ts);
    size_t num_mutations = tsk_treeseq_get_num_mutations(ts);

    ret = tsk_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(ts), num_trees);

    site_index = 0;
    mutation_index = 0;
    j = 0;
    for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
        CU_ASSERT_EQUAL(j, tree.index);
        tsk_tree_print_state(&tree, _devnull);
        /* tsk_tree_print_state(&tree, stdout); */
        for (u = 0; u < (tsk_id_t) num_nodes; u++) {
            ret = tsk_tree_get_parent(&tree, u, &v);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(v, parents[j * (tsk_id_t) num_nodes + u]);
        }
        ret = tsk_tree_get_sites(&tree, &sites, &tree_sites_length);
        CU_ASSERT_EQUAL(ret, 0);
        for (k = 0; k < tree_sites_length; k++) {
            CU_ASSERT_EQUAL(sites[k].id, site_index);
            for (l = 0; l < sites[k].mutations_length; l++) {
                CU_ASSERT_EQUAL(sites[k].mutations[l].id, mutation_index);
                CU_ASSERT_EQUAL(sites[k].mutations[l].site, site_index);
                mutation_index++;
            }
            site_index++;
        }
        j++;
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(site_index, num_sites);
    CU_ASSERT_EQUAL(mutation_index, num_mutations);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_tree_free(&tree);
}

static tsk_tree_t *
get_tree_list(tsk_treeseq_t *ts)
{
    int ret;
    tsk_tree_t t, *trees;
    size_t num_trees;

    num_trees = tsk_treeseq_get_num_trees(ts);
    ret = tsk_tree_alloc(&t, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    trees = malloc(num_trees * sizeof(tsk_tree_t));
    CU_ASSERT_FATAL(trees != NULL);
    for (ret = tsk_tree_first(&t); ret == 1; ret = tsk_tree_next(&t)) {
        CU_ASSERT_FATAL(t.index < num_trees);
        ret = tsk_tree_alloc(&trees[t.index], ts, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_copy(&trees[t.index], &t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tree_equal(&trees[t.index], &t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        /* Make sure the left and right coordinates are also OK */
        CU_ASSERT_DOUBLE_EQUAL(trees[t.index].left, t.left, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(trees[t.index].right, t.right, 1e-6);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    return trees;
}


static void
verify_tree_next_prev(tsk_treeseq_t *ts)
{
    int ret;
    tsk_tree_t *trees, t;
    size_t j;
    size_t num_trees = tsk_treeseq_get_num_trees(ts);

    trees = get_tree_list(ts);
    ret = tsk_tree_alloc(&t, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Single forward pass */
    j = 0;
    for (ret = tsk_tree_first(&t); ret == 1; ret = tsk_tree_next(&t)) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);

    /* Single reverse pass */
    j = num_trees;
    for (ret = tsk_tree_last(&t); ret == 1; ret = tsk_tree_prev(&t)) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        if (ret != 0) {
            printf("trees differ\n");
            printf("REVERSE tree::\n");
            tsk_tree_print_state(&t, stdout);
            printf("FORWARD tree::\n");
            tsk_tree_print_state(&trees[t.index], stdout);
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);

    /* Full forward, then reverse */
    j = 0;
    for (ret = tsk_tree_first(&t); ret == 1; ret = tsk_tree_next(&t)) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);
    j--;
    while ((ret = tsk_tree_prev(&t)) == 1) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);
    CU_ASSERT_EQUAL_FATAL(t.index, 0);
    /* Calling prev should return 0 and have no effect. */
    for (j = 0; j < 10; j++) {
        ret = tsk_tree_prev(&t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, 0);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Full reverse then forward */
    j = num_trees;
    for (ret = tsk_tree_last(&t); ret == 1; ret = tsk_tree_prev(&t)) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);
    j++;
    while ((ret = tsk_tree_next(&t)) == 1) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);
    CU_ASSERT_EQUAL_FATAL(t.index, num_trees - 1);
    /* Calling next should return 0 and have no effect. */
    for (j = 0; j < 10; j++) {
        ret = tsk_tree_next(&t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, num_trees - 1);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Do a zigzagging traversal */
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 1; j < TSK_MIN(10, num_trees / 2); j++) {
        while (t.index < num_trees - j) {
            ret = tsk_tree_next(&t);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        CU_ASSERT_EQUAL_FATAL(t.index, num_trees - j);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        while (t.index > j) {
            ret = tsk_tree_prev(&t);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        CU_ASSERT_EQUAL_FATAL(t.index, j);
        ret = tsk_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Free the trees. */
    ret = tsk_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < tsk_treeseq_get_num_trees(ts); j++) {
        ret = tsk_tree_free(&trees[j]);
    }
    free(trees);
}

static void
verify_tree_diffs(tsk_treeseq_t *ts)
{
    int ret;
    tsk_diff_iter_t iter;
    tsk_tree_t tree;
    tsk_edge_list_t *record, *records_out, *records_in;
    size_t num_nodes = tsk_treeseq_get_num_nodes(ts);
    size_t j, num_trees;
    double left, right;
    tsk_id_t *parent = malloc(num_nodes * sizeof(tsk_id_t));
    tsk_id_t *child = malloc(num_nodes * sizeof(tsk_id_t));
    tsk_id_t *sib = malloc(num_nodes * sizeof(tsk_id_t));
    tsk_id_t *samples;

    CU_ASSERT_FATAL(parent != NULL);
    CU_ASSERT_FATAL(child != NULL);
    CU_ASSERT_FATAL(sib != NULL);
    for (j = 0; j < num_nodes; j++) {
        parent[j] = TSK_NULL;
        child[j] = TSK_NULL;
        sib[j] = TSK_NULL;
    }
    ret = tsk_treeseq_get_samples(ts, &samples);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_diff_iter_alloc(&iter, ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    tsk_diff_iter_print_state(&iter, _devnull);

    num_trees = 0;
    while ((ret = tsk_diff_iter_next(
                &iter, &left, &right, &records_out, &records_in)) == 1) {
        tsk_diff_iter_print_state(&iter, _devnull);
        num_trees++;
        for (record = records_out; record != NULL; record = record->next) {
            parent[record->edge.child] = TSK_NULL;
        }
        for (record = records_in; record != NULL; record = record->next) {
            parent[record->edge.child] = record->edge.parent;
        }
        /* Now check against the sparse tree iterator. */
        for (j = 0; j < num_nodes; j++) {
            CU_ASSERT_EQUAL(parent[j], tree.parent[j]);
        }
        CU_ASSERT_EQUAL(tree.left, left);
        CU_ASSERT_EQUAL(tree.right, right);
        ret = tsk_tree_next(&tree);
        if (num_trees < tsk_treeseq_get_num_trees(ts)) {
            CU_ASSERT_EQUAL(ret, 1);
        } else {
            CU_ASSERT_EQUAL(ret, 0);
        }
    }
    CU_ASSERT_EQUAL(num_trees, tsk_treeseq_get_num_trees(ts));
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_diff_iter_free(&iter);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_free(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    free(parent);
    free(child);
    free(sib);
}

/* When we keep all sites in simplify, the genotypes for the subset of the
 * samples should be the same as the original */
static void
verify_simplify_genotypes(tsk_treeseq_t *ts, tsk_treeseq_t *subset,
        tsk_id_t *samples, size_t num_samples)
{
    int ret;
    size_t m = tsk_treeseq_get_num_sites(ts);
    tsk_vargen_t vargen, subset_vargen;
    tsk_variant_t *variant, *subset_variant;
    size_t j, k;
    tsk_id_t *all_samples;
    uint8_t a1, a2;
    tsk_id_t *sample_index_map;

    tsk_treeseq_get_sample_index_map(ts, &sample_index_map);

    /* tsk_treeseq_print_state(ts, stdout); */
    /* tsk_treeseq_print_state(subset, stdout); */

    ret = tsk_vargen_alloc(&vargen, ts, NULL, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_vargen_alloc(&subset_vargen, subset, NULL, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(m, tsk_treeseq_get_num_sites(subset));
    tsk_treeseq_get_samples(ts, &all_samples);

    for (j = 0; j < m; j++) {
        ret = tsk_vargen_next(&vargen, &variant);
        CU_ASSERT_EQUAL_FATAL(ret, 1);
        ret = tsk_vargen_next(&subset_vargen, &subset_variant);
        CU_ASSERT_EQUAL_FATAL(ret, 1);
        CU_ASSERT_EQUAL(variant->site->id, j)
        CU_ASSERT_EQUAL(subset_variant->site->id, j)
        CU_ASSERT_EQUAL(variant->site->position, subset_variant->site->position);
        for (k = 0; k < num_samples; k++) {
            CU_ASSERT_FATAL(sample_index_map[samples[k]] < (tsk_id_t) ts->num_samples);
            a1 = variant->genotypes.u8[sample_index_map[samples[k]]];
            a2 = subset_variant->genotypes.u8[k];
            /* printf("a1 = %d, a2 = %d\n", a1, a2); */
            /* printf("k = %d original node = %d " */
            /*         "original_index = %d a1=%.*s a2=%.*s\n", */
            /*         (int) k, samples[k], sample_index_map[samples[k]], */
            /*         variant->allele_lengths[a1], variant->alleles[a1], */
            /*         subset_variant->allele_lengths[a2], subset_variant->alleles[a2]); */
            CU_ASSERT_FATAL(a1 < variant->num_alleles);
            CU_ASSERT_FATAL(a2 < subset_variant->num_alleles);
            CU_ASSERT_EQUAL_FATAL(variant->allele_lengths[a1],
                    subset_variant->allele_lengths[a2]);
            CU_ASSERT_NSTRING_EQUAL_FATAL(
                variant->alleles[a1], subset_variant->alleles[a2],
                variant->allele_lengths[a1]);
        }
    }
    tsk_vargen_free(&vargen);
    tsk_vargen_free(&subset_vargen);
}


static void
verify_simplify_properties(tsk_treeseq_t *ts, tsk_treeseq_t *subset,
        tsk_id_t *samples, size_t num_samples, tsk_id_t *node_map)
{
    int ret;
    tsk_node_t n1, n2;
    tsk_tree_t full_tree, subset_tree;
    tsk_site_t *tree_sites;
    tsk_tbl_size_t tree_sites_length;
    uint32_t j, k;
    tsk_id_t u, mrca1, mrca2;
    size_t total_sites;

    CU_ASSERT_EQUAL(
        tsk_treeseq_get_sequence_length(ts),
        tsk_treeseq_get_sequence_length(subset));
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(subset), num_samples);
    CU_ASSERT(
        tsk_treeseq_get_num_nodes(ts) >= tsk_treeseq_get_num_nodes(subset));
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(subset), num_samples);

    /* Check the sample properties */
    for (j = 0; j < num_samples; j++) {
        ret = tsk_treeseq_get_node(ts, (size_t) samples[j], &n1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(node_map[samples[j]], j);
        ret = tsk_treeseq_get_node(subset, (size_t) node_map[samples[j]], &n2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(n1.population, n2.population);
        CU_ASSERT_EQUAL_FATAL(n1.time, n2.time);
        CU_ASSERT_EQUAL_FATAL(n1.flags, n2.flags);
        CU_ASSERT_EQUAL_FATAL(n1.metadata_length, n2.metadata_length);
        CU_ASSERT_NSTRING_EQUAL(n1.metadata, n2.metadata, n2.metadata_length);
    }
    /* Check that node mappings are correct */
    for (j = 0; j < tsk_treeseq_get_num_nodes(ts); j++) {
        ret = tsk_treeseq_get_node(ts, j, &n1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        if (node_map[j] != TSK_NULL) {
            ret = tsk_treeseq_get_node(subset, (size_t) node_map[j], &n2);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(n1.population, n2.population);
            CU_ASSERT_EQUAL_FATAL(n1.time, n2.time);
            CU_ASSERT_EQUAL_FATAL(n1.flags, n2.flags);
            CU_ASSERT_EQUAL_FATAL(n1.metadata_length, n2.metadata_length);
            CU_ASSERT_NSTRING_EQUAL(n1.metadata, n2.metadata, n2.metadata_length);
        }
    }
    if (num_samples == 0) {
        CU_ASSERT_EQUAL(tsk_treeseq_get_num_edges(subset), 0);
        CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(subset), 0);
    } else if (num_samples == 1) {
        CU_ASSERT_EQUAL(tsk_treeseq_get_num_edges(subset), 0);
        CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(subset), 1);
    }
    /* Check the pairwise MRCAs */
    ret = tsk_tree_alloc(&full_tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_alloc(&subset_tree, subset, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&full_tree);
    CU_ASSERT_EQUAL(ret, 1);
    ret = tsk_tree_first(&subset_tree);
    CU_ASSERT_EQUAL(ret, 1);

    total_sites = 0;
    while (1) {
        while (full_tree.right <= subset_tree.right) {
            for (j = 0; j < num_samples; j++) {
                for (k = j + 1; k < num_samples; k++) {
                    ret = tsk_tree_get_mrca(&full_tree, samples[j], samples[k], &mrca1);
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                    ret = tsk_tree_get_mrca(&subset_tree,
                            node_map[samples[j]], node_map[samples[k]], &mrca2);
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                    if (mrca1 == TSK_NULL) {
                        CU_ASSERT_EQUAL_FATAL(mrca2, TSK_NULL);
                    } else {
                        CU_ASSERT_EQUAL(node_map[mrca1], mrca2);
                    }
                }
            }
            ret = tsk_tree_next(&full_tree);
            CU_ASSERT_FATAL(ret >= 0);
            if (ret != 1) {
                break;
            }
        }
        /* Check the sites in this tree */
        ret = tsk_tree_get_sites(&subset_tree, &tree_sites, &tree_sites_length);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        for (j = 0; j < tree_sites_length; j++) {
            CU_ASSERT(subset_tree.left <= tree_sites[j].position);
            CU_ASSERT(tree_sites[j].position < subset_tree.right);
            for (k = 0; k < tree_sites[j].mutations_length; k++) {
                ret = tsk_tree_get_parent(&subset_tree,
                        tree_sites[j].mutations[k].node, &u);
                CU_ASSERT_EQUAL(ret, 0);
            }
            total_sites++;
        }
        ret = tsk_tree_next(&subset_tree);
        if (ret != 1) {
            break;
        }
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(subset), total_sites);

    tsk_tree_free(&subset_tree);
    tsk_tree_free(&full_tree);
}

static void
verify_simplify(tsk_treeseq_t *ts)
{
    int ret;
    size_t n = tsk_treeseq_get_num_samples(ts);
    size_t num_samples[] = {0, 1, 2, 3, n / 2, n - 1, n};
    size_t j;
    tsk_id_t *sample;
    tsk_id_t *node_map = malloc(tsk_treeseq_get_num_nodes(ts) * sizeof(tsk_id_t));
    tsk_treeseq_t subset;
    int flags = TSK_FILTER_SITES;

    CU_ASSERT_FATAL(node_map != NULL);
    ret = tsk_treeseq_get_samples(ts, &sample);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    if (tsk_treeseq_get_num_migrations(ts) > 0) {
        ret = tsk_treeseq_simplify(ts, sample, 2, 0, &subset, NULL);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SIMPLIFY_MIGRATIONS_NOT_SUPPORTED);
        /* Exiting early here because simplify isn't supported with migrations. */
        goto out;
    }

    for (j = 0; j < sizeof(num_samples) / sizeof(*num_samples); j++) {
        if (num_samples[j] <= n) {
            ret = tsk_treeseq_simplify(ts, sample, num_samples[j], flags, &subset,
                    node_map);
            /* printf("ret = %s\n", tsk_strerror(ret)); */
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            verify_simplify_properties(ts, &subset, sample, num_samples[j], node_map);
            tsk_treeseq_free(&subset);

            /* Keep all sites */
            ret = tsk_treeseq_simplify(ts, sample, num_samples[j], 0, &subset,
                    node_map);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            verify_simplify_properties(ts, &subset, sample, num_samples[j], node_map);
            verify_simplify_genotypes(ts, &subset, sample, num_samples[j]);
            tsk_treeseq_free(&subset);
        }
    }
out:
    free(node_map);
}

typedef struct {
    uint32_t tree_index;
    tsk_id_t node;
    uint32_t count;
} sample_count_test_t;

static void
verify_sample_counts(tsk_treeseq_t *ts, size_t num_tests, sample_count_test_t *tests)
{
    int ret;
    size_t j, num_samples, n, k;
    tsk_id_t stop, sample_index;
    tsk_tree_t tree;
    tsk_id_t *samples;

    n = tsk_treeseq_get_num_samples(ts);
    ret = tsk_treeseq_get_samples(ts, &samples);
    CU_ASSERT_EQUAL(ret, 0);

    /* First run without the TSK_SAMPLE_COUNTS feature */
    ret = tsk_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = tsk_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = tsk_tree_get_num_samples(&tree, tests[j].node, &num_samples);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_samples);
        /* all operations depending on tracked samples should fail. */
        ret = tsk_tree_get_num_tracked_samples(&tree, 0, &num_samples);
        CU_ASSERT_EQUAL(ret, TSK_ERR_UNSUPPORTED_OPERATION);
    }
    tsk_tree_free(&tree);

    /* Now run with TSK_SAMPLE_COUNTS but with no samples tracked. */
    ret = tsk_tree_alloc(&tree, ts, TSK_SAMPLE_COUNTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = tsk_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = tsk_tree_get_num_samples(&tree, tests[j].node, &num_samples);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_samples);
        /* all operations depending on tracked samples should fail. */
        ret = tsk_tree_get_num_tracked_samples(&tree, 0, &num_samples);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(num_samples, 0);
    }
    tsk_tree_free(&tree);

    /* Run with TSK_SAMPLE_LISTS, but without TSK_SAMPLE_COUNTS */
    ret = tsk_tree_alloc(&tree, ts, TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = tsk_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = tsk_tree_get_num_samples(&tree, tests[j].node, &num_samples);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_samples);
        /* all operations depending on tracked samples should fail. */
        ret = tsk_tree_get_num_tracked_samples(&tree, 0, &num_samples);
        CU_ASSERT_EQUAL(ret, TSK_ERR_UNSUPPORTED_OPERATION);

        sample_index = tree.left_sample[tests[j].node];
        k = 0;
        if (sample_index != TSK_NULL) {
            stop = tree.right_sample[tests[j].node];
            while (true) {
                k++;
                CU_ASSERT_FATAL(k <= tests[j].count);
                if (sample_index == stop) {
                    break;
                }
                sample_index = tree.next_sample[sample_index];
            }
        }
        CU_ASSERT_EQUAL(tests[j].count, k);
    }
    tsk_tree_free(&tree);

    /* Now use TSK_SAMPLE_COUNTS|TSK_SAMPLE_LISTS */
    ret = tsk_tree_alloc(&tree, ts, TSK_SAMPLE_COUNTS|TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_set_tracked_samples(&tree, n, samples);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = tsk_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = tsk_tree_get_num_samples(&tree, tests[j].node, &num_samples);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_samples);

        /* We're tracking all samples, so the count should be the same */
        ret = tsk_tree_get_num_tracked_samples(&tree, tests[j].node, &num_samples);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_samples);

        sample_index = tree.left_sample[tests[j].node];
        k = 0;
        if (sample_index != TSK_NULL) {
            stop = tree.right_sample[tests[j].node];
            while (true) {
                k++;
                if (sample_index == stop) {
                    break;
                }
                sample_index = tree.next_sample[sample_index];
            }
        }
        CU_ASSERT_EQUAL(tests[j].count, k);
    }
    tsk_tree_free(&tree);
}

static void
verify_sample_sets_for_tree(tsk_tree_t *tree)
{
    int ret, stack_top, j;
    tsk_id_t u, v;
    size_t tmp, n, num_nodes, num_samples;
    tsk_id_t *stack, *samples;
    tsk_treeseq_t *ts = tree->tree_sequence;
    tsk_id_t *sample_index_map = ts->sample_index_map;
    const tsk_id_t *list_left = tree->left_sample;
    const tsk_id_t *list_right = tree->right_sample;
    const tsk_id_t *list_next = tree->next_sample;
    tsk_id_t stop, sample_index;

    n = tsk_treeseq_get_num_samples(ts);
    num_nodes = tsk_treeseq_get_num_nodes(ts);
    stack = malloc(n * sizeof(tsk_id_t));
    samples = malloc(n * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(stack != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    for (u = 0; u < (tsk_id_t) num_nodes; u++) {
        if (tree->left_child[u] == TSK_NULL && !tsk_treeseq_is_sample(ts, u)) {
            CU_ASSERT_EQUAL(list_left[u], TSK_NULL);
            CU_ASSERT_EQUAL(list_right[u], TSK_NULL);
        } else {
            stack_top = 0;
            num_samples = 0;
            stack[stack_top] = u;
            while (stack_top >= 0) {
                v = stack[stack_top];
                stack_top--;
                if (tsk_treeseq_is_sample(ts, v)) {
                    samples[num_samples] = v;
                    num_samples++;
                }
                for (v = tree->right_child[v]; v != TSK_NULL; v = tree->left_sib[v]) {
                    stack_top++;
                    stack[stack_top] = v;
                }
            }
            ret = tsk_tree_get_num_samples(tree, u, &tmp);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL_FATAL(num_samples, tmp);

            j = 0;
            sample_index = list_left[u];
            if (sample_index != TSK_NULL) {
                stop = list_right[u];
                while (true) {
                    CU_ASSERT_TRUE_FATAL(j < (tsk_id_t) n);
                    CU_ASSERT_EQUAL_FATAL(sample_index, sample_index_map[samples[j]]);
                    j++;
                    if (sample_index == stop) {
                        break;
                    }
                    sample_index = list_next[sample_index];
                }
            }
            CU_ASSERT_EQUAL_FATAL(j, num_samples);
        }
    }
    free(stack);
    free(samples);
}

static void
verify_sample_sets(tsk_treeseq_t *ts)
{
    int ret;
    tsk_tree_t t;

    ret = tsk_tree_alloc(&t, ts, TSK_SAMPLE_COUNTS|TSK_SAMPLE_LISTS);
    CU_ASSERT_EQUAL(ret, 0);

    for (ret = tsk_tree_first(&t); ret == 1; ret = tsk_tree_next(&t)) {
        verify_sample_sets_for_tree(&t);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (ret = tsk_tree_last(&t); ret == 1; ret = tsk_tree_prev(&t)) {
        verify_sample_sets_for_tree(&t);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_tree_free(&t);
}


static void
verify_empty_tree_sequence(tsk_treeseq_t *ts, double sequence_length)
{
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_edges(ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_migrations(ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(ts), sequence_length);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(ts), 1);
}

/*=======================================================
 * Simplest test cases.
 *======================================================*/

static void
test_simplest_records(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0";
    const char *edges =
        "0  1   2   0,1\n";
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_nonbinary_records(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0";
    const char *edges =
        "0  1   4   0,1,2,3\n";
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_unary_records(void)
{
    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  1   0\n"
        "0  2   0";
    const char *edges =
        "0  1   2   0\n"
        "0  1   3   1\n"
        "0  1   4   2,3\n";
    tsk_treeseq_t ts, simplified;
    tsk_id_t sample_ids[] = {0, 1};

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_populations(&ts), 1);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&simplified), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&simplified), 1);

    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&simplified);
}

static void
test_simplest_non_sample_leaf_records(void)
{
    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  0   0\n"
        "0  0   0";
    const char *edges =
        "0  1   2   0,1,3,4\n";
    const char *sites =
        "0.1  0\n"
        "0.2  0\n"
        "0.3  0\n"
        "0.4  0\n";
    const char *mutations =
        "0    0     1\n"
        "1    1     1\n"
        "2    3     1\n"
        "3    4     1";
    tsk_treeseq_t ts, simplified;
    tsk_id_t sample_ids[] = {0, 1};
    tsk_hapgen_t hapgen;
    tsk_vargen_t vargen;
    char *haplotype;
    tsk_variant_t *var;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, sites, mutations, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_hapgen_get_haplotype(&hapgen, 0, &haplotype);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "1000");
    ret = tsk_hapgen_get_haplotype(&hapgen, 1, &haplotype);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "0100");
    tsk_hapgen_free(&hapgen);

    ret = tsk_vargen_alloc(&vargen, &ts, NULL, 0, 0);
    tsk_vargen_print_state(&vargen, _devnull);
    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 1);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_vargen_free(&vargen);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&simplified), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&simplified), 1);

    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&simplified);
}

static void
test_simplest_degenerate_multiple_root_records(void)
{

    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  1   0\n";
    const char *edges =
        "0  1   2   0\n"
        "0  1   3   1\n";
    tsk_treeseq_t ts, simplified;
    tsk_tree_t t;
    tsk_id_t sample_ids[] = {0, 1};

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_tree_alloc(&t, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&t), 2);
    CU_ASSERT_EQUAL(t.left_root, 2);
    CU_ASSERT_EQUAL(t.right_sib[2], 3);
    CU_ASSERT_EQUAL(t.right_sib[3], TSK_NULL);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 2);

    tsk_treeseq_free(&simplified);
    tsk_treeseq_free(&ts);
    tsk_tree_free(&t);
}

static void
test_simplest_multiple_root_records(void)
{
    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  1   0\n";
    const char *edges =
        "0  1   4   0,1\n"
        "0  1   5   2,3\n";
    tsk_treeseq_t ts, simplified;
    tsk_id_t sample_ids[] = {0, 1, 2, 3};

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 4, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&simplified), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&simplified), 1);
    tsk_treeseq_free(&simplified);

    /* Make one tree degenerate */
    ret = tsk_treeseq_simplify(&ts, sample_ids, 3, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&simplified), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&simplified), 1);
    tsk_treeseq_free(&simplified);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_zero_root_tree(void)
{
    int ret;
    const char *nodes =
        "0  0   0\n"
        "0  0   0\n"
        "0  0   0\n"
        "0  0   0\n"
        "0  1   0\n"
        "0  1   0\n";
    const char *edges =
        "0  1   4   0,1\n"
        "0  1   5   2,3\n";
    tsk_treeseq_t ts;
    tsk_tree_t t;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_tree_alloc(&t, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&t), 0);
    CU_ASSERT_EQUAL(t.left_root, TSK_NULL);
    CU_ASSERT_EQUAL(t.right_sib[2], 3);
    CU_ASSERT_EQUAL(t.right_sib[3], TSK_NULL);

    tsk_tree_free(&t);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_root_mutations(void)
{
    int ret;
    uint32_t j;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n";
    const char *edges =
        "0  1   2   0,1\n";
    const char *sites =
        "0.1 0";
    const char *mutations =
        "0    2     1";
    tsk_hapgen_t hapgen;
    char *haplotype;
    int flags = 0;
    tsk_id_t sample_ids[] = {0, 1};
    tsk_treeseq_t ts, simplified;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, sites, mutations, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < 2; j++) {
        ret = tsk_hapgen_get_haplotype(&hapgen, (tsk_id_t) j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, "1");
    }
    tsk_hapgen_free(&hapgen);

    ret = tsk_treeseq_simplify(&ts, sample_ids, 2, flags, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&simplified), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&simplified), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&simplified), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&simplified), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&simplified), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&simplified), 1);
    tsk_treeseq_free(&simplified);

    tsk_treeseq_free(&ts);
}

static void
test_simplest_back_mutations(void)
{
    int ret;
    uint32_t j;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  2   0\n";
    const char *edges =
        "0  1   3   0,1\n"
        "0  1   4   2,3\n";
    const char *sites =
        "0.5 0";
    const char *mutations =
        "0    3     1   -1\n"
        "0    0     0   0";
    tsk_hapgen_t hapgen;
    const char *haplotypes[] = {"0", "1", "0"};
    char *haplotype;
    tsk_treeseq_t ts;
    tsk_vargen_t vargen;
    tsk_variant_t *var;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, sites, mutations, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < 2; j++) {
        ret = tsk_hapgen_get_haplotype(&hapgen, (tsk_id_t) j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, haplotypes[j]);
    }
    tsk_hapgen_free(&hapgen);

    ret = tsk_vargen_alloc(&vargen, &ts, NULL, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_vargen_next(&vargen, &var);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(var->num_alleles, 2);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[0], "0", 1);
    CU_ASSERT_NSTRING_EQUAL(var->alleles[1], "1", 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[0], 0);
    CU_ASSERT_EQUAL(var->genotypes.u8[1], 1);
    CU_ASSERT_EQUAL(var->genotypes.u8[2], 0);
    CU_ASSERT_EQUAL(var->site->id, 0);
    CU_ASSERT_EQUAL(var->site->mutations_length, 2);
    tsk_vargen_free(&vargen);

    tsk_treeseq_free(&ts);
}

static void
test_simplest_general_samples(void)
{
    const char *nodes =
        "1  0   0\n"
        "0  1   0\n"
        "1  0   0";
    const char *edges =
        "0  1   1   0,2\n";
    const char *sites =
        "0.5  0\n"
        "0.75 0\n";
    const char *mutations =
        "0    2     1\n"
        "1    0     1";
    const char *haplotypes[] = {"01", "10"};
    char *haplotype;
    unsigned int j;
    tsk_id_t samples[2] = {0, 2};
    tsk_id_t *s;
    int ret;

    tsk_treeseq_t ts, simplified;
    tsk_hapgen_t hapgen;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, sites, mutations, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_treeseq_get_samples(&ts, &s);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_FATAL(s != NULL);
    CU_ASSERT_EQUAL(s[0], 0);
    CU_ASSERT_EQUAL(s[1], 2);

    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < 2; j++) {
        ret = tsk_hapgen_get_haplotype(&hapgen, (tsk_id_t) j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, haplotypes[j]);
    }
    tsk_hapgen_free(&hapgen);

    ret = tsk_treeseq_simplify(&ts, samples, 2, 0, &simplified, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_get_samples(&simplified, &s);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_FATAL(s != NULL);
    CU_ASSERT_EQUAL(s[0], 0);
    CU_ASSERT_EQUAL(s[1], 1);

    ret = tsk_hapgen_alloc(&hapgen, &simplified);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < 2; j++) {
        ret = tsk_hapgen_get_haplotype(&hapgen, (tsk_id_t) j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, haplotypes[j]);
    }
    tsk_hapgen_free(&hapgen);

    tsk_treeseq_free(&simplified);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_holey_tree_sequence(void)
{
    const char *nodes_txt =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0";
    const char *edges_txt =
        "0  1   2   0\n"
        "2  3   2   0\n"
        "0  1   2   1\n"
        "2  3   2   1\n";
    const char *sites_txt =
        "0.5  0\n"
        "1.5  0\n"
        "2.5  0\n";
    const char *mutations_txt =
        "0    0     1\n"
        "1    1     1\n"
        "2    2     1\n";
    const char *haplotypes[] = {"101", "011"};
    char *haplotype;
    unsigned int j;
    int ret;
    tsk_treeseq_t ts;
    tsk_hapgen_t hapgen;

    tsk_treeseq_from_text(&ts, 3, nodes_txt, edges_txt, NULL, sites_txt,
            mutations_txt, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 3.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 3);

    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < 2; j++) {
        ret = tsk_hapgen_get_haplotype(&hapgen, (tsk_id_t) j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, haplotypes[j]);
    }

    tsk_hapgen_free(&hapgen);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_holey_tsk_treeseq_mutation_parents(void)
{
    const char *nodes_txt =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0";
    const char *edges_txt =
        "0  1   2   0\n"
        "2  3   2   0\n"
        "0  1   2   1\n"
        "2  3   2   1\n";
    const char *sites_txt =
        "0.5  0\n"
        "1.5  0\n"
        "2.5  0\n";
    const char *mutations_txt =
        "0    0     1\n"
        "0    0     1\n"
        "1    1     1\n"
        "1    1     1\n"
        "2    2     1\n"
        "2    2     1\n";
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    int ret;

    tsk_treeseq_from_text(&ts, 3, nodes_txt, edges_txt, NULL, sites_txt,
            mutations_txt, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 3);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.mutations->parent[0], -1);
    CU_ASSERT_EQUAL(tables.mutations->parent[1], 0);
    CU_ASSERT_EQUAL(tables.mutations->parent[2], -1);
    CU_ASSERT_EQUAL(tables.mutations->parent[3], 2);
    CU_ASSERT_EQUAL(tables.mutations->parent[4], -1);
    CU_ASSERT_EQUAL(tables.mutations->parent[5], 4);
    tsk_tbl_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_initial_gap_tree_sequence(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0";
    const char *edges =
        "2  3   2   0,1\n";
    const char *sites =
        "0.5  0\n"
        "1.5  0\n"
        "2.5  0\n";
    const char *mutations =
        "0    0     1\n"
        "1    1     1\n"
        "2    2     1";
    const char *haplotypes[] = {"101", "011"};
    char *haplotype;
    unsigned int j;
    int ret;
    tsk_treeseq_t ts;
    tsk_hapgen_t hapgen;
    const tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
        z, z, z,
        2, 2, z,
    };
    uint32_t num_trees = 2;

    tsk_treeseq_from_text(&ts, 3, nodes, edges, NULL, sites, mutations, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 3.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 2);

    verify_trees(&ts, num_trees, parents);

    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < 2; j++) {
        ret = tsk_hapgen_get_haplotype(&hapgen, (tsk_id_t) j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, haplotypes[j]);
    }
    tsk_hapgen_free(&hapgen);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_initial_gap_zero_roots(void)
{
    const char *nodes =
        "0  0   0\n"
        "0  0   0\n"
        "0  1   0";
    const char *edges =
        "2  3   2   0,1\n";
    int ret;
    tsk_treeseq_t ts;
    const tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
        z, z, z,
        2, 2, z,
    };
    uint32_t num_trees = 2;
    tsk_tree_t tree;

    tsk_treeseq_from_text(&ts, 3, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 3.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 2);

    verify_trees(&ts, num_trees, parents);

    ret = tsk_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tree.left_root, TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 0);
    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tree.left_root, TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 0);
    CU_ASSERT_EQUAL(tree.parent[0], 2);
    CU_ASSERT_EQUAL(tree.parent[1], 2);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_holey_tsk_treeseq_zero_roots(void)
{
    const char *nodes_txt =
        "0  0   0\n"
        "0  0   0\n"
        "0  1   0";
    const char *edges_txt =
        "0  1   2   0\n"
        "2  3   2   0\n"
        "0  1   2   1\n"
        "2  3   2   1\n";
    int ret;
    tsk_treeseq_t ts;
    const tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
        2, 2, z,
        z, z, z,
        2, 2, z,
    };
    uint32_t num_trees = 3;
    tsk_tree_t tree;

    tsk_treeseq_from_text(&ts, 3, nodes_txt, edges_txt, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 3.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 3);

    verify_trees(&ts, num_trees, parents);

    ret = tsk_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tree.left_root, TSK_NULL);
    CU_ASSERT_EQUAL(tree.parent[0], 2);
    CU_ASSERT_EQUAL(tree.parent[1], 2);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 0);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tree.left_root, TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 0);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tree.left_root, TSK_NULL);
    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 0);
    CU_ASSERT_EQUAL(tree.parent[0], 2);
    CU_ASSERT_EQUAL(tree.parent[1], 2);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_initial_gap_tsk_treeseq_mutation_parents(void)
{
    const char *nodes_txt =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0";
    const char *edges_txt =
        "2  3   2   0,1\n";
    const char *sites_txt =
        "0.5  0\n"
        "1.5  0\n"
        "2.5  0\n";
    const char *mutations_txt =
        "0    0     1\n"
        "0    0     1\n"
        "1    1     1\n"
        "1    1     1\n"
        "2    2     1\n"
        "2    2     1\n";
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    int ret;

    tsk_treeseq_from_text(&ts, 3, nodes_txt, edges_txt, NULL, sites_txt,
            mutations_txt, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 2);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.mutations->parent[0], -1);
    CU_ASSERT_EQUAL(tables.mutations->parent[1], 0);
    CU_ASSERT_EQUAL(tables.mutations->parent[2], -1);
    CU_ASSERT_EQUAL(tables.mutations->parent[3], 2);
    CU_ASSERT_EQUAL(tables.mutations->parent[4], -1);
    CU_ASSERT_EQUAL(tables.mutations->parent[5], 4);
    tsk_tbl_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_final_gap_tree_sequence(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0";
    const char *edges =
        "0  2   2   0,1\n";
    const char *sites =
        "0.5  0\n"
        "1.5  0\n"
        "2.5  0\n";
    const char *mutations =
        "0    0     1\n"
        "1    1     1\n"
        "2    0     1";
    const char *haplotypes[] = {"101", "010"};
    char *haplotype;
    unsigned int j;
    int ret;
    tsk_treeseq_t ts;
    tsk_hapgen_t hapgen;
    const tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
        2, 2, z,
        z, z, z,
    };
    uint32_t num_trees = 2;

    tsk_treeseq_from_text(&ts, 3, nodes, edges, NULL, sites, mutations, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 3.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 2);

    verify_trees(&ts, num_trees, parents);

    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < 2; j++) {
        ret = tsk_hapgen_get_haplotype(&hapgen, (tsk_id_t) j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, haplotypes[j]);
    }
    tsk_hapgen_free(&hapgen);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_final_gap_tsk_treeseq_mutation_parents(void)
{
    const char *nodes_txt =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0";
    const char *edges_txt =
        "0  2   2   0,1\n";
    const char *sites_txt =
        "0.5  0\n"
        "1.5  0\n"
        "2.5  0\n";
    const char *mutations_txt =
        "0    0     1\n"
        "0    0     1\n"
        "1    1     1\n"
        "1    1     1\n"
        "2    0     1\n"
        "2    0     1\n";
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    int ret;

    tsk_treeseq_from_text(&ts, 3, nodes_txt, edges_txt, NULL, sites_txt,
            mutations_txt, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 2);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.mutations->parent[0], -1);
    CU_ASSERT_EQUAL(tables.mutations->parent[1], 0);
    CU_ASSERT_EQUAL(tables.mutations->parent[2], -1);
    CU_ASSERT_EQUAL(tables.mutations->parent[3], 2);
    CU_ASSERT_EQUAL(tables.mutations->parent[4], -1);
    CU_ASSERT_EQUAL(tables.mutations->parent[5], 4);
    tsk_tbl_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_individuals(void)
{
    const char *individuals =
        "1      0.2\n"
        "2      0.5,0.6\n";
    const char *nodes =
        "1  0   -1  -1\n"
        "1  0   -1  1\n"
        "0  0   -1  -1\n"
        "1  0   -1  0\n"
        "0  0   -1  1\n";
    tsk_tbl_collection_t tables;
    tsk_treeseq_t ts;
    tsk_node_t node;
    tsk_individual_t individual;
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1.0;
    parse_individuals(individuals, tables.individuals);
    CU_ASSERT_EQUAL_FATAL(tables.individuals->num_rows, 2);
    parse_nodes(nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 5);

    ret = tsk_treeseq_alloc(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_get_node(&ts, 0, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(node.individual, TSK_NULL);

    ret = tsk_treeseq_get_node(&ts, 1, &node);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(node.individual, 1);

    ret = tsk_treeseq_get_individual(&ts, 0, &individual);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(individual.id, 0);
    CU_ASSERT_EQUAL_FATAL(individual.flags, 1);
    CU_ASSERT_EQUAL_FATAL(individual.location_length, 1);
    CU_ASSERT_EQUAL_FATAL(individual.location[0], 0.2);
    CU_ASSERT_EQUAL_FATAL(individual.nodes_length, 1);
    CU_ASSERT_EQUAL_FATAL(individual.nodes[0], 3);

    ret = tsk_treeseq_get_individual(&ts, 1, &individual);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(individual.id, 1);
    CU_ASSERT_EQUAL_FATAL(individual.flags, 2);
    CU_ASSERT_EQUAL_FATAL(individual.location_length, 2);
    CU_ASSERT_EQUAL_FATAL(individual.location[0], 0.5);
    CU_ASSERT_EQUAL_FATAL(individual.location[1], 0.6);
    CU_ASSERT_EQUAL_FATAL(individual.nodes_length, 2);
    CU_ASSERT_EQUAL_FATAL(individual.nodes[0], 1);
    CU_ASSERT_EQUAL_FATAL(individual.nodes[1], 4);

    ret = tsk_treeseq_get_individual(&ts, 3, &individual);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);

    tsk_tbl_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_simplest_bad_individuals(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "1  0   0\n"
        "0  1   0\n";
    const char *edges =
        "0  1   2   0\n"
        "0  1   2   1\n"
        "0  1   4   3\n";
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    int load_flags = TSK_BUILD_INDEXES;
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1.0;
    parse_nodes(nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 5);
    parse_edges(edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 3);
    ret = tsk_population_tbl_add_row(tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Make sure we have a good set of records */
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    /* Bad individual ID */
    tables.nodes->individual[0] = -2;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes->individual[0] = TSK_NULL;

    /* Bad individual ID */
    tables.nodes->individual[0] = 0;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes->individual[0] = TSK_NULL;

    /* add two individuals */
    ret = tsk_individual_tbl_add_row(tables.individuals, 0, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_individual_tbl_add_row(tables.individuals, 0, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL(ret, 1);

    /* Bad individual ID */
    tables.nodes->individual[0] = 2;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes->individual[0] = TSK_NULL;

    tsk_treeseq_free(&ts);
    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_bad_edges(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "1  0   0\n"
        "0  1   0\n";
    const char *edges =
        "0  1   2   0\n"
        "0  1   2   1\n"
        "0  1   4   3\n";
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    int ret;
    int load_flags = TSK_BUILD_INDEXES;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1.0;
    parse_nodes(nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 5);
    parse_edges(edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 3);
    ret = tsk_population_tbl_add_row(tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Make sure we have a good set of records */
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts);

    /* NULL for tables should be an error */
    ret = tsk_treeseq_alloc(&ts, NULL, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    tsk_treeseq_free(&ts);

    /* Bad population ID */
    tables.nodes->population[0] = -2;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes->population[0] = 0;

    /* Bad population ID */
    tables.nodes->population[0] = 1;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes->population[0] = 0;

    /* Bad interval */
    tables.edges->right[0] = 0.0;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_EDGE_INTERVAL);
    tsk_treeseq_free(&ts);
    tables.edges->right[0]= 1.0;

    /* Left coordinate < 0. */
    tables.edges->left[0] = -1;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_LEFT_LESS_ZERO);
    tsk_treeseq_free(&ts);
    tables.edges->left[0]= 0.0;

    /* Right coordinate > sequence length. */
    tables.edges->right[0] = 2.0;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_RIGHT_GREATER_SEQ_LENGTH);
    tsk_treeseq_free(&ts);
    tables.edges->right[0]= 1.0;

    /* Duplicate records */
    tables.edges->child[0] = 1;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_DUPLICATE_EDGES);
    tsk_treeseq_free(&ts);
    tables.edges->child[0] = 0;

    /* Duplicate records */
    tables.edges->child[0] = 1;
    tables.edges->left[0] = 0.5;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_EDGES_NOT_SORTED_LEFT);
    tsk_treeseq_free(&ts);
    tables.edges->child[0] = 0;
    tables.edges->left[0] = 0.0;

    /* child node == parent */
    tables.edges->child[1] = 2;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_NODE_TIME_ORDERING);
    tsk_treeseq_free(&ts);
    tables.edges->child[1] = 1;

    /* Unsorted child nodes */
    tables.edges->child[0] = 1;
    tables.edges->child[1] = 0;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_EDGES_NOT_SORTED_CHILD);
    tsk_treeseq_free(&ts);
    tables.edges->child[0] = 0;
    tables.edges->child[1] = 1;

    /* discontinuous parent nodes */
    /* Swap rows 1 and 2 */
    tables.edges->parent[1] = 4;
    tables.edges->child[1] = 3;
    tables.edges->parent[2] = 2;
    tables.edges->child[2] = 1;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_EDGES_NONCONTIGUOUS_PARENTS);
    tsk_treeseq_free(&ts);
    tables.edges->parent[2] = 4;
    tables.edges->child[2] = 3;
    tables.edges->parent[1] = 2;
    tables.edges->child[1] = 1;

    /* Null parent */
    tables.edges->parent[0] = TSK_NULL;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NULL_PARENT);
    tsk_treeseq_free(&ts);
    tables.edges->parent[0] = 2;

    /* parent not in nodes list */
    tables.nodes->num_rows = 2;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.nodes->num_rows = 5;

    /* parent negative */
    tables.edges->parent[0] = -2;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.edges->parent[0] = 2;

    /* Null child */
    tables.edges->child[0] = TSK_NULL;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NULL_CHILD);
    tsk_treeseq_free(&ts);
    tables.edges->child[0] = 0;

    /* child node reference out of bounds */
    tables.edges->child[0] = 100;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.edges->child[0] = 0;

    /* child node reference negative */
    tables.edges->child[0] = -2;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.edges->child[0] = 0;

    /* Make sure we've preserved a good tree sequence */
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    tsk_treeseq_free(&ts);

    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_bad_indexes(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "1  0   0\n"
        "0  1   0\n";
    const char *edges =
        "0  1   2   0\n"
        "0  1   2   1\n"
        "0  1   4   3\n";
    tsk_tbl_collection_t tables;
    tsk_id_t bad_indexes[] = {-1, 3, 4, 1000};
    size_t j;
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1.0;
    parse_nodes(nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 5);
    parse_edges(edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 3);
    ret = tsk_population_tbl_add_row(tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Make sure we have a good set of records */
    ret = tsk_tbl_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_tbl_collection_check_integrity(&tables, TSK_CHECK_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLES_NOT_INDEXED);
    ret = tsk_tbl_collection_build_indexes(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_check_integrity(&tables, TSK_CHECK_ALL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < sizeof(bad_indexes) / sizeof(*bad_indexes); j++) {
        tables.indexes.edge_insertion_order[0] = bad_indexes[j];
        ret = tsk_tbl_collection_check_integrity(&tables, TSK_CHECK_ALL);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
        tables.indexes.edge_insertion_order[0] = 0;

        tables.indexes.edge_removal_order[0] = bad_indexes[j];
        ret = tsk_tbl_collection_check_integrity(&tables, TSK_CHECK_ALL);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
        tables.indexes.edge_removal_order[0] = 0;
    }

    ret = tsk_tbl_collection_drop_indexes(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_check_integrity(&tables, TSK_CHECK_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLES_NOT_INDEXED);

    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_bad_migrations(void)
{
    tsk_tbl_collection_t tables;
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    /* insert two populations and one node to refer to. */
    ret = tsk_node_tbl_add_row(tables.nodes, 0, 0.0, TSK_NULL,
            TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_tbl_add_row(tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_tbl_add_row(tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    /* One migration, node 0 goes from population 0 to 1. */
    ret = tsk_migration_tbl_add_row(tables.migrations, 0, 1, 0, 0, 1, 1.0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* We only need basic intregity checks for migrations */
    ret = tsk_tbl_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Bad node reference */
    tables.migrations->node[0] = -1;
    ret = tsk_tbl_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tables.migrations->node[0] = 0;

    /* Bad node reference */
    tables.migrations->node[0] = 1;
    ret = tsk_tbl_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tables.migrations->node[0] = 0;

    /* Bad population reference */
    tables.migrations->source[0] = -1;
    ret = tsk_tbl_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tables.migrations->source[0] = 0;

    /* Bad population reference */
    tables.migrations->source[0] = 2;
    ret = tsk_tbl_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tables.migrations->source[0] = 0;

    /* Bad population reference */
    tables.migrations->dest[0] = -1;
    ret = tsk_tbl_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tables.migrations->dest[0] = 1;

    /* Bad population reference */
    tables.migrations->dest[0] = 2;
    ret = tsk_tbl_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tables.migrations->dest[0] = 1;

    /* Bad left coordinate */
    tables.migrations->left[0] = -1;
    ret = tsk_tbl_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_LEFT_LESS_ZERO);
    tables.migrations->left[0] = 0;

    /* Bad right coordinate */
    tables.migrations->right[0] = 2;
    ret = tsk_tbl_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_RIGHT_GREATER_SEQ_LENGTH);
    tables.migrations->right[0] = 1;

    /* Bad interval coordinate */
    tables.migrations->right[0] = 0;
    ret = tsk_tbl_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_EDGE_INTERVAL);
    tables.migrations->right[0] = 1;

    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_migration_simplify(void)
{
    tsk_tbl_collection_t tables;
    int ret;
    tsk_id_t samples[] = {0, 1};

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    /* insert two populations and one node to refer to. */
    ret = tsk_node_tbl_add_row(tables.nodes, TSK_NODE_IS_SAMPLE, 0.0,
            TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_tbl_add_row(tables.nodes, TSK_NODE_IS_SAMPLE, 0.0,
            TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    ret = tsk_population_tbl_add_row(tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_tbl_add_row(tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    /* One migration, node 0 goes from population 0 to 1. */
    ret = tsk_migration_tbl_add_row(tables.migrations, 0, 1, 0, 0, 1, 1.0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_tbl_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SIMPLIFY_MIGRATIONS_NOT_SUPPORTED);

    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_overlapping_parents(void)
{
    const char *nodes =
        "1  0   -1\n"
        "1  0   -1\n"
        "0  1   -1\n";
    const char *edges =
        "0  1   2   0\n"
        "0  1   2   1\n";
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    tsk_tree_t tree;
    int ret;
    int load_flags = TSK_BUILD_INDEXES;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 3);
    parse_edges(edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 2);

    tables.edges->left[0] = 0;
    tables.edges->parent[0] = 2;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tree.parent[0], 2);
    CU_ASSERT_EQUAL(tree.parent[1], 2);
    CU_ASSERT_EQUAL(tree.left_sib[2], TSK_NULL);
    CU_ASSERT_EQUAL(tree.right_sib[2], TSK_NULL);
    CU_ASSERT_EQUAL(tree.left_child[2], 0);
    CU_ASSERT_EQUAL(tree.right_child[2], 1);
    CU_ASSERT_EQUAL(tree.left_sib[0], TSK_NULL);
    CU_ASSERT_EQUAL(tree.right_sib[0], 1);
    CU_ASSERT_EQUAL(tree.left_sib[1], 0);
    CU_ASSERT_EQUAL(tree.right_sib[1], TSK_NULL);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_contradictory_children(void)
{
    const char *nodes =
        "1  0   -1\n"
        "1  1   -1\n"
        "0  1   -1\n";
    const char *edges =
        "0  1   1   0\n"
        "0  1   2   0\n";
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    tsk_tree_t tree;
    int ret;
    int load_flags = TSK_BUILD_INDEXES;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    parse_nodes(nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 3);
    parse_edges(edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 2);
    tables.sequence_length = 1.0;

    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_EDGES_CONTRADICTORY_CHILDREN);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_overlapping_edges_simplify(void)
{
    const char *nodes =
        "1  0   -1\n"
        "1  0   -1\n"
        "1  0   -1\n"
        "0  1   -1";
    const char *edges =
        "0  2   3   0\n"
        "1  3   3   1\n"
        "0  3   3   2\n";
    tsk_id_t samples[] = {0, 1, 2};
    tsk_tbl_collection_t tables;
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 3;
    parse_nodes(nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 4);
    parse_edges(edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 3);

    ret = tsk_tbl_collection_simplify(&tables, samples, 3, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(tables.nodes->num_rows, 4);
    CU_ASSERT_EQUAL(tables.edges->num_rows, 3);

    /* Identical to the input.
    0  2   3   0
    1  3   3   1
    0  3   3   2
    */
    CU_ASSERT_EQUAL(tables.edges->left[0], 0);
    CU_ASSERT_EQUAL(tables.edges->left[1], 1);
    CU_ASSERT_EQUAL(tables.edges->left[2], 0);
    CU_ASSERT_EQUAL(tables.edges->right[0], 2);
    CU_ASSERT_EQUAL(tables.edges->right[1], 3);
    CU_ASSERT_EQUAL(tables.edges->right[2], 3);
    CU_ASSERT_EQUAL(tables.edges->parent[0], 3);
    CU_ASSERT_EQUAL(tables.edges->parent[1], 3);
    CU_ASSERT_EQUAL(tables.edges->parent[2], 3);
    CU_ASSERT_EQUAL(tables.edges->child[0], 0);
    CU_ASSERT_EQUAL(tables.edges->child[1], 1);
    CU_ASSERT_EQUAL(tables.edges->child[2], 2);

    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_overlapping_unary_edges_simplify(void)
{
    const char *nodes =
        "1  0   -1\n"
        "1  0   -1\n"
        "0  1   -1";
    const char *edges =
        "0  2   2   0\n"
        "1  3   2   1\n";
    tsk_id_t samples[] = {0, 1};
    tsk_tbl_collection_t tables;
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 3;
    parse_nodes(nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 3);
    parse_edges(edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 2);

    ret = tsk_tbl_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(tables.nodes->num_rows, 3);
    CU_ASSERT_EQUAL(tables.edges->num_rows, 2);

    /* Because we only sample 0 and 1, the flanking unary edges are removed
     1       2       2       0
     1       2       2       1
     */
    CU_ASSERT_EQUAL(tables.edges->left[0], 1);
    CU_ASSERT_EQUAL(tables.edges->right[0], 2);
    CU_ASSERT_EQUAL(tables.edges->parent[0], 2);
    CU_ASSERT_EQUAL(tables.edges->child[0], 0);
    CU_ASSERT_EQUAL(tables.edges->left[1], 1);
    CU_ASSERT_EQUAL(tables.edges->right[1], 2);
    CU_ASSERT_EQUAL(tables.edges->parent[1], 2);
    CU_ASSERT_EQUAL(tables.edges->child[1], 1);

    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_overlapping_unary_edges_internal_samples_simplify(void)
{
    const char *nodes =
        "1  0   -1\n"
        "1  0   -1\n"
        "1  1   -1";
    const char *edges =
        "0  2   2   0\n"
        "1  3   2   1\n";
    tsk_id_t samples[] = {0, 1, 2};
    tsk_tbl_collection_t tables;
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 3;
    parse_nodes(nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 3);
    parse_edges(edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 2);

    ret = tsk_tbl_collection_simplify(&tables, samples, 3, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(tables.nodes->num_rows, 3);
    CU_ASSERT_EQUAL(tables.edges->num_rows, 2);
    /* Identical to the input.
        0  2   2   0
        1  3   2   1
     */
    CU_ASSERT_EQUAL(tables.edges->left[0], 0);
    CU_ASSERT_EQUAL(tables.edges->left[1], 1);
    CU_ASSERT_EQUAL(tables.edges->right[0], 2);
    CU_ASSERT_EQUAL(tables.edges->right[1], 3);
    CU_ASSERT_EQUAL(tables.edges->parent[0], 2);
    CU_ASSERT_EQUAL(tables.edges->parent[1], 2);
    CU_ASSERT_EQUAL(tables.edges->child[0], 0);
    CU_ASSERT_EQUAL(tables.edges->child[1], 1);

    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_reduce_site_topology(void)
{
    /* Two trees side by side, with a site on the second one. The first
     * tree should disappear. */
    const char *nodes =
        "1  0   -1\n"
        "1  0   -1\n"
        "0  1   -1\n"
        "0  2   -1\n";
    const char *edges =
        "0  1   2   0\n"
        "0  1   2   1\n"
        "1  2   3   0\n"
        "1  2   3   1\n";
    const char *sites =
        "1.0  0\n";
    tsk_id_t samples[] = {0, 1};
    tsk_tbl_collection_t tables;
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 2;
    parse_nodes(nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 4);
    parse_edges(edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 4);
    parse_sites(sites, tables.sites);
    CU_ASSERT_EQUAL_FATAL(tables.sites->num_rows, 1);

    ret = tsk_tbl_collection_simplify(&tables, samples, 2,
            TSK_REDUCE_TO_SITE_TOPOLOGY, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(tables.nodes->num_rows, 3);
    CU_ASSERT_EQUAL(tables.edges->num_rows, 2);
    CU_ASSERT_EQUAL(tables.edges->left[0], 0);
    CU_ASSERT_EQUAL(tables.edges->left[1], 0);
    CU_ASSERT_EQUAL(tables.edges->right[0], 2);
    CU_ASSERT_EQUAL(tables.edges->right[1], 2);
    CU_ASSERT_EQUAL(tables.edges->parent[0], 2);
    CU_ASSERT_EQUAL(tables.edges->parent[1], 2);
    CU_ASSERT_EQUAL(tables.edges->child[0], 0);
    CU_ASSERT_EQUAL(tables.edges->child[1], 1);

    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_population_filter(void)
{
    tsk_tbl_collection_t tables;
    tsk_id_t samples[] = {0, 1};
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    tsk_population_tbl_add_row(tables.populations, "0", 1);
    tsk_population_tbl_add_row(tables.populations, "1", 1);
    tsk_population_tbl_add_row(tables.populations, "2", 1);
    /* Two nodes referring to population 1 */
    tsk_node_tbl_add_row(tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 1, TSK_NULL,
            NULL, 0);
    tsk_node_tbl_add_row(tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 1, TSK_NULL,
            NULL, 0);

    ret = tsk_tbl_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes->num_rows, 2);
    CU_ASSERT_EQUAL(tables.populations->num_rows, 3);
    CU_ASSERT_EQUAL(tables.populations->metadata[0], '0');
    CU_ASSERT_EQUAL(tables.populations->metadata[1], '1');
    CU_ASSERT_EQUAL(tables.populations->metadata[2], '2');

    ret = tsk_tbl_collection_simplify(&tables, samples, 2, TSK_FILTER_POPULATIONS, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes->num_rows, 2);
    CU_ASSERT_EQUAL(tables.nodes->population[0], 0);
    CU_ASSERT_EQUAL(tables.nodes->population[1], 0);
    CU_ASSERT_EQUAL(tables.populations->num_rows, 1);
    CU_ASSERT_EQUAL(tables.populations->metadata[0], '1');

    tsk_tbl_collection_free(&tables);
}

static void
test_simplest_individual_filter(void)
{
    tsk_tbl_collection_t tables;
    tsk_id_t samples[] = {0, 1};
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    tsk_individual_tbl_add_row(tables.individuals, 0, NULL, 0, "0", 1);
    tsk_individual_tbl_add_row(tables.individuals, 0, NULL, 0, "1", 1);
    tsk_individual_tbl_add_row(tables.individuals, 0, NULL, 0, "2", 1);
    /* Two nodes referring to individual 1 */
    tsk_node_tbl_add_row(tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, 1,
            NULL, 0);
    tsk_node_tbl_add_row(tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, 1,
            NULL, 0);

    ret = tsk_tbl_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes->num_rows, 2);
    CU_ASSERT_EQUAL(tables.individuals->num_rows, 3);
    CU_ASSERT_EQUAL(tables.individuals->metadata[0], '0');
    CU_ASSERT_EQUAL(tables.individuals->metadata[1], '1');
    CU_ASSERT_EQUAL(tables.individuals->metadata[2], '2');

    ret = tsk_tbl_collection_simplify(&tables, samples, 2, TSK_FILTER_INDIVIDUALS, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes->num_rows, 2);
    CU_ASSERT_EQUAL(tables.nodes->individual[0], 0);
    CU_ASSERT_EQUAL(tables.nodes->individual[1], 0);
    CU_ASSERT_EQUAL(tables.individuals->num_rows, 1);
    CU_ASSERT_EQUAL(tables.individuals->metadata[0], '1');

    tsk_tbl_collection_free(&tables);
}

/*=======================================================
 * Single tree tests.
 *======================================================*/

static void
test_single_tree_good_records(void)
{
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges,
            NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 7);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_treeseq_free(&ts);
}


static void
test_single_nonbinary_tree_good_records(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  2   0\n"
        "0  3   0\n";
    const char *edges =
        "0 1 7 0,1,2,3\n"
        "0 1 8 4,5\n"
        "0 1 9 6,7,8";
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 7);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 10);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_bad_records(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    int load_flags = TSK_BUILD_INDEXES;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(single_tree_ex_nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 7);
    parse_edges(single_tree_ex_edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 6);

    /* Not sorted in time order */
    tables.nodes->time[5] = 0.5;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_EDGES_NOT_SORTED_PARENT_TIME);
    tsk_treeseq_free(&ts);
    tables.nodes->time[5] = 2.0;

    /* Left value greater than sequence right */
    tables.edges->left[2] = 2.0;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_EDGE_INTERVAL);
    tsk_treeseq_free(&ts);
    tables.edges->left[2] = 0.0;

    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    tsk_treeseq_free(&ts);
    tsk_tbl_collection_free(&tables);
}


static void
test_single_tree_good_mutations(void)
{
    tsk_treeseq_t ts;
    size_t j;
    size_t num_sites = 3;
    size_t num_mutations = 7;
    tsk_site_t other_sites[num_sites];
    tsk_mutation_t other_mutations[num_mutations];
    int ret;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges,
            NULL, single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 7);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), num_sites);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), num_mutations);

    for (j = 0; j < num_sites; j++) {
        ret = tsk_treeseq_get_site(&ts, j, other_sites + j);
        CU_ASSERT_EQUAL(ret, 0);
    }
    for (j = 0; j < num_mutations; j++) {
        ret = tsk_treeseq_get_mutation(&ts, j, other_mutations + j);
        CU_ASSERT_EQUAL(ret, 0);
    }
    CU_ASSERT_EQUAL(other_sites[0].position, 0.1);
    CU_ASSERT_NSTRING_EQUAL(other_sites[0].ancestral_state, "0", 1);
    CU_ASSERT_EQUAL(other_sites[1].position, 0.2);
    CU_ASSERT_NSTRING_EQUAL(other_sites[1].ancestral_state, "0", 1);
    CU_ASSERT_EQUAL(other_sites[2].position, 0.3);
    CU_ASSERT_NSTRING_EQUAL(other_sites[2].ancestral_state, "0", 1);

    CU_ASSERT_EQUAL(other_mutations[0].id, 0);
    CU_ASSERT_EQUAL(other_mutations[0].node, 2);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[0].derived_state, "1", 1);
    CU_ASSERT_EQUAL(other_mutations[1].id, 1);
    CU_ASSERT_EQUAL(other_mutations[1].node, 4);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[1].derived_state, "1", 1);
    CU_ASSERT_EQUAL(other_mutations[2].id, 2);
    CU_ASSERT_EQUAL(other_mutations[2].node, 0);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[2].derived_state, "0", 1);
    CU_ASSERT_EQUAL(other_mutations[3].id, 3);
    CU_ASSERT_EQUAL(other_mutations[3].node, 0);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[3].derived_state, "1", 1);
    CU_ASSERT_EQUAL(other_mutations[4].id, 4);
    CU_ASSERT_EQUAL(other_mutations[4].node, 1);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[4].derived_state, "1", 1);
    CU_ASSERT_EQUAL(other_mutations[5].id, 5);
    CU_ASSERT_EQUAL(other_mutations[5].node, 2);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[5].derived_state, "1", 1);
    CU_ASSERT_EQUAL(other_mutations[6].id, 6);
    CU_ASSERT_EQUAL(other_mutations[6].node, 3);
    CU_ASSERT_NSTRING_EQUAL(other_mutations[6].derived_state, "1", 1);

    tsk_treeseq_free(&ts);
}

static void
test_single_tree_bad_mutations(void)
{
    int ret = 0;
    const char *sites =
        "0       0\n"
        "0.1     0\n"
        "0.2     0\n";
    const char *mutations =
        "0   0  1  -1\n"
        "1   1  1  -1\n"
        "2   4  1  -1\n"
        "2   1  0  2\n"
        "2   1  1  3\n"
        "2   2  1  -1\n";
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    int load_flags = TSK_BUILD_INDEXES;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(single_tree_ex_nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 7);
    parse_edges(single_tree_ex_edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 6);
    parse_sites(sites, tables.sites);
    parse_mutations(mutations, tables.mutations);
    CU_ASSERT_EQUAL_FATAL(tables.sites->num_rows, 3);
    CU_ASSERT_EQUAL_FATAL(tables.mutations->num_rows, 6);
    tables.sequence_length = 1.0;

    /* Check to make sure we have legal mutations */
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    tsk_treeseq_free(&ts);

    /* negative coordinate */
    tables.sites->position[0] = -1.0;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SITE_POSITION);
    tsk_treeseq_free(&ts);
    tables.sites->position[0] = 0.0;

    /* coordinate == sequence length */
    tables.sites->position[2] = 1.0;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SITE_POSITION);
    tsk_treeseq_free(&ts);
    tables.sites->position[2] = 0.2;

    /* coordinate > sequence length */
    tables.sites->position[2] = 1.1;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SITE_POSITION);
    tsk_treeseq_free(&ts);
    tables.sites->position[2] = 0.2;

    /* Duplicate positions */
    tables.sites->position[0] = 0.1;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_DUPLICATE_SITE_POSITION);
    tsk_treeseq_free(&ts);
    tables.sites->position[0] = 0.0;

    /* Unsorted positions */
    tables.sites->position[0] = 0.3;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSORTED_SITES);
    tsk_treeseq_free(&ts);
    tables.sites->position[0] = 0.0;

    /* site < 0 */
    tables.mutations->site[0] = -2;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations->site[0] = 0;

    /* site == num_sites */
    tables.mutations->site[0] = 3;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations->site[0] = 0;

    /* node = NULL */
    tables.mutations->node[0] = TSK_NULL;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations->node[0] = 0;

    /* node >= num_nodes */
    tables.mutations->node[0] = 7;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations->node[0] = 0;

    /* parent < -1 */
    tables.mutations->parent[0] = -2;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations->parent[0] = TSK_NULL;

    /* parent >= num_mutations */
    tables.mutations->parent[0] = 7;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts);
    tables.mutations->parent[0] = TSK_NULL;

    /* parent on a different site */
    tables.mutations->parent[1] = 0;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_PARENT_DIFFERENT_SITE);
    tsk_treeseq_free(&ts);
    tables.mutations->parent[1] = TSK_NULL;

    /* parent is the same mutation */
    tables.mutations->parent[0] = 0;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_PARENT_EQUAL);
    tsk_treeseq_free(&ts);
    tables.mutations->parent[0] = TSK_NULL;

    /* parent_id > mutation id */
    tables.mutations->parent[2] = 3;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_PARENT_AFTER_CHILD);
    tsk_treeseq_free(&ts);
    tables.mutations->parent[2] = TSK_NULL;

    /* Check to make sure we've maintained legal mutations */
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    tsk_treeseq_free(&ts);

    tsk_tbl_collection_free(&tables);
}

static void
test_single_tree_iter(void)
{
    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  2   0\n"
        "0  3   0\n";
    const char *edges =
        "0  6   4   0,1\n"
        "0  6   5   2,3\n"
        "0  6   6   4,5\n";
    tsk_id_t parents[] = {4, 4, 5, 5, 6, 6, TSK_NULL};
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_id_t u, v, w;
    size_t num_samples;
    uint32_t num_nodes = 7;

    tsk_treeseq_from_text(&ts, 6, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    ret = tsk_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_tree_print_state(&tree, _devnull);

    for (u = 0; u < (tsk_id_t) num_nodes; u++) {
        ret = tsk_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
    }
    ret = tsk_tree_get_num_samples(&tree, 0, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 1);
    ret = tsk_tree_get_num_samples(&tree, 4, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 2);
    ret = tsk_tree_get_num_samples(&tree, 6, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 4);
    ret = tsk_tree_get_mrca(&tree, 0, 1, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 4);
    ret = tsk_tree_get_mrca(&tree, 0, 2, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 6);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_single_nonbinary_tree_iter(void)
{
    int ret;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  1   0\n"
        "0  2   0\n"
        "0  3   0\n";
    const char *edges =
        "0  1   7   0,1,2,3\n"
        "0  1   8   4,5\n"
        "0  1   9   6,7,8\n";
    tsk_id_t parents[] = {7, 7, 7, 7, 8, 8, 9, 9, 9, TSK_NULL};
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_id_t u, v, w;
    size_t num_samples;
    size_t num_nodes = 10;
    size_t total_samples = 7;

    tsk_treeseq_from_text(&ts, 1, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    ret = tsk_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_tree_print_state(&tree, _devnull);

    for (u = 0; u < (tsk_id_t) num_nodes; u++) {
        ret = tsk_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
    }
    for (u = 0; u < (tsk_id_t) total_samples; u++) {
        ret = tsk_tree_get_num_samples(&tree, u, &num_samples);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(num_samples, 1);
        CU_ASSERT_EQUAL(tree.left_child[u], TSK_NULL);
    }

    u = 7;
    ret = tsk_tree_get_num_samples(&tree, u, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 4);
    CU_ASSERT_EQUAL(tree.right_child[u], 3);
    CU_ASSERT_EQUAL(tree.left_sib[3], 2);
    CU_ASSERT_EQUAL(tree.left_sib[2], 1);
    CU_ASSERT_EQUAL(tree.left_sib[1], 0);
    CU_ASSERT_EQUAL(tree.left_sib[0], TSK_NULL);

    u = 8;
    ret = tsk_tree_get_num_samples(&tree, u, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 2);
    CU_ASSERT_EQUAL(tree.right_child[u], 5);
    CU_ASSERT_EQUAL(tree.left_sib[5], 4);
    CU_ASSERT_EQUAL(tree.left_sib[4], TSK_NULL);

    u = 9;
    ret = tsk_tree_get_num_samples(&tree, u, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 7);
    CU_ASSERT_EQUAL(tree.right_child[u], 8);
    CU_ASSERT_EQUAL(tree.left_sib[8], 7);
    CU_ASSERT_EQUAL(tree.left_sib[7], 6);
    CU_ASSERT_EQUAL(tree.left_sib[6], TSK_NULL);

    CU_ASSERT_EQUAL(tsk_tree_get_num_roots(&tree), 1);
    CU_ASSERT_EQUAL(tree.left_root, 9);

    ret = tsk_tree_get_mrca(&tree, 0, 1, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 7);
    ret = tsk_tree_get_mrca(&tree, 0, 4, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 9);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_general_samples_iter(void)
{
    int ret;
    const char *nodes =
        "0  3   0\n"
        "0  2   0\n"
        "0  1   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n";
    const char *edges =
        "0  6   2   3,4\n"
        "0  6   1   5,6\n"
        "0  6   0   1,2\n";
    tsk_id_t parents[] = {TSK_NULL, 0, 0, 2, 2, 1, 1};
    tsk_id_t *samples;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_id_t u, v, w;
    size_t num_samples;
    uint32_t num_nodes = 7;

    tsk_treeseq_from_text(&ts, 6, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    ret = tsk_treeseq_get_samples(&ts, &samples);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(samples[0], 3);
    CU_ASSERT_EQUAL(samples[1], 4);
    CU_ASSERT_EQUAL(samples[2], 5);
    CU_ASSERT_EQUAL(samples[3], 6);

    ret = tsk_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_tree_print_state(&tree, _devnull);

    for (u = 0; u < (tsk_id_t) num_nodes; u++) {
        ret = tsk_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
    }
    ret = tsk_tree_get_num_samples(&tree, 3, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 1);
    ret = tsk_tree_get_num_samples(&tree, 2, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 2);
    ret = tsk_tree_get_num_samples(&tree, 0, &num_samples);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_samples, 4);
    ret = tsk_tree_get_mrca(&tree, 3, 4, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 2);
    ret = tsk_tree_get_mrca(&tree, 3, 6, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 0);

    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

static void
test_single_tree_iter_times(void)
{
    int ret = 0;
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  2   0\n"
        "1  3   0\n"
        "0  1   0\n"
        "0  4   0\n"
        "0  5   0\n";
    const char *edges =
        "0  6   4   0,1\n"
        "0  6   5   2,3\n"
        "0  6   6   4,5\n";
    tsk_id_t parents[] = {4, 4, 5, 5, 6, 6, TSK_NULL};
    double times[] = {0.0, 0.0, 2.0, 3.0, 1.0, 4.0, 5.0};
    double t;
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    tsk_id_t u, v;
    uint32_t num_nodes = 7;

    tsk_treeseq_from_text(&ts, 6, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    ret = tsk_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_tree_print_state(&tree, _devnull);

    for (u = 0; u < (tsk_id_t) num_nodes; u++) {
        ret = tsk_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
        ret = tsk_tree_get_time(&tree, u, &t);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(t, times[u]);
    }
    ret = tsk_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}


static void
test_single_tree_simplify(void)
{
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    tsk_id_t samples[] = {0, 1};
    int ret;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
            single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL);
    verify_simplify(&ts);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_tbl_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tables.nodes->num_rows, 3);
    CU_ASSERT_EQUAL(tables.edges->num_rows, 2);

    /* Make sure we detect unsorted edges */
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    unsort_edges(tables.edges, 0);
    ret = tsk_tbl_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGES_NOT_SORTED_CHILD);

    /* detect bad parents */
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.edges->parent[0] = -1;
    ret = tsk_tbl_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NULL_PARENT);

    /* detect bad children */
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.edges->child[0] = -1;
    ret = tsk_tbl_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NULL_CHILD);

    /* detect loops */
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.edges->child[0] = tables.edges->parent[0];
    ret = tsk_tbl_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_NODE_TIME_ORDERING);

    /* Test the interface for NULL inputs */
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_simplify(&tables, NULL, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    tsk_treeseq_free(&ts);
    tsk_tbl_collection_free(&tables);
}

static void
test_single_tree_compute_mutation_parents(void)
{
    int ret = 0;
    const char *sites =
        "0       0\n"
        "0.1     0\n"
        "0.2     0\n";
    const char *mutations =
        "0   0  1  -1\n"
        "1   1  1  -1\n"
        "2   4  1  -1\n"
        "2   1  0  2\n"
        "2   1  1  3\n"
        "2   2  1  -1\n";
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(single_tree_ex_nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 7);
    parse_edges(single_tree_ex_edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 6);
    parse_sites(sites, tables.sites);
    parse_mutations(mutations, tables.mutations);
    CU_ASSERT_EQUAL_FATAL(tables.sites->num_rows, 3);
    CU_ASSERT_EQUAL_FATAL(tables.mutations->num_rows, 6);
    tables.sequence_length = 1.0;

    ret = tsk_tbl_collection_build_indexes(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Check to make sure we have legal mutations */
    ret = tsk_treeseq_alloc(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);

    /* Compute the mutation parents */
    verify_compute_mutation_parents(&ts);

    /* Verify consistency of individuals */
    verify_individual_nodes(&ts);
    tsk_treeseq_free(&ts);

    /* Bad site reference */
    tables.mutations->site[0] = -1;
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tables.mutations->site[0] = 0;

    /* Bad site reference */
    tables.mutations->site[0] = -1;
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tables.mutations->site[0] = 0;

    /* mutation sites out of order */
    tables.mutations->site[0] = 2;
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSORTED_MUTATIONS);
    tables.mutations->site[0] = 0;

    /* sites out of order */
    tables.sites->position[0] = 0.11;
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSORTED_SITES);
    tables.sites->position[0] = 0;

    /* Bad node reference */
    tables.mutations->node[0] = -1;
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tables.mutations->node[0] = 0;

    /* Bad node reference */
    tables.mutations->node[0] = (tsk_id_t) tables.nodes->num_rows;
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tables.mutations->node[0] = 0;

    /* Mutations not ordered by tree */
    tables.mutations->node[2] = 1;
    tables.mutations->node[3] = 4;
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_MUTATION_PARENT_AFTER_CHILD);
    tables.mutations->node[2] = 4;
    tables.mutations->node[3] = 1;

    /* Need to reset the parent field here */
    memset(tables.mutations->parent, 0xff,
            tables.mutations->num_rows * sizeof(tsk_id_t));
    /* Mutations not ordered by site */
    tables.mutations->site[3] = 1;
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_MUTATIONS);
    tables.mutations->site[3] = 2;

    /* Check to make sure we still have legal mutations */
    ret = tsk_tbl_collection_compute_mutation_parents(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_alloc(&ts, &tables, 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 3);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 6);
    tsk_treeseq_free(&ts);

    tsk_treeseq_free(&ts);
    tsk_tbl_collection_free(&tables);
}


/*=======================================================
 * Multi tree tests.
 *======================================================*/

static void
test_simple_multi_tree(void)
{
    tsk_id_t parents[] = {
        6, 5, 8, 5, TSK_NULL, 6, 8, TSK_NULL, TSK_NULL,
        6, 5, 4, 4, 5, 6, TSK_NULL, TSK_NULL, TSK_NULL,
        7, 5, 4, 4, 5, 7, TSK_NULL, TSK_NULL, TSK_NULL,
    };
    uint32_t num_trees = 3;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL,
            paper_ex_sites, paper_ex_mutations, paper_ex_individuals, NULL);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);
}

static void
test_unary_multi_tree(void)
{
    tsk_id_t parents[] = {
        6, 5, 7, 5, TSK_NULL, 6, 8, 8, TSK_NULL,
        6, 5, 4, 4, 5, 6, 8, TSK_NULL, TSK_NULL,
        7, 5, 4, 4, 5, 7, TSK_NULL, TSK_NULL, TSK_NULL,
    };
    tsk_treeseq_t ts;
    uint32_t num_trees = 3;

    tsk_treeseq_from_text(&ts, 10, unary_ex_nodes, unary_ex_edges, NULL,
            unary_ex_sites, unary_ex_mutations, NULL, NULL);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);
}

static void
test_internal_sample_multi_tree(void)
{
    tsk_id_t parents[] = {
        7, 5, 4, 4, 5, 7, TSK_NULL, TSK_NULL, TSK_NULL,
        4, 5, 4, 8, 5, 8, TSK_NULL, TSK_NULL, TSK_NULL,
        6, 5, 4, 4, 5, 6, TSK_NULL, TSK_NULL, TSK_NULL,
    };
    tsk_treeseq_t ts;
    uint32_t num_trees = 3;

    tsk_treeseq_from_text(&ts, 10, internal_sample_ex_nodes, internal_sample_ex_edges, NULL,
            internal_sample_ex_sites, internal_sample_ex_mutations, NULL, NULL);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);
}

static void
test_internal_sample_simplified_multi_tree(void)
{
    int ret;
    tsk_treeseq_t ts, simplified;
    tsk_id_t samples[] = {2, 3, 5};
    tsk_id_t node_map[9];
    tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
    /*  0  1  2  3  4 */
        3, 3, z, 2, z,
        2, 4, 4, z, z,
        3, 3, z, 2, z,
    };
    uint32_t num_trees = 3;

    tsk_treeseq_from_text(&ts, 10, internal_sample_ex_nodes, internal_sample_ex_edges, NULL,
            internal_sample_ex_sites, internal_sample_ex_mutations, NULL, NULL);
    ret = tsk_treeseq_simplify(&ts, samples, 3, 0, &simplified, node_map);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(node_map[2], 0);
    CU_ASSERT_EQUAL(node_map[3], 1);
    CU_ASSERT_EQUAL(node_map[5], 2);

    verify_trees(&simplified, num_trees, parents);
    tsk_treeseq_free(&simplified);
    tsk_treeseq_free(&ts);
}

static void
test_nonbinary_multi_tree(void)
{
    /* We make one mutation for each tree */
    tsk_id_t parents[] = {
        8, 8, 8, 8, 10, 10, 9, 10, 9, 12, 12, TSK_NULL, TSK_NULL,
        8, 8, 8, 8, 10, 11, 9, 10, 9, 11, 12, 12, TSK_NULL,
    };

    tsk_treeseq_t ts;
    uint32_t num_trees = 2;

    tsk_treeseq_from_text(&ts, 100, nonbinary_ex_nodes, nonbinary_ex_edges, NULL,
            nonbinary_ex_sites, nonbinary_ex_mutations, NULL, NULL);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);
}

static void
test_left_to_right_multi_tree(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  0.090   0\n"
        "0  0.170   0\n"
        "0  0.253   0\n"
        "0  0.071   0\n"
        "0  0.202   0\n";
    const char *edges =
        "2 10 7 2,3\n"
        "0 2  4 1\n"
        "2 10 4 1\n"
        "0 2  4 3\n"
        "2 10 4 7\n"
        "0 7  5 0,4\n"
        "7 10 8 0,4\n"
        "0 2  6 2,5\n";
    const char *sites =
        "1      0\n"
        "4.5    0\n"
        "8.5    0\n";
    const char *mutations =
        "0    2    1\n"
        "1    0    1\n"
        "2    4    1\n";

    tsk_id_t parents[] = {
        5, 4, 6, 4, 5, 6, TSK_NULL, TSK_NULL, TSK_NULL,
        5, 4, 7, 7, 5, TSK_NULL, TSK_NULL, 4, TSK_NULL,
        8, 4, 7, 7, 8, TSK_NULL, TSK_NULL, 4, TSK_NULL,
    };
    tsk_treeseq_t ts;
    uint32_t num_trees = 3;

    tsk_treeseq_from_text(&ts, 10, nodes, edges, NULL, sites, mutations, NULL, NULL);
    verify_trees(&ts, num_trees, parents);
    verify_tree_next_prev(&ts);
    tsk_treeseq_free(&ts);
}

static void
test_gappy_multi_tree(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "1  0   0\n"
        "0  0.090   0\n"
        "0  0.170   0\n"
        "0  0.253   0\n"
        "0  0.071   0\n"
        "0  0.202   0\n";
    const char *edges =
        "2 7  7 2\n"
        "8 10 7 2\n"
        "2 7  7 3\n"
        "8 10 7 3\n"
        "1 2  4 1\n"
        "2 7  4 1\n"
        "8 10 4 1\n"
        "1 2  4 3\n"
        "2 7  4 7\n"
        "8 10 4 7\n"
        "1 7  5 0,4\n"
        "8 10 8 0,4\n"
        "1 2  6 2,5\n";
    tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
        z, z, z, z, z, z, z, z, z,
        5, 4, 6, 4, 5, 6, z, z, z,
        5, 4, 7, 7, 5, z, z, 4, z,
        z, z, z, z, z, z, z, z, z,
        8, 4, 7, 7, 8, z, z, 4, z,
        z, z, z, z, z, z, z, z, z,
    };
    tsk_treeseq_t ts;
    uint32_t num_trees = 6;

    tsk_treeseq_from_text(&ts, 12, nodes, edges, NULL, NULL, NULL, NULL, NULL);
    verify_trees(&ts, num_trees, parents);
    verify_tree_next_prev(&ts);
    tsk_treeseq_free(&ts);
}

static void
test_tsk_treeseq_bad_records(void)
{
    int ret = 0;
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    uint32_t num_trees = 3;
    tsk_id_t parents[] = {
        6, 5, 8, 5, TSK_NULL, 6, 8, TSK_NULL, TSK_NULL,
        6, 5, 4, 4, 5, 6, TSK_NULL, TSK_NULL, TSK_NULL,
        7, 5, 4, 4, 5, 7, TSK_NULL, TSK_NULL, TSK_NULL,
    };
    int load_flags = TSK_BUILD_INDEXES;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 10;
    parse_nodes(paper_ex_nodes, tables.nodes);
    parse_edges(paper_ex_edges, tables.edges);
    parse_individuals(paper_ex_individuals, tables.individuals);

    /* Make sure we have a good set of records */
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ts.num_trees, 3);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);

    /* Left value greater than right */
    tables.edges->left[0] = 10.0;
    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_EDGE_INTERVAL);
    tsk_treeseq_free(&ts);
    tables.edges->left[0] = 2.0;

    ret = tsk_treeseq_alloc(&ts, &tables, load_flags);
    CU_ASSERT_EQUAL(ret, 0);
    verify_trees(&ts, num_trees, parents);
    tsk_treeseq_free(&ts);

    tsk_tbl_collection_free(&tables);
}

/*=======================================================
 * Diff iter tests.
 *======================================================*/

static void
test_simple_diff_iter(void)
{
    int ret;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, NULL, NULL,
            paper_ex_individuals, NULL);
    verify_tree_diffs(&ts);

    ret = tsk_treeseq_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_nonbinary_diff_iter(void)
{
    int ret;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 100, nonbinary_ex_nodes, nonbinary_ex_edges, NULL,
            NULL, NULL, NULL, NULL);
    verify_tree_diffs(&ts);

    ret = tsk_treeseq_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_unary_diff_iter(void)
{
    int ret;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, unary_ex_nodes, unary_ex_edges, NULL,
            NULL, NULL, NULL, NULL);
    verify_tree_diffs(&ts);

    ret = tsk_treeseq_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_internal_sample_diff_iter(void)
{
    int ret;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, internal_sample_ex_nodes, internal_sample_ex_edges, NULL,
            NULL, NULL, NULL, NULL);
    verify_tree_diffs(&ts);

    ret = tsk_treeseq_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
}

/*=======================================================
 * Sample sets
 *======================================================*/

static void
test_simple_sample_sets(void)
{
    sample_count_test_t tests[] = {
        {0, 0, 1}, {0, 5, 2}, {0, 6, 3},
        {1, 4, 2}, {1, 5, 3}, {1, 6, 4}};
    uint32_t num_tests = 6;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges,
            NULL, NULL, NULL, paper_ex_individuals, NULL);
    verify_sample_counts(&ts, num_tests, tests);
    verify_sample_sets(&ts);

    tsk_treeseq_free(&ts);
}

static void
test_nonbinary_sample_sets(void)
{
    sample_count_test_t tests[] = {
        {0, 0, 1}, {0, 8, 4}, {0, 9, 5}, {0, 10, 3}, {0, 12, 8},
        {1, 5, 1}, {1, 8, 4}, {1, 9, 5}, {0, 10, 2}, {0, 11, 1}};
    uint32_t num_tests = 8;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 100, nonbinary_ex_nodes, nonbinary_ex_edges, NULL,
            NULL, NULL, NULL, NULL);
    verify_sample_counts(&ts, num_tests, tests);
    verify_sample_sets(&ts);

    tsk_treeseq_free(&ts);
}

static void
test_internal_sample_sample_sets(void)
{
    sample_count_test_t tests[] = {
        {0, 0, 1}, {0, 5, 4}, {0, 4, 2}, {0, 7, 5},
        {1, 4, 2}, {1, 5, 4}, {1, 8, 5},
        {2, 5, 4}, {2, 6, 5}};
    uint32_t num_tests = 9;
    tsk_treeseq_t ts;

    tsk_treeseq_from_text(&ts, 10, internal_sample_ex_nodes, internal_sample_ex_edges,
            NULL, NULL, NULL, NULL, NULL);
    verify_sample_counts(&ts, num_tests, tests);
    verify_sample_sets(&ts);

    tsk_treeseq_free(&ts);
}


/*=======================================================
 * Miscellaneous tests.
 *======================================================*/

static void
test_genealogical_nearest_neighbours_errors(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_id_t *reference_sets[2];
    tsk_id_t reference_set_0[4], reference_set_1[4];
    tsk_id_t focal[] = {0, 1, 2, 3};
    size_t reference_set_size[2];
    size_t num_focal = 4;
    double *A = malloc(2 * num_focal * sizeof(double));
    CU_ASSERT_FATAL(A != NULL);

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges,
            NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 4);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);

    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, 0, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, INT16_MAX, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Overlapping sample sets */
    reference_sets[0] = focal;
    reference_set_size[0] = 1;
    reference_sets[1] = focal;
    reference_set_size[1] = num_focal;
    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);

    /* bad values in the sample sets */
    reference_set_0[0] = 0;
    reference_set_0[1] = 1;
    reference_set_1[0] = 2;
    reference_set_1[1] = 3;
    reference_set_size[0] = 2;
    reference_set_size[1] = 2;
    reference_sets[0] = reference_set_0;
    reference_sets[1] = reference_set_1;
    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    reference_set_0[0] = -1;
    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    reference_set_0[0] = (tsk_id_t) tsk_treeseq_get_num_nodes(&ts);
    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    reference_set_0[0] = (tsk_id_t) tsk_treeseq_get_num_nodes(&ts) + 1;
    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    /* Duplicate values in the focal sets */
    reference_set_0[0] = 1;
    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);
    reference_set_0[0] = 3;
    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);

    /* Bad sample ID */
    reference_sets[0] = focal;
    reference_set_size[0] = 1;
    reference_sets[1] = focal + 1;
    reference_set_size[1] = num_focal - 1;
    focal[0] = -1;
    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    focal[0] = (tsk_id_t) tsk_treeseq_get_num_nodes(&ts);
    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    focal[0] = (tsk_id_t) tsk_treeseq_get_num_nodes(&ts) + 100;
    ret = tsk_treeseq_genealogical_nearest_neighbours(&ts,
        focal, num_focal, reference_sets, reference_set_size, 2, 0, A);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    tsk_treeseq_free(&ts);
    free(A);
}

static void
test_tree_errors(void)
{
    int ret;
    size_t j;
    tsk_id_t num_nodes = 9;
    tsk_id_t u;
    tsk_node_t node;
    tsk_treeseq_t ts, other_ts;
    tsk_tree_t t, other_t;
    tsk_id_t bad_nodes[] = {num_nodes, num_nodes + 1, -1};
    tsk_id_t tracked_samples[] = {0, 0, 0};

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, NULL, NULL,
            paper_ex_individuals, NULL);

    ret = tsk_tree_alloc(&t, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_tree_alloc(&t, &ts, TSK_SAMPLE_COUNTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 1);

    /* Out-of-bounds queries */
    for (j = 0; j < sizeof(bad_nodes) / sizeof(tsk_id_t); j++) {
        u = bad_nodes[j];
        ret = tsk_tree_get_parent(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_tree_get_time(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_tree_get_mrca(&t, u, 0, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_tree_get_mrca(&t, 0, u, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_tree_get_num_samples(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_tree_get_num_tracked_samples(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        /* Also check tree sequence methods */
        ret = tsk_treeseq_get_node(&ts, (tsk_tbl_size_t) u, &node);
        CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        CU_ASSERT(!tsk_treeseq_is_sample(&ts, u));
        CU_ASSERT(!tsk_tree_is_sample(&t, u));
    }

    tracked_samples[0] = 0;
    tracked_samples[1] = (tsk_id_t) tsk_treeseq_get_num_samples(&ts);
    ret = tsk_tree_set_tracked_samples(&t, 2, tracked_samples);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SAMPLES);
    tracked_samples[1] = (tsk_id_t) tsk_treeseq_get_num_nodes(&ts);
    ret = tsk_tree_set_tracked_samples(&t, 2, tracked_samples);
    CU_ASSERT_EQUAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tracked_samples[1] = 0;
    ret = tsk_tree_set_tracked_samples(&t, 2, tracked_samples);
    CU_ASSERT_EQUAL(ret, TSK_ERR_DUPLICATE_SAMPLE);

    tsk_treeseq_from_text(&other_ts, 10, paper_ex_nodes, paper_ex_edges, NULL, NULL, NULL,
            paper_ex_individuals, NULL);
    ret = tsk_tree_alloc(&other_t, &other_ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_tree_copy(&t, &t);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_tree_copy(&t, &other_t);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* TODO run checks for the various unsupported operations with flags */

    tsk_tree_free(&t);
    tsk_tree_free(&other_t);
    tsk_treeseq_free(&other_ts);
    tsk_treeseq_free(&ts);
}

static void
test_deduplicate_sites(void)
{
    int ret;
    // Modified from paper_ex
    const char *tidy_sites =
        "1      0\n"
        "4.5    0\n"
        "8.5    0\n";
    const char *tidy_mutations =
        "0      2   1\n"
        "0      1   2\n"
        "0      6   3\n"
        "0      3   4\n"
        "1      0   1\n"
        "1      2   2\n"
        "1      4   3\n"
        "1      5   4\n"
        "2      5   1\n"
        "2      7   2\n"
        "2      1   3\n"
        "2      0   4\n";
    const char *messy_sites =
        "1      0\n"
        "1      0\n"
        "1      0\n"
        "1      0\n"
        "4.5    0\n"
        "4.5    0\n"
        "4.5    0\n"
        "4.5    0\n"
        "8.5    0\n"
        "8.5    0\n"
        "8.5    0\n"
        "8.5    0\n";
    const char *messy_mutations =
        "0      2   1\n"
        "1      1   2\n"
        "2      6   3\n"
        "3      3   4\n"
        "4      0   1\n"
        "5      2   2\n"
        "6      4   3\n"
        "7      5   4\n"
        "8      5   1\n"
        "9      7   2\n"
        "10     1   3\n"
        "11     0   4\n";
    tsk_tbl_collection_t tidy, messy;

    ret = tsk_tbl_collection_alloc(&tidy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_alloc(&messy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    messy.sequence_length = 10;
    tidy.sequence_length = 10;
    parse_individuals(paper_ex_individuals, tidy.individuals);
    parse_nodes(paper_ex_nodes, tidy.nodes);
    parse_sites(tidy_sites, tidy.sites);
    parse_mutations(tidy_mutations, tidy.mutations);
    // test cleaning doesn't mess up the tidy one
    parse_individuals(paper_ex_individuals, messy.individuals);
    parse_nodes(paper_ex_nodes, messy.nodes);
    parse_sites(tidy_sites, messy.sites);
    parse_mutations(tidy_mutations, messy.mutations);

    ret = tsk_tbl_collection_deduplicate_sites(&messy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_site_tbl_equals(tidy.sites, messy.sites));
    CU_ASSERT_TRUE(tsk_mutation_tbl_equals(tidy.mutations, messy.mutations));

    tsk_site_tbl_clear(messy.sites);
    tsk_mutation_tbl_clear(messy.mutations);

    // test with the actual messy one
    parse_sites(messy_sites, messy.sites);
    parse_mutations(messy_mutations, messy.mutations);

    ret = tsk_tbl_collection_deduplicate_sites(&messy, 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(tsk_site_tbl_equals(tidy.sites, messy.sites));
    CU_ASSERT_TRUE(tsk_mutation_tbl_equals(tidy.mutations, messy.mutations));

    tsk_tbl_collection_free(&tidy);
    tsk_tbl_collection_free(&messy);
}

static void
test_deduplicate_sites_errors(void)
{
    int ret;
    tsk_tbl_collection_t tables;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 10;
    ret = tsk_site_tbl_add_row(tables.sites, 2, "A", 1, "m", 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_tbl_add_row(tables.sites, 2, "TT", 2, "MM", 2);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    ret = tsk_mutation_tbl_add_row(tables.mutations, 0, 0, -1,
            "T", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_tbl_add_row(tables.nodes, 0, 0, TSK_NULL,
            TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Negative position */
    tables.sites->position[0] = -1;
    ret = tsk_tbl_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_SITE_POSITION);
    tables.sites->position[0] = 2;

    /* unsorted position */
    tables.sites->position[1] = 0.5;
    ret = tsk_tbl_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_UNSORTED_SITES);
    tables.sites->position[1] = 2;

    /* negative site ID */
    tables.mutations->site[0] = -1;
    ret = tsk_tbl_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tables.mutations->site[0] = 0;

     /* site ID out of bounds */
    tables.mutations->site[0] = 2;
    ret = tsk_tbl_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    tables.mutations->site[0] = 0;

    /* Bad offset in metadata */
    tables.sites->metadata_offset[0] = 2;
    ret = tsk_tbl_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    tables.sites->metadata_offset[0] = 0;

    /* Bad length in metadata */
    tables.sites->metadata_offset[2] = 100;
    ret = tsk_tbl_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    tables.sites->metadata_offset[2] = 3;

    /* Bad offset in ancestral_state */
    tables.sites->ancestral_state_offset[0] = 2;
    ret = tsk_tbl_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    tables.sites->ancestral_state_offset[0] = 0;

    /* Bad length in ancestral_state */
    tables.sites->ancestral_state_offset[2] = 100;
    ret = tsk_tbl_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    tables.sites->ancestral_state_offset[2] = 3;

    ret = tsk_tbl_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL(ret, 0);

    tsk_tbl_collection_free(&tables);
}

static void
test_deduplicate_sites_multichar(void)
{
    int ret;
    tsk_tbl_collection_t tables;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 10;
    ret = tsk_site_tbl_add_row(tables.sites, 0, "AA", 1, "M", 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_tbl_add_row(tables.sites, 0, "0", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    ret = tsk_site_tbl_add_row(tables.sites, 1, "BBBBB", 5, "NNNNN", 5);
    CU_ASSERT_EQUAL_FATAL(ret, 2);
    ret = tsk_site_tbl_add_row(tables.sites, 1, "0", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 3);

    ret = tsk_tbl_collection_deduplicate_sites(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites->num_rows, 2);
    CU_ASSERT_EQUAL_FATAL(tables.sites->position[0], 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites->position[1], 1);
    CU_ASSERT_EQUAL_FATAL(tables.sites->ancestral_state[0], 'A');
    CU_ASSERT_EQUAL_FATAL(tables.sites->ancestral_state_offset[1], 1);
    CU_ASSERT_EQUAL_FATAL(tables.sites->metadata[0], 'M');
    CU_ASSERT_EQUAL_FATAL(tables.sites->metadata_offset[1], 1);

    CU_ASSERT_NSTRING_EQUAL(tables.sites->ancestral_state + 1, "BBBBB", 5);
    CU_ASSERT_EQUAL_FATAL(tables.sites->ancestral_state_offset[2], 6);
    CU_ASSERT_NSTRING_EQUAL(tables.sites->metadata + 1, "NNNNN", 5);
    CU_ASSERT_EQUAL_FATAL(tables.sites->metadata_offset[2], 6);

    tsk_tbl_collection_free(&tables);
}

static void
test_empty_tree_sequence(void)
{
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    tsk_tree_t t;
    tsk_id_t v;
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_alloc(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SEQUENCE_LENGTH);
    tsk_treeseq_free(&ts);
    tables.sequence_length = 1.0;
    ret = tsk_treeseq_alloc(&ts, &tables, TSK_BUILD_INDEXES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    verify_empty_tree_sequence(&ts, 1.0);

    ret = tsk_tree_alloc(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL_FATAL(t.left_root, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(t.left, 0);
    CU_ASSERT_EQUAL_FATAL(t.right, 1);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_parent(&t, 0, &v), TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_tree_free(&t);

    ret = tsk_tree_alloc(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_last(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL_FATAL(t.left_root, TSK_NULL);
    CU_ASSERT_EQUAL_FATAL(t.left, 0);
    CU_ASSERT_EQUAL_FATAL(t.right, 1);
    CU_ASSERT_EQUAL_FATAL(tsk_tree_get_parent(&t, 0, &v), TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_tree_free(&t);

    tsk_treeseq_free(&ts);
    tsk_tbl_collection_free(&tables);
}

static void
test_zero_edges(void)
{
    const char *nodes =
        "1  0   0\n"
        "1  0   0\n";
    const char *edges = "";
    const char *sites =
        "0.1  0\n"
        "0.2  0\n";
    const char *mutations =
        "0    0     1\n"
        "1    1     1\n";
    tsk_treeseq_t ts, tss;
    tsk_tree_t t;
    const char *haplotypes[] = {"10", "01"};
    char *haplotype;
    tsk_hapgen_t hapgen;
    unsigned int j;
    tsk_id_t samples, node_map;
    const tsk_id_t z = TSK_NULL;
    tsk_id_t parents[] = {
        z, z,
    };
    int ret;

    tsk_treeseq_from_text(&ts, 2, nodes, edges, NULL, sites, mutations, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&ts), 2.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&ts), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&ts), 1);
    tsk_treeseq_print_state(&ts, _devnull);

    verify_trees(&ts, 1, parents);

    ret = tsk_tree_alloc(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(t.left, 0);
    CU_ASSERT_EQUAL(t.right, 2);
    CU_ASSERT_EQUAL(t.parent[0], TSK_NULL);
    CU_ASSERT_EQUAL(t.parent[1], TSK_NULL);
    CU_ASSERT_EQUAL(t.left_root, 0);
    CU_ASSERT_EQUAL(t.left_sib[0], TSK_NULL);
    CU_ASSERT_EQUAL(t.right_sib[0], 1);
    tsk_tree_print_state(&t, _devnull);
    tsk_tree_free(&t);

    ret = tsk_tree_alloc(&t, &ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tree_last(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(t.left, 0);
    CU_ASSERT_EQUAL(t.right, 2);
    CU_ASSERT_EQUAL(t.parent[0], TSK_NULL);
    CU_ASSERT_EQUAL(t.parent[1], TSK_NULL);
    CU_ASSERT_EQUAL(t.left_root, 0);
    CU_ASSERT_EQUAL(t.left_sib[0], TSK_NULL);
    CU_ASSERT_EQUAL(t.right_sib[0], 1);
    tsk_tree_print_state(&t, _devnull);
    tsk_tree_free(&t);

    /* We give pointers ot samples and node_map here as they must be non null */
    ret = tsk_treeseq_simplify(&ts, &samples, 0, 0, &tss, &node_map);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_samples(&tss), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_sequence_length(&tss), 2.0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_nodes(&tss), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_sites(&tss), 2);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_mutations(&tss), 0);
    CU_ASSERT_EQUAL(tsk_treeseq_get_num_trees(&tss), 1);
    tsk_treeseq_print_state(&ts, _devnull);

    ret = tsk_hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < 2; j++) {
        ret = tsk_hapgen_get_haplotype(&hapgen, (tsk_id_t) j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, haplotypes[j]);
    }
    tsk_hapgen_free(&hapgen);
    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&tss);
}




int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        /* simplest example tests */
        {"test_simplest_records", test_simplest_records},
        {"test_simplest_nonbinary_records", test_simplest_nonbinary_records},
        {"test_simplest_unary_records", test_simplest_unary_records},
        {"test_simplest_non_sample_leaf_records", test_simplest_non_sample_leaf_records},
        {"test_simplest_degenerate_multiple_root_records",
            test_simplest_degenerate_multiple_root_records},
        {"test_simplest_multiple_root_records", test_simplest_multiple_root_records},
        {"test_simplest_zero_root_tree", test_simplest_zero_root_tree},
        {"test_simplest_root_mutations", test_simplest_root_mutations},
        {"test_simplest_back_mutations", test_simplest_back_mutations},
        {"test_simplest_general_samples", test_simplest_general_samples},
        {"test_simplest_holey_tree_sequence", test_simplest_holey_tree_sequence},
        {"test_simplest_holey_tsk_treeseq_zero_roots",
            test_simplest_holey_tsk_treeseq_zero_roots},
        {"test_simplest_holey_tsk_treeseq_mutation_parents",
            test_simplest_holey_tsk_treeseq_mutation_parents},
        {"test_simplest_initial_gap_tree_sequence", test_simplest_initial_gap_tree_sequence},
        {"test_simplest_initial_gap_zero_roots", test_simplest_initial_gap_zero_roots},
        {"test_simplest_initial_gap_tsk_treeseq_mutation_parents",
            test_simplest_initial_gap_tsk_treeseq_mutation_parents},
        {"test_simplest_final_gap_tree_sequence", test_simplest_final_gap_tree_sequence},
        {"test_simplest_final_gap_tsk_treeseq_mutation_parents",
            test_simplest_final_gap_tsk_treeseq_mutation_parents},
        {"test_simplest_individuals", test_simplest_individuals},
        {"test_simplest_bad_individuals", test_simplest_bad_individuals},
        {"test_simplest_bad_edges", test_simplest_bad_edges},
        {"test_simplest_bad_indexes", test_simplest_bad_indexes},
        {"test_simplest_bad_migrations", test_simplest_bad_migrations},
        {"test_simplest_migration_simplify", test_simplest_migration_simplify},
        {"test_simplest_overlapping_parents", test_simplest_overlapping_parents},
        {"test_simplest_contradictory_children", test_simplest_contradictory_children},
        {"test_simplest_overlapping_edges_simplify",
            test_simplest_overlapping_edges_simplify},
        {"test_simplest_overlapping_unary_edges_simplify",
            test_simplest_overlapping_unary_edges_simplify},
        {"test_simplest_overlapping_unary_edges_internal_samples_simplify",
            test_simplest_overlapping_unary_edges_internal_samples_simplify},
        {"test_simplest_reduce_site_topology", test_simplest_reduce_site_topology},
        {"test_simplest_population_filter", test_simplest_population_filter},
        {"test_simplest_individual_filter", test_simplest_individual_filter},

        /* Single tree tests */
        {"test_single_tree_good_records", test_single_tree_good_records},
        {"test_single_nonbinary_tree_good_records",
            test_single_nonbinary_tree_good_records},
        {"test_single_tree_bad_records", test_single_tree_bad_records},
        {"test_single_tree_good_mutations", test_single_tree_good_mutations},
        {"test_single_tree_bad_mutations", test_single_tree_bad_mutations},
        {"test_single_tree_iter", test_single_tree_iter},
        {"test_single_tree_general_samples_iter", test_single_tree_general_samples_iter},
        {"test_single_nonbinary_tree_iter", test_single_nonbinary_tree_iter},
        {"test_single_tree_iter_times", test_single_tree_iter_times},
        {"test_single_tree_simplify", test_single_tree_simplify},
        {"test_single_tree_compute_mutation_parents", test_single_tree_compute_mutation_parents},

        /* Multi tree tests */
        {"test_simple_multi_tree", test_simple_multi_tree},
        {"test_nonbinary_multi_tree", test_nonbinary_multi_tree},
        {"test_unary_multi_tree", test_unary_multi_tree},
        {"test_internal_sample_multi_tree", test_internal_sample_multi_tree},
        {"test_internal_sample_simplified_multi_tree",
            test_internal_sample_simplified_multi_tree},
        {"test_left_to_right_multi_tree", test_left_to_right_multi_tree},
        {"test_gappy_multi_tree", test_gappy_multi_tree},
        {"test_tsk_treeseq_bad_records", test_tsk_treeseq_bad_records},

        /* Diff iter tests */
        {"test_simple_diff_iter", test_simple_diff_iter},
        {"test_nonbinary_diff_iter", test_nonbinary_diff_iter},
        {"test_unary_diff_iter", test_unary_diff_iter},
        {"test_internal_sample_diff_iter", test_internal_sample_diff_iter},

        /* Sample sets */
        {"test_simple_sample_sets", test_simple_sample_sets},
        {"test_nonbinary_sample_sets", test_nonbinary_sample_sets},
        {"test_internal_sample_sample_sets", test_internal_sample_sample_sets},

        /* Misc */
        {"test_tree_errors", test_tree_errors},
        {"test_genealogical_nearest_neighbours_errors",
            test_genealogical_nearest_neighbours_errors},
        {"test_deduplicate_sites", test_deduplicate_sites},
        {"test_deduplicate_sites_errors", test_deduplicate_sites_errors},
        {"test_deduplicate_sites_multichar", test_deduplicate_sites_multichar},
        {"test_empty_tree_sequence", test_empty_tree_sequence},
        {"test_zero_edges", test_zero_edges},

        {NULL},
    };

    return test_main(tests, argc, argv);
}
