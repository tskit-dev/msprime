/*
** Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
**
** This file is part of msprime.
**
** msprime is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** msprime is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with msprime.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdarg.h>
#include <float.h>

#include <libconfig.h>
#include <gsl/gsl_math.h>

#include "msprime.h"
#include "err.h"

/* This file defines a crude CLI for msprime. It is intended for development
 * use only.
 */

typedef struct {
    double mutation_rate;
} mutation_params_t;

static void
fatal_error(const char *msg, ...)
{
    va_list argp;
    fprintf(stderr, "main:");
    va_start(argp, msg);
    vfprintf(stderr, msg, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}

static int
read_population_configuration(msp_t *msp, config_t *config)
{
    int ret = 0;
    int j;
    double growth_rate, initial_size;
    int sample_size;
    int num_populations;
    config_setting_t *s, *t;
    config_setting_t *setting = config_lookup(config, "population_configuration");

    if (setting == NULL) {
        fatal_error("population_configuration is a required parameter");
    }
    if (config_setting_is_list(setting) == CONFIG_FALSE) {
        fatal_error("population_configuration must be a list");
    }
    num_populations = config_setting_length(setting);
    ret = msp_set_num_populations(msp, (size_t) num_populations);
    if (ret != 0) {
        fatal_error("Error reading number of populations");
    }
    for (j = 0; j < num_populations; j++) {
        s = config_setting_get_elem(setting, (unsigned int) j);
        if (s == NULL) {
            fatal_error("error reading population_configurations[%d]", j);
        }
        if (config_setting_is_group(s) == CONFIG_FALSE) {
            fatal_error("population_configurations[%d] not a group", j);
        }
        t = config_setting_get_member(s, "sample_size");
        if (t == NULL) {
            fatal_error("sample_size not specified");
        }
        sample_size = config_setting_get_int(t);
        t = config_setting_get_member(s, "growth_rate");
        if (t == NULL) {
            fatal_error("growth_rate not specified");
        }
        growth_rate = config_setting_get_float(t);
        t = config_setting_get_member(s, "initial_size");
        if (t == NULL) {
            fatal_error("initial_size not specified");
        }
        initial_size = config_setting_get_float(t);
        ret = msp_set_population_configuration(msp, j, (size_t) sample_size,
                initial_size, growth_rate);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

static int
read_demographic_events(msp_t *msp, config_t *config)
{
    int ret = 0;
    int j;
    const char *type;
    double time, growth_rate, initial_size, migration_rate, proportion;
    int num_demographic_events, population_id, matrix_index, source, dest;
    config_setting_t *s, *t;
    config_setting_t *setting = config_lookup(config, "demographic_events");

    if (setting == NULL) {
        fatal_error("demographic_events is a required parameter");
    }
    if (config_setting_is_list(setting) == CONFIG_FALSE) {
        fatal_error("demographic_events must be a list");
    }
    num_demographic_events = config_setting_length(setting);
    for (j = 0; j < num_demographic_events; j++) {
        s = config_setting_get_elem(setting, (unsigned int) j);
        if (s == NULL) {
            fatal_error("error reading demographic_events[%d]", j);
        }
        if (config_setting_is_group(s) == CONFIG_FALSE) {
            fatal_error("demographic_events[%d] not a group", j);
        }
        t = config_setting_get_member(s, "time");
        if (t == NULL) {
            fatal_error("time not specified");
        }
        time = config_setting_get_float(t);
        if (time < 0.0) {
            fatal_error("demographic event time must be > 0");
        }
        t = config_setting_get_member(s, "type");
        if (t == NULL) {
            fatal_error("type not specified");
        }
        type = config_setting_get_string(t);
        if (strcmp(type, "population_parameters_change") == 0) {
            growth_rate = GSL_NAN;
            t = config_setting_get_member(s, "growth_rate");
            if (t != NULL) {
                growth_rate = config_setting_get_float(t);
            }
            initial_size = GSL_NAN;
            t = config_setting_get_member(s, "initial_size");
            if (t != NULL) {
                initial_size = config_setting_get_float(t);
            }
            t = config_setting_get_member(s, "population_id");
            if (t == NULL) {
                fatal_error("population_id not specified");
            }
            population_id = config_setting_get_int(t);
            ret = msp_add_population_parameters_change(msp, time,
                    population_id, initial_size, growth_rate);
        } else if (strcmp(type, "migration_rate_change") == 0) {
            t = config_setting_get_member(s, "migration_rate");
            if (t == NULL) {
                fatal_error("migration_rate not specified");
            }
            migration_rate = config_setting_get_float(t);
            t = config_setting_get_member(s, "matrix_index");
            if (t == NULL) {
                fatal_error("matrix_index not specified");
            }
            matrix_index = config_setting_get_int(t);
            ret = msp_add_migration_rate_change(msp, time,
                    matrix_index, migration_rate);
        } else if (strcmp(type, "mass_migration") == 0) {
            t = config_setting_get_member(s, "proportion");
            if (t == NULL) {
                fatal_error("proportion not specified");
            }
            proportion = config_setting_get_float(t);
            t = config_setting_get_member(s, "source");
            if (t == NULL) {
                fatal_error("matrix_index not specified");
            }
            source = config_setting_get_int(t);
            t = config_setting_get_member(s, "dest");
            if (t == NULL) {
                fatal_error("matrix_index not specified");
            }
            dest = config_setting_get_int(t);
            ret = msp_add_mass_migration(msp, time, source, dest, proportion);
        } else {
            fatal_error("unknown demographic event type '%s'", type);
        }
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

static int
read_migration_matrix(msp_t *msp, config_t *config)
{
    int ret = 0;
    size_t j, size;
    double *migration_matrix = NULL;
    config_setting_t *s;
    config_setting_t *setting = config_lookup(config, "migration_matrix");

    if (setting == NULL) {
        fatal_error("migration_matrix is a required parameter");
    }
    if (config_setting_is_array(setting) == CONFIG_FALSE) {
        fatal_error("migration_matrix must be an array");
    }
    size = (size_t) config_setting_length(setting);
    migration_matrix = malloc(size * sizeof(double));
    if (migration_matrix == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < size; j++) {
        s = config_setting_get_elem(setting, (unsigned int) j);
        if (s == NULL) {
            fatal_error("error reading migration_matrix[%d]", j);
        }
        migration_matrix[j] = (double) config_setting_get_float(s);
    }
    ret = msp_set_migration_matrix(msp, size, migration_matrix);
out:
    if (migration_matrix != NULL) {
        free(migration_matrix);
    }
    return ret;
}


static int
read_recomb_map(uint32_t num_loci, recomb_map_t *recomb_map, config_t *config)
{
    int ret = 0;
    size_t j, size;
    double *rates = NULL;
    double *coordinates = NULL;
    config_setting_t *row, *s;
    config_setting_t *setting = config_lookup(config, "recombination_map");

    if (setting == NULL) {
        fatal_error("recombination_map is a required parameter");
    }
    if (config_setting_is_list(setting) == CONFIG_FALSE) {
        fatal_error("recombination_map must be a list");
    }
    size = (size_t) config_setting_length(setting);
    rates = malloc(size * sizeof(double));
    coordinates = malloc(size * sizeof(double));
    if (rates == NULL || coordinates == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < size; j++) {
        row = config_setting_get_elem(setting, (unsigned int) j);
        if (row == NULL) {
            fatal_error("error reading recombination_map[%d]", j);
        }
        if (config_setting_is_array(row) == CONFIG_FALSE) {
            fatal_error("recombination_map[%d] not an array", j);
        }
        if (config_setting_length(row) != 2) {
            fatal_error(
                "recombination_map entries must be [coord, rate] pairs");
        }
        s = config_setting_get_elem(row, 0);
        if (!config_setting_is_number(s)) {
            fatal_error("recombination_map entries must be numbers");
        }
        coordinates[j] = config_setting_get_float(s);
        s = config_setting_get_elem(row, 1);
        if (!config_setting_is_number(s)) {
            fatal_error("recombination_map entries must be numbers");
        }
        rates[j] = config_setting_get_float(s);
    }
    ret = recomb_map_alloc(recomb_map, num_loci, coordinates[size - 1],
            coordinates, rates, size);
out:
    if (rates != NULL) {
        free(rates);
    }
    if (coordinates != NULL) {
        free(coordinates);
    }
    return ret;
}

static int
get_configuration(gsl_rng *rng, msp_t *msp, mutation_params_t *mutation_params,
        recomb_map_t *recomb_map, char **output_file, const char *filename)
{
    int ret = 0;
    int err;
    int int_tmp;
    char *str;
    const char *str_tmp;
    double rho;
    config_t *config = malloc(sizeof(config_t));

    if (config == NULL) {
        fatal_error("no memory");
    }
    config_init(config);
    err = config_read_file(config, filename);
    if (err == CONFIG_FALSE) {
        fatal_error("configuration error:%s at line %d in file %s\n",
                config_error_text(config), config_error_line(config),
                filename);
    }
    if (config_lookup_int(config, "random_seed", &int_tmp) == CONFIG_FALSE) {
        fatal_error("random_seed is a required parameter");
    }
    gsl_rng_set(rng,  (unsigned long) int_tmp);
    if (config_lookup_int(config, "sample_size", &int_tmp) == CONFIG_FALSE) {
        fatal_error("sample_size is a required parameter");
    }
    ret = msp_alloc(msp, (size_t) int_tmp, rng);
    if (config_lookup_int(config, "num_loci", &int_tmp) == CONFIG_FALSE) {
        fatal_error("num_loci is a required parameter");
    }
    ret = msp_set_num_loci(msp, (size_t) int_tmp);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    /* Set the mutation rate */
    if (config_lookup_float(config,
            "mutation_rate", &mutation_params->mutation_rate)
            == CONFIG_FALSE) {
        fatal_error("mutation_rate is a required parameter");
    }
    if (config_lookup_int(config, "avl_node_block_size", &int_tmp)
            == CONFIG_FALSE) {
        fatal_error("avl_node_block_size is a required parameter");
    }
    msp_set_avl_node_block_size(msp, (size_t) int_tmp);
    if (config_lookup_int(config, "segment_block_size", &int_tmp)
            == CONFIG_FALSE) {
        fatal_error("segment_block_size is a required parameter");
    }
    msp_set_segment_block_size(msp, (size_t) int_tmp);
    if (config_lookup_int(config, "node_mapping_block_size", &int_tmp)
            == CONFIG_FALSE) {
        fatal_error("node_mapping_block_size is a required parameter");
    }
    msp_set_node_mapping_block_size(msp, (size_t) int_tmp);
    if (config_lookup_int(config, "coalescence_record_block_size", &int_tmp)
            == CONFIG_FALSE) {
        fatal_error("coalescence_record_block_size is a required parameter");
    }
    msp_set_coalescence_record_block_size(msp, (size_t) int_tmp);
    if (config_lookup_int(config, "max_memory", &int_tmp)
            == CONFIG_FALSE) {
        fatal_error("max_memory is a required parameter");
    }
    msp_set_max_memory(msp, (size_t) int_tmp * 1024 * 1024);

    ret = read_population_configuration(msp, config);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    ret = read_migration_matrix(msp, config);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    ret = read_demographic_events(msp, config);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    if (config_lookup_string(config, "output_file", &str_tmp)
            == CONFIG_FALSE) {
        fatal_error("output_file is a required parameter");
    }
    ret = read_recomb_map((uint32_t) msp_get_num_loci(msp),
            recomb_map, config);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    rho = recomb_map_get_per_locus_recombination_rate(recomb_map);
    ret = msp_set_scaled_recombination_rate(msp, rho);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    /* Create a copy of output_file to return */
    str = malloc(strlen(str_tmp) + 1);
    if (str == NULL) {
        fatal_error("Out of memory");
    }
    strcpy(str, str_tmp);
    *output_file = str;
    config_destroy(config);
    free(config);
    return ret;
}

static void
print_variants(tree_sequence_t *ts)
{
    int ret = 0;
    vargen_t *vg = calloc(1, sizeof(vargen_t));
    uint32_t j;
    double x;
    char *variant;

    printf("variants (%d) \n", (int) ts->num_mutations);
    if (vg == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = vargen_alloc(vg, ts);
    if (ret != 0) {
        goto out;
    }
    j = 0;
    while ((ret = vargen_next(vg, &x, &variant)) == 1) {
        printf("%d\t%f\t%s\n", j, x, variant);
        j++;
    }
    if (ret != 0) {
        goto out;
    }
out:
    if (vg != NULL) {
        vargen_free(vg);
        free(vg);
    }
    if (ret != 0) {
        printf("error occured:%d:%s\n", ret, msp_strerror(ret));
    }
}

static void
print_haplotypes(tree_sequence_t *ts)
{
    int ret = 0;
    hapgen_t *hg = calloc(1, sizeof(hapgen_t));
    size_t num_mutations = tree_sequence_get_num_mutations(ts);
    mutation_t *mutations = malloc(num_mutations * sizeof(mutation_t));
    uint32_t j;
    char *haplotype;

    printf("haplotypes \n");
    if (hg == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = hapgen_alloc(hg, ts);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < ts->sample_size; j++) {
        ret = hapgen_get_haplotype(hg, j, &haplotype);
        if (ret < 0) {
            goto out;
        }
        printf("%d\t%s\n", j, haplotype);
    }
    /* Get the mutations, reset them, redo the same thing to check */
    ret = tree_sequence_get_mutations(ts, mutations);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_set_mutations(ts, num_mutations, mutations,
            "{}", "{}");
    if (ret != 0) {
        goto out;
    }
    hapgen_free(hg);
    printf("checking set_mutations\n");
    ret = hapgen_alloc(hg, ts);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < ts->sample_size; j++) {
        ret = hapgen_get_haplotype(hg, j, &haplotype);
        if (ret < 0) {
            goto out;
        }
        printf("%d\t%s\n", j, haplotype);
    }

out:
    if (hg != NULL) {
        hapgen_free(hg);
        free(hg);
    }
    if (mutations != NULL) {
        free(mutations);
    }
    if (ret != 0) {
        printf("error occured:%d:%s\n", ret, msp_strerror(ret));
    }
}

static void
print_stats(tree_sequence_t *ts)
{
    int ret = 0;
    uint32_t j;
    uint32_t sample_size = tree_sequence_get_sample_size(ts) / 2;
    uint32_t *sample = malloc(sample_size * sizeof(uint32_t));
    double pi;

    if (sample == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    for (j = 0; j < sample_size; j++) {
        sample[j] = j;
    }
    ret = tree_sequence_get_pairwise_diversity(ts, sample,
        sample_size, &pi);
    if (ret != 0) {
        goto out;
    }
    printf("pi = %f\n", pi);
out:
    if (sample != NULL) {
        free(sample);
    }
    if (ret != 0) {
        printf("error occured:%d:%s\n", ret, msp_strerror(ret));
    }
}



static void
print_newick_trees(tree_sequence_t *ts)
{
    int ret = 0;
    newick_converter_t *nc = calloc(1, sizeof(newick_converter_t));
    double length;
    char *tree;

    printf("converting newick trees\n");
    if (nc == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    /* We're using an Ne of 0.25 here throughout to cancel 4Ne conversions */
    ret = newick_converter_alloc(nc, ts, 4, 0.25);
    if (ret != 0) {
        goto out;
    }
    while ((ret = newick_converter_next(nc, &length, &tree)) == 1) {
        printf("Tree: %f: %s\n", length, tree);
        newick_converter_print_state(nc);
    }
    if (ret != 0) {
        goto out;
    }
out:
    if (nc != NULL) {
        newick_converter_free(nc);
        free(nc);
    }
    if (ret != 0) {
        printf("error occured:%d:%s\n", ret, msp_strerror(ret));
    }
}

static void
print_tree_sequence(tree_sequence_t *ts)
{
    int ret = 0;
    size_t j;
    size_t num_records = tree_sequence_get_num_coalescence_records(ts);
    uint32_t mrca;
    double length;
    sparse_tree_t tree;
    node_record_t *records_in, *records_out, *record;
    coalescence_record_t cr;
    tree_diff_iterator_t *iter = calloc(1, sizeof(tree_diff_iterator_t));
    sparse_tree_iterator_t *sparse_iter = calloc(1, sizeof(sparse_tree_iterator_t));
    uint32_t tracked_leaves[] = {1, 2};

    if (iter == NULL || sparse_iter == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    printf("Records:\n");
    for (j = 0; j < num_records; j++) {
        if (tree_sequence_get_record(ts, j, &cr, MSP_ORDER_TIME) != 0) {
            fatal_error("tree sequence out of bounds\n");
        }
        printf("\t%f\t%f\t%d\t%d\t%d\t%f\n", cr.left, cr.right, cr.children[0],
                cr.children[1], cr.node, cr.time);
    }
    ret = tree_diff_iterator_alloc(iter, ts);
    if (ret != 0) {
        goto out;
    }
    printf("Tree diffs:\n");
    tree_diff_iterator_print_state(iter);
    while ((ret = tree_diff_iterator_next(
                    iter, &length, &records_out, &records_in)) == 1) {
        tree_diff_iterator_print_state(iter);
        printf("New tree: %f\n", length);
        printf("Nodes In:\n");
        record = records_in;
        while (record != NULL) {
            printf("\t(%d\t%d)\t%d\n", record->children[0],
                    record->children[1], record->node);
            record = record->next;
        }
        printf("Nodes Out:\n");
        record = records_out;
        while (record != NULL) {
            printf("\t(%d\t%d)\t%d\n", record->children[0],
                    record->children[1], record->node);
            record = record->next;
        }
    }
    if (ret != 0) {
        goto out;
    }
    ret = tree_diff_iterator_free(iter);
    if (ret != 0) {
        goto out;
    }

    /* sparse trees */
    ret = tree_sequence_alloc_sparse_tree(ts, &tree, tracked_leaves,
            sizeof(tracked_leaves) / sizeof(uint32_t), MSP_COUNT_LEAVES);
    if (ret != 0) {
        goto out;
    }
    ret = sparse_tree_iterator_alloc(sparse_iter, ts, &tree);
    if (ret != 0) {
        goto out;
    }
    printf("Sparse trees:\n");
    while ((ret = sparse_tree_iterator_next(sparse_iter)) == 1) {
        printf("New tree: %f (%d)\n", tree.right - tree.left,
                (int) tree.num_nodes);
        sparse_tree_iterator_print_state(sparse_iter);
        /* print some mrcas */
        printf("MRCAS:\n");
        for (j = 0; j < tree.num_nodes; j++) {
            ret = sparse_tree_get_mrca(&tree, 1, (uint32_t) j, &mrca);
            if (ret != 0) {
                goto out;
            }
            printf("\t%d %d -> %d\n", 1, (int) j, mrca);
        }
    }
    sparse_tree_iterator_free(sparse_iter);
    sparse_tree_free(&tree);
out:
    if (iter != NULL) {
        free(iter);
    }
    if (sparse_iter != NULL) {
        free(sparse_iter);
    }
    if (ret != 0) {
        fatal_error("ERROR: %d: %s\n", ret, msp_strerror(ret));
    }
}

static void
run_simulate(char *conf_file)
{
    int ret = -1;
    int result, j;
    double start_time, end_time;
    mutation_params_t mutation_params;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    msp_t *msp = calloc(1, sizeof(msp_t));
    char *output_file = NULL;
    tree_sequence_t *tree_seq = calloc(1, sizeof(tree_sequence_t));
    recomb_map_t *recomb_map = calloc(1, sizeof(recomb_map_t));
    mutgen_t *mutgen = calloc(1, sizeof(mutgen_t));

    if (rng == NULL || msp == NULL || tree_seq == NULL || recomb_map == NULL
            || mutgen == NULL) {
        goto out;
    }
    ret = get_configuration(rng, msp, &mutation_params, recomb_map,
            &output_file, conf_file);
    if (ret != 0) {
        goto out;
    }
    /* recomb_map_print_state(recomb_map); */
    ret = msp_initialise(msp);
    if (ret != 0) {
        goto out;
    }
    /* print out the demographic event debug state */
    start_time = 0;
    do {

        ret = msp_debug_demography(msp, &end_time);
        printf("interval %f - %f\n", start_time, end_time);
        msp_print_state(msp);
        start_time = end_time;
    } while (! gsl_isinf(end_time));
    if (ret != 0) {
        goto out;
    }
    ret = msp_reset(msp);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < 4; j++) {
        ret = msp_reset(msp);
        if (ret != 0) {
            goto out;
        }
        printf("Simulation run %d::\n", j);
        result = 1;
        while (result == 1) {
            result = msp_run(msp, DBL_MAX, 1);
            if (result < 0) {
                ret = result;
                goto out;
            }
            msp_verify(msp);
            /* ret = msp_print_state(msp); */
        }
        ret = msp_print_state(msp);
        if (ret != 0) {
            goto out;
        }
    }

    recomb_map_print_state(recomb_map);
    /* Create the tree_sequence from the state of the simulator.
     * We want to use coalescent time here, so use an Ne of 1/4
     * to cancel scaling factor. */
    ret = tree_sequence_create(tree_seq, msp, recomb_map, 0.25);
    if (ret != 0) {
        goto out;
    }

    print_tree_sequence(tree_seq);

    ret = mutgen_alloc(mutgen, tree_seq, mutation_params.mutation_rate, rng);
    if (ret != 0) {
        goto out;
    }
    ret = mutgen_generate(mutgen);
    if (ret != 0) {
        goto out;
    }
    mutgen_print_state(mutgen);
    ret = tree_sequence_set_mutations(tree_seq, mutgen->num_mutations,
            mutgen->mutations, mutgen->parameters, mutgen->environment);
    if (ret != 0) {
        goto out;
    }
    tree_sequence_print_state(tree_seq);
    print_stats(tree_seq);

    if (0) {
        for (j = 0; j < 1; j++) {
            ret = tree_sequence_dump(tree_seq, output_file, 0);
            if (ret != 0) {
                goto out;
            }
            tree_sequence_free(tree_seq);
            memset(tree_seq, 0, sizeof(tree_sequence_t));
            ret = tree_sequence_load(tree_seq, output_file, 0);
            if (ret != 0) {
                goto out;
            }
        }
        tree_sequence_print_state(tree_seq);

        print_newick_trees(tree_seq);

        print_haplotypes(tree_seq);
        print_variants(tree_seq);
        print_tree_sequence(tree_seq);

        tree_sequence_print_state(tree_seq);
        print_tree_sequence(tree_seq);
        print_haplotypes(tree_seq);
        tree_sequence_print_state(tree_seq);
        tree_sequence_print_state(tree_seq);
    }
out:
    if (msp != NULL) {
        msp_free(msp);
        free(msp);
    }
    if (tree_seq != NULL) {
        tree_sequence_free(tree_seq);
        free(tree_seq);
    }
    if (output_file != NULL) {
        free(output_file);
    }
    if (recomb_map != NULL) {
        recomb_map_free(recomb_map);
        free(recomb_map);
    }
    if (mutgen != NULL) {
        mutgen_free(mutgen);
        free(mutgen);
    }
    if (rng != NULL) {
        gsl_rng_free(rng);
    }
    if (ret != 0) {
        printf("error occured:%d:%s\n", ret, msp_strerror(ret));
    }
}

int
main(int argc, char** argv)
{
    if (argc != 2) {
        fatal_error("usage: %s CONFIG_FILE", argv[0]);
    }
    run_simulate(argv[1]);
    return EXIT_SUCCESS;
}
