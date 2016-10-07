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

static void
fatal_library_error(int err, const char *msg, ...)
{
    va_list argp;
    fprintf(stderr, "error:");
    va_start(argp, msg);
    vfprintf(stderr, msg, argp);
    va_end(argp);
    fprintf(stderr, ":%d:'%s'\n", err, msp_strerror(err));
    exit(EXIT_FAILURE);
}

static int
read_samples(config_t *config, size_t *sample_size, sample_t **samples)
{
    int ret = 0;
    size_t j, n;
    sample_t *ret_samples = NULL;
    config_setting_t *s, *t;
    config_setting_t *setting = config_lookup(config, "samples");

    if (setting == NULL) {
        fatal_error("samples is a required parameter");
    }
    if (config_setting_is_list(setting) == CONFIG_FALSE) {
        fatal_error("samples must be a list");
    }
    n = (size_t) config_setting_length(setting);
    ret_samples = malloc(n * sizeof(sample_t));
    if (ret_samples == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < n; j++) {
        s = config_setting_get_elem(setting, (unsigned int) j);
        if (s == NULL) {
            fatal_error("error reading samples[%d]", j);
        }
        if (config_setting_is_group(s) == CONFIG_FALSE) {
            fatal_error("samples[%d] not a group", j);
        }
        t = config_setting_get_member(s, "population");
        if (t == NULL) {
            fatal_error("population not specified");
        }
        ret_samples[j].population_id = (uint8_t) config_setting_get_int(t);
        t = config_setting_get_member(s, "time");
        if (t == NULL) {
            fatal_error("population not specified");
        }
        ret_samples[j].time = config_setting_get_float(t);
    }
    *samples = ret_samples;
    *sample_size = n;
out:
    return ret;
}



static int
read_population_configuration(msp_t *msp, config_t *config)
{
    int ret = 0;
    int j;
    double growth_rate, initial_size;
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
        ret = msp_set_population_configuration(msp, j, initial_size,
                growth_rate);
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
    double time, growth_rate, initial_size, migration_rate, proportion,
           intensity;
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
                fatal_error("source not specified");
            }
            source = config_setting_get_int(t);
            t = config_setting_get_member(s, "dest");
            if (t == NULL) {
                fatal_error("dest not specified");
            }
            dest = config_setting_get_int(t);
            ret = msp_add_mass_migration(msp, time, source, dest, proportion);
        } else if (strcmp(type, "bottleneck") == 0) {
            t = config_setting_get_member(s, "intensity");
            if (t == NULL) {
                fatal_error("intensity not specified");
            }
            intensity = config_setting_get_float(t);
            t = config_setting_get_member(s, "population_id");
            if (t == NULL) {
                fatal_error("population_id not specified");
            }
            population_id = config_setting_get_int(t);
            ret = msp_add_bottleneck(msp, time, population_id, intensity);
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
        recomb_map_t *recomb_map, const char *filename)
{
    int ret = 0;
    int err;
    int int_tmp;
    double rho;
    size_t sample_size;
    sample_t *samples = NULL;
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

    ret = read_samples(config, &sample_size, &samples);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    ret = msp_alloc(msp, sample_size, samples, rng);
    if (config_lookup_int(config, "num_loci", &int_tmp) == CONFIG_FALSE) {
        fatal_error("num_loci is a required parameter");
    }
    ret = msp_set_num_loci(msp, (size_t) int_tmp);
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
    config_destroy(config);
    free(config);
    free(samples);
    return ret;
}

static void
print_variants(tree_sequence_t *ts)
{
    int ret = 0;
    vargen_t vg;
    uint32_t j, k;
    mutation_t *mut;
    char *genotypes = malloc(tree_sequence_get_sample_size(ts) * sizeof(char));

    if (genotypes == NULL) {
        fatal_error("no memory");
    }
    printf("variants (%d) \n", (int) ts->num_mutations);
    ret = vargen_alloc(&vg, ts, 0);
    if (ret != 0) {
        fatal_library_error(ret, "vargen_alloc");
    }
    j = 0;
    while ((ret = vargen_next(&vg, &mut, genotypes)) == 1) {
        printf("%d\t%f\t", j, mut->position);
        for (k = 0; k < tree_sequence_get_sample_size(ts); k++) {
            printf("%d", genotypes[k]);
        }
        printf("\n");
        j++;
    }
    if (ret != 0) {
        fatal_library_error(ret, "vargen_next");
    }
    if (j != ts->num_mutations) {
        printf("ERROR!! missing variants %d %d\n", j, (int) ts->num_mutations);
    }

    while ((ret = vargen_next(&vg, &mut, genotypes)) == 1) {
        /* this should never happen as the iterators should always
         * fail after they finish. */
        assert(0);
    }
    vargen_free(&vg);
    free(genotypes);
}

static void
print_haplotypes(tree_sequence_t *ts)
{
    int ret = 0;
    hapgen_t hg;
    uint32_t j;
    char *haplotype;

    printf("haplotypes \n");
    ret = hapgen_alloc(&hg, ts);
    if (ret != 0) {
        fatal_library_error(ret, "hapgen_alloc");
    }
    for (j = 0; j < ts->sample_size; j++) {
        ret = hapgen_get_haplotype(&hg, j, &haplotype);
        if (ret < 0) {
            fatal_library_error(ret, "hapgen_get_haplotype");
        }
        printf("%d\t%s\n", j, haplotype);
    }
    hapgen_free(&hg);
}

static void
print_ld_matrix(tree_sequence_t *ts)
{
    int ret;
    size_t num_mutations = tree_sequence_get_num_mutations(ts);
    mutation_t *mutations;
    double *r2 = malloc(num_mutations * sizeof(double));
    size_t j, k, num_r2_values;
    ld_calc_t ld_calc;

    if (r2 == NULL) {
        fatal_error("no memory");
    }
    ret = tree_sequence_get_mutations(ts, &mutations);
    if (ret != 0) {
        fatal_library_error(ret, "get_mutations");
    }
    ret = ld_calc_alloc(&ld_calc, ts);
    printf("alloc: ret = %d\n", ret);
    if (ret != 0) {
        fatal_library_error(ret, "ld_calc_alloc");
    }
    ld_calc_print_state(&ld_calc, stdout);
    for (j = 0; j < num_mutations; j++) {
        ret = ld_calc_get_r2_array(&ld_calc, j, MSP_DIR_FORWARD, num_mutations,
                DBL_MAX, r2, &num_r2_values);
        if (ret != 0) {
            fatal_library_error(ret, "ld_calc_get_r2_array");
        }
        for (k = 0; k < num_r2_values; k++) {
            printf("%d\t%f\t%d\t%f\t%.3f\n",
                (int) mutations[j].index, mutations[j].position,
                (int) mutations[j + k + 1].index,
                mutations[j + k + 1].position, r2[k]);
        }
    }
    free(r2);
    ret = ld_calc_free(&ld_calc);
    if (ret != 0) {
        fatal_library_error(ret, "ld_calc_write_table");
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
        fatal_error("no memory");
    }
    for (j = 0; j < sample_size; j++) {
        sample[j] = j;
    }
    ret = tree_sequence_get_pairwise_diversity(ts, sample,
        sample_size, &pi);
    if (ret != 0) {
        fatal_library_error(ret, "get_pairwise_diversity");
    }
    printf("pi = %f\n", pi);
    free(sample);
}

static void
print_vcf(tree_sequence_t *ts, unsigned int ploidy)
{
    int ret = 0;
    char *record = NULL;
    char *header = NULL;
    vcf_converter_t vc;

    ret = vcf_converter_alloc(&vc, ts, ploidy);
    if (ret != 0) {
        fatal_library_error(ret, "vcf alloc");
    }
    vcf_converter_print_state(&vc, stdout);
    printf("START VCF\n");
    ret = vcf_converter_get_header(&vc, &header);
    if (ret != 0) {
        fatal_library_error(ret, "vcf get header");
    }
    printf("%s", header);
    while ((ret = vcf_converter_next(&vc, &record)) == 1) {
        printf("%s", record);
    }
    if (ret != 0) {
        fatal_library_error(ret, "vcf next");
    }
    vcf_converter_free(&vc);
}

static void
print_newick_trees(tree_sequence_t *ts)
{
    int ret = 0;
    newick_converter_t nc;
    double length;
    char *tree;

    printf("converting newick trees\n");
    /* We're using an Ne of 0.25 here throughout to cancel 4Ne conversions */
    ret = newick_converter_alloc(&nc, ts, 4, 0.25);
    if (ret != 0) {
        fatal_library_error(ret, "newick alloc");
    }
    while ((ret = newick_converter_next(&nc, &length, &tree)) == 1) {
        printf("Tree: %f: %s\n", length, tree);
        newick_converter_print_state(&nc, stdout);
    }
    if (ret != 0) {
        fatal_library_error(ret, "newick next");
    }
    newick_converter_free(&nc);
}

static void
print_tree_sequence(tree_sequence_t *ts)
{
    int ret = 0;
    size_t j, k; //, k, l;
    size_t num_records = tree_sequence_get_num_coalescence_records(ts);
    sparse_tree_t tree;
    coalescence_record_t *cr;
    uint32_t tracked_leaves[] = {1, 2};

    // TODO tidy this up to make it more specific to the task of examining the
    // tree sequence itself.
    printf("Records:\n");
    for (j = 0; j < num_records; j++) {
        if (tree_sequence_get_record(ts, j, &cr, MSP_ORDER_TIME) != 0) {
            fatal_error("tree sequence out of bounds\n");
        }
        printf("\t%f\t%f\t%d\t%d\t%d\t%f\n", cr->left, cr->right, cr->children[0],
                cr->children[1], cr->node, cr->time);
    }

    tree_sequence_print_state(ts, stdout);
    /* sparse trees */
    ret = sparse_tree_alloc(&tree, ts, MSP_LEAF_COUNTS);
    if (ret != 0) {
        goto out;
    }
    ret = sparse_tree_set_tracked_leaves(&tree,
            sizeof(tracked_leaves) / sizeof(uint32_t), tracked_leaves);
    if (ret != 0) {
        goto out;
    }
    for (ret = sparse_tree_first(&tree); ret == 1;
            ret = sparse_tree_next(&tree)) {
        printf("New tree: %d: %f (%d)\n", (int) tree.index,
                tree.right - tree.left, (int) tree.num_nodes);
        /* sparse_tree_print_state(&tree, stdout); */
        k = tree.index;
        if (k > 3) {
            printf("Rewinding\n");
            for (j = 0; j < 2; j++) {
                ret = sparse_tree_prev(&tree);
                if (ret != 1) {
                    goto out;
                }
                printf("\twent back to %d\n", (int) tree.index);
            }
            printf("Going forward again\n");
            for (j = 0; j < 2; j++) {
                ret = sparse_tree_next(&tree);
                if (ret != 1) {
                    goto out;
                }
                printf("\twent forward to %d\n", (int) tree.index);
            }
        }
    }
    if (ret < 0) {
        goto out;
    }
    sparse_tree_free(&tree);
out:
    if (ret != 0) {
        fatal_error("ERROR: %d: %s\n", ret, msp_strerror(ret));
    }
}

static void
run_simulate(char *conf_file, char *output_file)
{
    int ret = -1;
    int result, j;
    double start_time, end_time;
    mutation_params_t mutation_params;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    msp_t *msp = calloc(1, sizeof(msp_t));
    tree_sequence_t *tree_seq = calloc(1, sizeof(tree_sequence_t));
    recomb_map_t *recomb_map = calloc(1, sizeof(recomb_map_t));
    mutgen_t *mutgen = calloc(1, sizeof(mutgen_t));

    /* TODO tidy this up a bit and make the behaviour configurable. */
    if (rng == NULL || msp == NULL || tree_seq == NULL || recomb_map == NULL
            || mutgen == NULL) {
        goto out;
    }
    ret = get_configuration(rng, msp, &mutation_params, recomb_map, conf_file);
    if (ret != 0) {
        goto out;
    }
    /* recomb_map_print_state(recomb_map); */
    ret = msp_initialise(msp);
    if (ret != 0) {
        goto out;
    }
    if (0) {
        /* print out the demographic event debug state */
        start_time = 0;
        do {

            ret = msp_debug_demography(msp, &end_time);
            printf("interval %f - %f\n", start_time, end_time);
            msp_print_state(msp, stdout);
            start_time = end_time;
        } while (! gsl_isinf(end_time));
        if (ret != 0) {
            goto out;
        }
        ret = msp_reset(msp);
        if (ret != 0) {
            goto out;
        }
    }
    for (j = 0; j < 1; j++) {
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
            /* ret = msp_print_state(msp, stdout); */
        }
        ret = msp_print_state(msp, stdout);
        if (ret != 0) {
            goto out;
        }
    }
    /* Create the tree_sequence from the state of the simulator.
     * We want to use coalescent time here, so use an Ne of 1/4
     * to cancel scaling factor. */
    ret = tree_sequence_create(tree_seq, msp, recomb_map, 0.25);
    if (ret != 0) {
        goto out;
    }

    ret = tree_sequence_add_provenance_string(tree_seq, "Tree Provenance!!!");
    if (ret != 0) {
        goto out;
    }
    ret = mutgen_alloc(mutgen, tree_seq, mutation_params.mutation_rate, rng);
    if (ret != 0) {
        goto out;
    }
    ret = mutgen_generate(mutgen);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_set_mutations(tree_seq, mutgen->num_mutations,
            mutgen->mutations);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_dump(tree_seq, output_file, 0);
    if (ret != 0) {
        goto out;
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

static void
load_tree_sequence(tree_sequence_t *ts, char *filename)
{
    int ret = tree_sequence_load(ts, filename, 0);
    if (ret != 0) {
        fatal_library_error(ret, "Load error");
    }
}

static void
run_ld(char *filename)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_ld_matrix(&ts);
    tree_sequence_free(&ts);
}

static void
run_haplotypes(char *filename)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_haplotypes(&ts);
    tree_sequence_free(&ts);
}

static void
run_variants(char *filename)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_variants(&ts);
    tree_sequence_free(&ts);
}

static void
run_vcf(char *filename)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_vcf(&ts, 1);
    tree_sequence_free(&ts);
}

static void
run_print(char *filename)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_tree_sequence(&ts);
    tree_sequence_free(&ts);
}

static void
run_newick(char *filename)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_newick_trees(&ts);
    tree_sequence_free(&ts);
}

static void
run_stats(char *filename)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_stats(&ts);
    tree_sequence_free(&ts);
}

int
main(int argc, char** argv)
{
    char *cmd;
    if (argc < 2) {
        fatal_error("usage: %s <cmd>", argv[0]);
    }
    cmd = argv[1];
    if (strncmp(cmd, "simulate", strlen(cmd)) == 0) {
        if (argc < 4) {
            fatal_error("usage: %s simulate CONFIG_FILE OUTPUT_FILE", argv[0]);
        }
        run_simulate(argv[2], argv[3]);
    } else if (strncmp(cmd, "ld", strlen(cmd)) == 0) {
        if (argc < 3) {
            fatal_error("usage: %s ld INPUT_FILE", argv[0]);
        }
        run_ld(argv[2]);
    } else if (strncmp(cmd, "haplotypes", strlen(cmd)) == 0) {
        if (argc < 3) {
            fatal_error("usage: %s haplotypes INPUT_FILE", argv[0]);
        }
        run_haplotypes(argv[2]);
    } else if (strncmp(cmd, "variants", strlen(cmd)) == 0) {
        if (argc < 3) {
            fatal_error("usage: %s variants INPUT_FILE", argv[0]);
        }
        run_variants(argv[2]);
    } else if (strncmp(cmd, "vcf", strlen(cmd)) == 0) {
        if (argc < 3) {
            fatal_error("usage: %s vcf INPUT_FILE", argv[0]);
        }
        run_vcf(argv[2]);
    } else if (strncmp(cmd, "print", strlen(cmd)) == 0) {
        if (argc < 3) {
            fatal_error("usage: %s print INPUT_FILE", argv[0]);
        }
        run_print(argv[2]);
    } else if (strncmp(cmd, "newick", strlen(cmd)) == 0) {
        if (argc < 3) {
            fatal_error("usage: %s newick INPUT_FILE", argv[0]);
        }
        run_newick(argv[2]);
    } else if (strncmp(cmd, "stats", strlen(cmd)) == 0) {
        if (argc < 3) {
            fatal_error("usage: %s stats INPUT_FILE", argv[0]);
        }
        run_stats(argv[2]);
    } else {
        fatal_error("Unknown command '%s'", cmd);
    }
    return EXIT_SUCCESS;
}
