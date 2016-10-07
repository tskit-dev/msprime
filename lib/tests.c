/*
** Copyright (C) 2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

/*
 * Unit tests for the low-level msprime API.
 */

#include "msprime.h"

#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <unistd.h>

#include <hdf5.h>
#include <gsl/gsl_math.h>
#include <CUnit/Basic.h>

/* Global variables used for test in state in the test suite */

char * _tmp_file_name;
FILE * _devnull;

typedef struct {
    double time;
    uint32_t population_id;
    double intensity;
} bottleneck_desc_t;

/* Simple utility to parse records so we can write declaritive
 * tests. This is not intended as a robust general input mechanism.
 */
static void
parse_text_records(size_t num_records, const char **text_records,
        coalescence_record_t **returned_records)
{
    size_t j, k, num_children;
    size_t MAX_LINE = 1024;
    char line[MAX_LINE];
    char sub_line[MAX_LINE];
    const char *whitespace = " \t";
    char *p, *q;
    coalescence_record_t *records = malloc(
        num_records * sizeof(coalescence_record_t));

    CU_ASSERT_FATAL(records != NULL);

    for (j = 0; j < num_records; j++) {
        CU_ASSERT_FATAL(text_records[j] != NULL);
        CU_ASSERT_FATAL(strlen(text_records[j]) < MAX_LINE - 1);
        strncpy(line, text_records[j], MAX_LINE);
        p = strtok(line, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        records[j].left = atof(p);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        records[j].right = atof(p);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        records[j].node = atoi(p);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        num_children = 1;
        q = p;
        while (*q != '\0') {
            if (*q == ',') {
                num_children++;
            }
            q++;
        }
        CU_ASSERT_FATAL(num_children >= 2);
        records[j].children = malloc(num_children * sizeof(uint32_t));
        records[j].num_children = num_children;
        CU_ASSERT_FATAL(records[j].children != NULL);
        strncpy(sub_line, p, MAX_LINE);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        records[j].time = atof(p);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p != NULL);
        records[j].population_id = atoi(p);
        p = strtok(NULL, whitespace);
        CU_ASSERT_FATAL(p == NULL);
        q = strtok(sub_line, ",");
        for (k = 0; k < num_children; k++) {
            CU_ASSERT_FATAL(q != NULL);
            records[j].children[k] = atoi(q);
            q = strtok(NULL, ",");
        }
        CU_ASSERT_FATAL(q == NULL);
    }
    *returned_records = records;
}

/* Frees records locally alloced by (e.g.) parse_text_records function.
 */
static void
free_local_records(size_t num_records, coalescence_record_t *records)
{
    size_t j;

    for (j = 0; j < num_records; j++) {
        free(records[j].children);
    }
    free(records);
}

static void
copy_record(coalescence_record_t *dest, coalescence_record_t *source)
{
    size_t j;

    CU_ASSERT_FATAL(dest != NULL);
    CU_ASSERT_FATAL(source != NULL);
    CU_ASSERT_FATAL(source->children != NULL);
    CU_ASSERT_FATAL(source->num_children >= 2);
    memcpy(dest, source, sizeof(coalescence_record_t));
    dest->children = malloc(source->num_children * sizeof(uint32_t));
    CU_ASSERT_FATAL(dest->children != NULL);
    for (j = 0; j < source->num_children; j++) {
        dest->children[j] = source->children[j];
    }
}

static void
verify_coalescence_records_equal(coalescence_record_t *r1,
        coalescence_record_t *r2, double scale)
{
    uint32_t j;
    double eps = 1e-6;

    /* if (r1->left * scale != r2->left) { */
    /*     printf("error: %f %f, %f %f\n", r1->left * scale, r2->left, r1->right, r2->right); */
    /*     printf("%d\n, ", r1->left * scale == r2->left); */
    /* } */
    CU_ASSERT_DOUBLE_EQUAL_FATAL(r1->left * scale, r2->left, eps);
    CU_ASSERT_DOUBLE_EQUAL_FATAL(r1->right * scale, r2->right, eps);
    CU_ASSERT_EQUAL_FATAL(r1->num_children, r2->num_children);
    for (j = 0; j < r1->num_children; j++) {
        CU_ASSERT_EQUAL(r1->children[j], r2->children[j]);
    }
    CU_ASSERT_EQUAL(r1->time, r2->time);
    CU_ASSERT_EQUAL(r1->population_id, r2->population_id);
}

/* Utility function to return a tree sequence for testing. It is the
 * callers responsilibility to free all memory.
 */
static tree_sequence_t *
get_example_tree_sequence(uint32_t sample_size,
        uint32_t num_historical_samples, uint32_t num_loci,
        double sequence_length, double scaled_recombination_rate,
        double mutation_rate, uint32_t num_bottlenecks,
        bottleneck_desc_t *bottlenecks)
{
    int ret;
    msp_t *msp = malloc(sizeof(msp_t));
    sample_t *samples = malloc(sample_size * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    tree_sequence_t *tree_seq = malloc(sizeof(tree_sequence_t));
    recomb_map_t *recomb_map = malloc(sizeof(recomb_map_t));
    mutgen_t *mutgen = malloc(sizeof(mutgen_t));
    coalescence_record_t *sim_records, *ts_record;
    uint32_t j;
    size_t num_records;
    sample_t sample;
    double positions[] = {0.0, 0.0};
    double rates[] = {0.0, 0.0};

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);
    CU_ASSERT_FATAL(recomb_map != NULL);
    gsl_rng_set(rng, 1);

    /* initialise the samples to zero for the default configuration */
    memset(samples, 0, sample_size * sizeof(sample_t));
    for (j = 0; j < num_historical_samples; j++) {
        samples[j].time = 0.1 * (j + 1);
    }
    ret = msp_alloc(msp, sample_size, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(msp, num_loci);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(msp, scaled_recombination_rate);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_bottlenecks; j++) {
        ret = msp_add_bottleneck(msp, bottlenecks[j].time,
                bottlenecks[j].population_id, bottlenecks[j].intensity);
        CU_ASSERT_EQUAL(ret, 0);
    }
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);

    rates[0] = scaled_recombination_rate;
    positions[1] = sequence_length;
    ret = recomb_map_alloc(recomb_map, num_loci, sequence_length,
            positions, rates, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Create the tree_sequence from the state of the simulator.
     * We want to use coalescent time here, so use an Ne of 1/4
     * to cancel scaling factor. */
    ret = tree_sequence_create(tree_seq, msp, recomb_map, 0.25);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(
            tree_sequence_get_num_coalescence_records(tree_seq),
            msp_get_num_coalescence_records(msp));
    CU_ASSERT_EQUAL_FATAL(
            tree_sequence_get_sample_size(tree_seq),
            msp_get_sample_size(msp));
    CU_ASSERT_FATAL(
            tree_sequence_get_num_nodes(tree_seq) >= sample_size);
    ret = msp_get_coalescence_records(msp, &sim_records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    num_records = msp_get_num_coalescence_records(msp);
    for (j = 0; j < num_records; j++) {
        ret = tree_sequence_get_record(tree_seq, j, &ts_record,
                MSP_ORDER_TIME);
        CU_ASSERT_EQUAL(ret, 0);
        verify_coalescence_records_equal(&sim_records[j], ts_record,
                sequence_length / num_loci);
    }
    for (j = num_records; j < num_records + 10; j++) {
        ret = tree_sequence_get_record(tree_seq, j, &ts_record,
                MSP_ORDER_TIME);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    }
    for (j = 0; j < sample_size; j++) {
        ret = tree_sequence_get_sample(tree_seq, j, &sample);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(sample.population_id, samples[j].population_id);
        CU_ASSERT_EQUAL(sample.time, samples[j].time);
    }
    for (j = sample_size; j < sample_size + 10; j++) {
        ret = tree_sequence_get_sample(tree_seq, j, &sample);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    }

    ret = mutgen_alloc(mutgen, tree_seq, mutation_rate, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = mutgen_generate(mutgen);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(tree_seq, mutgen->num_mutations,
            mutgen->mutations);
    CU_ASSERT_EQUAL(ret, 0);

    gsl_rng_free(rng);
    free(samples);
    msp_free(msp);
    free(msp);
    recomb_map_free(recomb_map);
    free(recomb_map);
    mutgen_free(mutgen);
    free(mutgen);
    return tree_seq;
}

tree_sequence_t **
get_example_nonbinary_tree_sequences(void)
{
    size_t max_examples = 1024;
    tree_sequence_t **ret = malloc(max_examples * sizeof(tree_sequence_t *));
    bottleneck_desc_t bottlenecks[] = {
        {0.1, 0, 0.5},
        {0.4, 0, 1.0},
    };
    bottleneck_desc_t other_bottlenecks[] = {
        {0.1, 0, 0.1},
        {0.1, 0, 0.9},
    };

    CU_ASSERT_FATAL(ret != NULL);
    ret[0] = get_example_tree_sequence(100, 0, 100, 100.0, 10.0, 1.0,
            1, bottlenecks);
    ret[1] = get_example_tree_sequence(10, 2, 100, 10.0, 1.0, 2.0,
            1, bottlenecks);
    ret[2] = get_example_tree_sequence(500, 10, 10, 1000.0, 0.5, 3.0,
            2, bottlenecks);
    ret[3] = get_example_tree_sequence(100, 0, 100, 1.0, 1.0, 0.0,
            2, other_bottlenecks);
    ret[4] = NULL;
    return ret;
}

tree_sequence_t **
get_example_tree_sequences(int include_nonbinary)
{
    size_t max_examples = 1024;
    int err, j, k;
    tree_sequence_t **ret = malloc(max_examples * sizeof(tree_sequence_t *));
    tree_sequence_t **nonbinary = NULL;
    CU_ASSERT_FATAL(ret != NULL);

    ret[0] = get_example_tree_sequence(10, 0, 100, 100.0, 1.0, 1.0, 0, NULL);
    ret[1] = get_example_tree_sequence(2, 0, 1, 0.1, 1.0, 1.0, 0, NULL);
    ret[2] = get_example_tree_sequence(3, 0, 3, 10.0, 10.0, 0.0, 0, NULL);
    ret[3] = get_example_tree_sequence(10, 0, UINT32_MAX, 10.0,
            9.31322575049e-08, 10.0, 0, NULL);
    k = 4;
    if (include_nonbinary) {
        nonbinary = get_example_nonbinary_tree_sequences();
        for (j = 0; nonbinary[j] != NULL; j++) {
            ret[k] = nonbinary[j];
            k++;
        }
        free(nonbinary);
    }
    ret[k] = NULL;
    for (j = 0; ret[j] != NULL; j++) {
        for (k = 0; k < j; k++) {
            err = tree_sequence_add_provenance_string(ret[j], "test provenance");
            CU_ASSERT_EQUAL_FATAL(err, 0);
        }
    }
    return ret;
}

/* Simple unit tests for the Fenwick tree API. */
static void
test_fenwick(void)
{
    fenwick_t t;
    int64_t s;
    size_t j, n;
    for (n = 1; n < 100; n++) {
        s = 0;
        CU_ASSERT(fenwick_alloc(&t, n) == 0);
        for (j = 1; j <= n; j++) {
            fenwick_increment(&t, j, (int64_t) j);
            s = s + (int64_t) j;
            CU_ASSERT(fenwick_get_value(&t, j) == (int64_t) j);
            CU_ASSERT(fenwick_get_cumulative_sum(&t, j) == s);
            CU_ASSERT(fenwick_get_total(&t) == s);
            CU_ASSERT(fenwick_find(&t, s) == j);
            fenwick_set_value(&t, j, 0);
            CU_ASSERT(fenwick_get_value(&t, j) == 0);
            CU_ASSERT(fenwick_get_cumulative_sum(&t, j) == s - (int64_t) j);
            fenwick_set_value(&t, j, (int64_t) j);
            CU_ASSERT(fenwick_get_value(&t, j) == (int64_t) j);
            /* Just make sure that we're seeing the same values even when
             * we expand.
             */
            CU_ASSERT(fenwick_expand(&t, 1) == 0);
        }
        CU_ASSERT(fenwick_free(&t) == 0);
    }
}

static void
test_vcf(void)
{
    int ret;
    char *str = NULL;
    unsigned int ploidy, num_variants;
    vcf_converter_t *vc = malloc(sizeof(vcf_converter_t));
    tree_sequence_t *ts = get_example_tree_sequence(10, 0, 100, 100.0, 1.0, 1.0,
            0, NULL);

    CU_ASSERT_FATAL(ts != NULL);
    CU_ASSERT_FATAL(vc != NULL);

    ret = vcf_converter_alloc(vc, ts, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = vcf_converter_alloc(vc, ts, 3);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = vcf_converter_alloc(vc, ts, 11);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    for (ploidy = 1; ploidy < 3; ploidy++) {
        ret = vcf_converter_alloc(vc, ts, ploidy);
        CU_ASSERT_FATAL(ret ==  0);
        vcf_converter_print_state(vc, _devnull);
        ret = vcf_converter_get_header(vc, &str);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_NSTRING_EQUAL("##", str, 2);
        num_variants = 0;
        while ((ret = vcf_converter_next(vc, &str)) == 1) {
            CU_ASSERT_NSTRING_EQUAL("1\t", str, 2);
            num_variants++;
        }
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_TRUE(num_variants > 0);
        vcf_converter_free(vc);
    }

    free(vc);
    tree_sequence_free(ts);
    free(ts);
}

static void
test_simple_recomb_map(void)
{
    int ret;
    recomb_map_t recomb_map;
    double positions[] = {0.0, 1.0};
    double rates[] = {0.0, 1.0};
    uint32_t num_loci[] = {1, 100, UINT32_MAX};
    size_t j;

    for (j = 0; j < sizeof(num_loci) / sizeof(uint32_t); j++) {
        ret = recomb_map_alloc(&recomb_map, num_loci[j], 1.0,
                positions, rates, 2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        recomb_map_print_state(&recomb_map, _devnull);
        CU_ASSERT_EQUAL(recomb_map_get_num_loci(&recomb_map), num_loci[j]);
        CU_ASSERT_EQUAL(recomb_map_get_size(&recomb_map), 2);
        CU_ASSERT_EQUAL(recomb_map_get_sequence_length(&recomb_map), 1.0);
        CU_ASSERT_EQUAL(
                recomb_map_genetic_to_phys(&recomb_map, num_loci[j]), 1.0);
        CU_ASSERT_EQUAL(
                recomb_map_phys_to_genetic(&recomb_map, 1.0), num_loci[j]);
        ret = recomb_map_free(&recomb_map);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
test_recomb_map_errors(void)
{
    int ret;
    recomb_map_t recomb_map;
    double positions[] = {0.0, 1.0, 2.0};
    double rates[] = {1.0, 2.0, 0.0};
    double values[] = {0.0, 1.0};

    ret = recomb_map_alloc(&recomb_map, 10, 1.0, positions, rates, 0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 10, 1.0, positions, rates, 1);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 0, 1.0, positions, rates, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    ret = recomb_map_alloc(&recomb_map, 10, 2.0, positions, rates, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);

    positions[0] = 1;
    ret = recomb_map_alloc(&recomb_map, 10, 1.0, positions, rates, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);
    positions[0] = 0;

    positions[1] = 3.0;
    ret = recomb_map_alloc(&recomb_map, 10, 2.0, positions, rates, 3);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_RECOMBINATION_MAP);
    recomb_map_free(&recomb_map);
    positions[1] = 1.0;

    ret = recomb_map_alloc(&recomb_map, 10, 2.0, positions, rates, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = recomb_map_genetic_to_phys_bulk(&recomb_map, values, 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* values must be increasing */
    values[0] = 2.0;
    ret = recomb_map_genetic_to_phys_bulk(&recomb_map, values, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_GENERIC);

    /* values must be <= num_loci */
    values[0] = 1000;
    ret = recomb_map_genetic_to_phys_bulk(&recomb_map, values, 2);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_GENERIC);

    recomb_map_free(&recomb_map);
}

static void
verify_recomb_map(uint32_t num_loci, double length, double *positions,
        double *rates, size_t size)
{

    int ret;
    recomb_map_t recomb_map;
    double total_rate, x, y, z;
    size_t j;
    size_t num_checks = 1000;
    double eps = 1e-6;
    double *ret_rates, *ret_positions, *bulk_x;

    ret_rates = malloc(size * sizeof(double));
    ret_positions = malloc(size * sizeof(double));
    bulk_x = malloc(num_checks * sizeof(double));

    CU_ASSERT_FATAL(ret_rates != NULL);
    CU_ASSERT_FATAL(ret_positions != NULL);

    ret = recomb_map_alloc(&recomb_map, num_loci, length,
           positions, rates, size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    recomb_map_print_state(&recomb_map, _devnull);

    CU_ASSERT_EQUAL(recomb_map_get_num_loci(&recomb_map), num_loci);
    CU_ASSERT_EQUAL(recomb_map_get_size(&recomb_map), size);
    CU_ASSERT_EQUAL(recomb_map_genetic_to_phys(&recomb_map, 0), 0);
    CU_ASSERT_EQUAL(
            recomb_map_genetic_to_phys(&recomb_map, num_loci), length);
    total_rate = recomb_map_get_total_recombination_rate(&recomb_map);
    CU_ASSERT_DOUBLE_EQUAL(
            total_rate / (num_loci - 1),
            recomb_map_get_per_locus_recombination_rate(&recomb_map), eps)

    for (j = 0; j < num_checks; j++) {
        x = j * (num_loci / num_checks);
        bulk_x[j] = x;
        y = recomb_map_genetic_to_phys(&recomb_map, x);
        CU_ASSERT_TRUE(0 <= y && y <= length);
        z = recomb_map_phys_to_genetic(&recomb_map, y);
        CU_ASSERT_DOUBLE_EQUAL(x, z, eps);
    }
    ret = recomb_map_genetic_to_phys_bulk(&recomb_map, bulk_x, num_checks);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_checks; j++) {
        x = j * (num_loci / num_checks);
        y = recomb_map_genetic_to_phys(&recomb_map, x);
        CU_ASSERT_DOUBLE_EQUAL(bulk_x[j], y, 1e-12);
    }
    ret = recomb_map_get_positions(&recomb_map, ret_positions);
    CU_ASSERT_EQUAL(ret, 0);
    ret = recomb_map_get_rates(&recomb_map, ret_rates);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < size; j++) {
        CU_ASSERT_EQUAL(ret_rates[j], rates[j]);
        CU_ASSERT_EQUAL(ret_positions[j], positions[j]);
    }
    ret = recomb_map_free(&recomb_map);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    free(ret_rates);
    free(ret_positions);
    free(bulk_x);
}

static void
test_recomb_map_examples(void)
{
    double p1[] =  {0, 0.05019838393314813, 0.36933662489552865, 1};
    double r1[] = {3.5510784169955434, 4.184964179610539, 3.800808140657212, 0};
    double p2[] =  {0, 0.125, 0.875, 1, 4, 8, 16};
    double r2[] = {0.1, 6.0, 3.333, 2.1, 0.0, 2.2, 0};

    verify_recomb_map(2, 1.0, p1, r1, 4);
    verify_recomb_map(1000, 1.0, p1, r1, 4);
    verify_recomb_map(UINT32_MAX, 1.0, p1, r1, 4);

    verify_recomb_map(2, 16.0, p2, r2, 7);
    verify_recomb_map(100, 16.0, p2, r2, 7);
    verify_recomb_map(UINT32_MAX, 16.0, p2, r2, 7);
}

static void
test_single_locus_two_populations(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = {{0, 0.0}, {0, 0.0}, {1, 40.0}};
    coalescence_record_t *records;
    size_t num_records;
    uint32_t n = 3;

    CU_ASSERT_FATAL(rng != NULL);

    ret = msp_alloc(&msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_mass_migration(&msp, 30.0, 0, 1, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_mass_migration(&msp, 30.5, 1, 0, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_mass_migration(&msp, 40.5, 1, 0, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp);
    msp_print_state(&msp, _devnull);
    num_records = msp_get_num_coalescence_records(&msp);
    CU_ASSERT_EQUAL_FATAL(num_records, 2);
    ret = msp_get_coalescence_records(&msp, &records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(records[0].node, 3);
    CU_ASSERT_TRUE(records[0].time < 40.0);
    CU_ASSERT_EQUAL(records[0].population_id, 0);
    CU_ASSERT_EQUAL(records[1].node, 4);
    CU_ASSERT_TRUE(records[1].time > 40.5);
    CU_ASSERT_EQUAL(records[1].population_id, 0);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
}

static void
test_single_locus_historical_sample(void)
{
    int ret;
    msp_t msp;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    sample_t samples[] = {{0, 0.0}, {0, 10.0}};
    coalescence_record_t *record;
    size_t num_records;
    uint32_t n = 2;

    CU_ASSERT_FATAL(rng != NULL);

    ret = msp_alloc(&msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(&msp);
    msp_print_state(&msp, _devnull);
    num_records = msp_get_num_coalescence_records(&msp);
    CU_ASSERT_EQUAL_FATAL(num_records, 1);
    ret = msp_get_coalescence_records(&msp, &record);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(record->left, 0);
    CU_ASSERT_EQUAL(record->right, 1);
    CU_ASSERT_EQUAL(record->node, 2);
    CU_ASSERT_TRUE(record->time > 10.0);

    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
}

static void
test_simulator_getters_setters(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10;
    uint32_t m = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double migration_matrix[] = {0, 0, 0, 0};
    double matrix[4], growth_rate, initial_size;
    size_t migration_events[4];
    size_t breakpoints[m];
    population_t *population;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    for (j = 0; j < n; j++) {
        samples[j].time = j;
        samples[j].population_id = j % 2;
    }
    CU_ASSERT_EQUAL(msp_alloc(&msp, 0, NULL, rng), MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    samples[0].time = 1.0;
    CU_ASSERT_EQUAL(msp_alloc(&msp, n, samples, rng), MSP_ERR_BAD_SAMPLES);
    msp_free(&msp);
    samples[0].time = -1.0;
    CU_ASSERT_EQUAL(msp_alloc(&msp, n, samples, rng), MSP_ERR_BAD_PARAM_VALUE);
    msp_free(&msp);
    samples[0].time = 0.0;

    ret = msp_alloc(&msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(msp_set_max_memory(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_node_mapping_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_segment_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_avl_node_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_coalescence_record_block_size(&msp, 0),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_num_loci(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(msp_set_num_populations(&msp, 0), MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
            msp_set_scaled_recombination_rate(&msp, -1),
            MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
            msp_set_population_configuration(&msp, -1, 0, 0),
            MSP_ERR_BAD_POPULATION_ID);
    CU_ASSERT_EQUAL(
            msp_set_population_configuration(&msp, 3, 0, 0),
            MSP_ERR_BAD_POPULATION_ID);

    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
            msp_get_population_configuration(&msp, 3, NULL, NULL),
            MSP_ERR_BAD_POPULATION_ID);
    ret = msp_set_population_configuration(&msp, 0, 2, 0.5);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 0, NULL),
            MSP_ERR_BAD_MIGRATION_MATRIX);
    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 3, migration_matrix),
            MSP_ERR_BAD_MIGRATION_MATRIX);
    migration_matrix[0] = 1;
    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 4, migration_matrix),
            MSP_ERR_BAD_MIGRATION_MATRIX);
    migration_matrix[0] = 0;
    migration_matrix[1] = -1;
    CU_ASSERT_EQUAL(
            msp_set_migration_matrix(&msp, 4, migration_matrix),
            MSP_ERR_BAD_MIGRATION_MATRIX);

    migration_matrix[1] = 1;
    migration_matrix[2] = 1;
    ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(&msp, m);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(&msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_get_population_configuration(&msp, 0,
            &initial_size, &growth_rate);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(initial_size, 2);
    CU_ASSERT_EQUAL(growth_rate, 0.5);
    CU_ASSERT_EQUAL(
            msp_get_population(&msp, 3, NULL),
            MSP_ERR_BAD_POPULATION_ID);
    ret = msp_get_population(&msp, 0, &population);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(population->initial_size, 2);
    CU_ASSERT_EQUAL(population->growth_rate, 0.5);
    CU_ASSERT_EQUAL(population->start_time, 0.0);

    CU_ASSERT_EQUAL(msp_get_num_avl_node_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_node_mapping_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_segment_blocks(&msp), 1);
    CU_ASSERT_EQUAL(msp_get_num_coalescence_record_blocks(&msp), 1);
    CU_ASSERT(msp_get_used_memory(&msp) > 0);
    CU_ASSERT_EQUAL(msp_get_num_populations(&msp), 2);

    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(msp_get_num_breakpoints(&msp), m - 1);
    ret = msp_get_breakpoints(&msp, breakpoints);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < m - 1; j++) {
        CU_ASSERT_EQUAL(breakpoints[j], j + 1);
    }
    ret = msp_get_num_migration_events(&msp, migration_events);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(migration_events[0], 0);
    CU_ASSERT(migration_events[1] > 0);
    CU_ASSERT(migration_events[2] > 0);
    CU_ASSERT_EQUAL(migration_events[3], 0);
    ret = msp_get_migration_matrix(&msp, matrix);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < 4; j++) {
        CU_ASSERT_EQUAL(matrix[j], migration_matrix[j]);
    }
    CU_ASSERT(msp_get_num_common_ancestor_events(&msp) > 0);
    CU_ASSERT(msp_get_num_recombination_events(&msp) > 0);

    free(samples);
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
}

static void
test_simulator_demographic_events(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10;
    uint32_t m = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double migration_matrix[] = {0, 1, 1, 0};
    double time;
    msp_t msp;

    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    for (j = 0; j < n; j++) {
        samples[j].time = j;
        samples[j].population_id = j % 2;
    }
    ret = msp_alloc(&msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_populations(&msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 0, 1, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_population_configuration(&msp, 1, 2, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_migration_matrix(&msp, 4, migration_matrix);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(&msp, m);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(&msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, -1, 0, 1),
        MSP_ERR_BAD_POPULATION_ID);
    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, 2, 0, 1),
        MSP_ERR_BAD_POPULATION_ID);
    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, 0, 0, 1),
        MSP_ERR_SOURCE_DEST_EQUAL);
    CU_ASSERT_EQUAL(
        msp_add_mass_migration(&msp, 10, 0, 1, -5),
        MSP_ERR_BAD_PARAM_VALUE);

    CU_ASSERT_EQUAL(
        msp_add_migration_rate_change(&msp, 10, -2, 2.0),
        MSP_ERR_BAD_MIGRATION_MATRIX_INDEX);
    CU_ASSERT_EQUAL(
        msp_add_migration_rate_change(&msp, 10, -1, -2.0),
        MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
        msp_add_migration_rate_change(&msp, 10, 3, 2.0),
        MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX);

    CU_ASSERT_EQUAL(
        msp_add_population_parameters_change(&msp, 10, -2, 0, 0),
        MSP_ERR_BAD_POPULATION_ID);
    CU_ASSERT_EQUAL(
        msp_add_population_parameters_change(&msp, 10, -1, -1, 0),
        MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
        msp_add_population_parameters_change(&msp, 10, -1, GSL_NAN, GSL_NAN),
        MSP_ERR_BAD_PARAM_VALUE);

    CU_ASSERT_EQUAL(
        msp_add_bottleneck(&msp, 10, -1, 0),
        MSP_ERR_BAD_POPULATION_ID);
    CU_ASSERT_EQUAL(
        msp_add_bottleneck(&msp, 10, 0, -1),
        MSP_ERR_BAD_PARAM_VALUE);
    CU_ASSERT_EQUAL(
        msp_add_bottleneck(&msp, 10, 0, 1.1),
        MSP_ERR_BAD_PARAM_VALUE);

    ret = msp_add_mass_migration(&msp, 0.1, 0, 1, 0.5);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_migration_rate_change(&msp, 0.2, 1, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_migration_rate_change(&msp, 0.3, -1, 3.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.4, 0, 0.5, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.5, -1, 0.5, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.6, 0, GSL_NAN, 2.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_population_parameters_change(&msp, 0.7, 1, 1, GSL_NAN);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_add_bottleneck(&msp, 0.8, 0, 1.0);
    CU_ASSERT_EQUAL(ret, 0);

    CU_ASSERT_EQUAL(
            msp_add_mass_migration(&msp, 0.1, 0, 1, 0.5),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
    CU_ASSERT_EQUAL(
            msp_add_migration_rate_change(&msp, 0.2, 1, 2.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
    CU_ASSERT_EQUAL(
            msp_add_population_parameters_change(&msp, 0.4, 0, 0.5, 1.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);
    CU_ASSERT_EQUAL(
            msp_add_bottleneck(&msp, 0.7, 0, 1.0),
            MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS);

    CU_ASSERT_EQUAL(
        msp_debug_demography(&msp, &time),
        MSP_ERR_BAD_STATE);

    ret = msp_initialise(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    j = 0;
    do {
        ret = msp_debug_demography(&msp, &time);
        CU_ASSERT_EQUAL(ret, 0);
        msp_print_state(&msp, _devnull);
        j++;
    } while (! gsl_isinf(time));
    CU_ASSERT_EQUAL(j, 9);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
        msp_run(&msp, DBL_MAX, ULONG_MAX),
        MSP_ERR_BAD_STATE);
    ret = msp_reset(&msp);
    CU_ASSERT_EQUAL(ret, 0);

    msp_print_state(&msp, _devnull);
    ret = msp_run(&msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);

    free(samples);
    ret = msp_free(&msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
}

static void
test_single_locus_simulation(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    /* For the single locus sim we should have exactly n - 1 events */
    for (j = 0; j < n - 2; j++) {
        ret = msp_run(msp, DBL_MAX, 1);
        CU_ASSERT_EQUAL(ret, 1);
        msp_verify(msp);
    }
    ret = msp_run(msp, DBL_MAX, 1);
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
}

static void
test_simulation_memory_limit(void)
{
    int ret;
    uint32_t n = 1000;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_max_memory(msp, 1024 * 1024);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(msp, 10000);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, MSP_ERR_NO_MEMORY);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
}

static void
test_multi_locus_simulation(void)
{
    int ret;
    uint32_t num_events;
    uint32_t n = 100;
    uint32_t m = 100;
    long seed = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    gsl_rng_set(rng, seed);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    /* set all the block sizes to something small to provoke the memory
     * expansions. */
    ret = msp_set_avl_node_block_size(msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_node_mapping_block_size(msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_segment_block_size(msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_coalescence_record_block_size(msp, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(msp, m);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    num_events = 0;
    while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
        msp_verify(msp);
        num_events++;
    }
    CU_ASSERT_EQUAL(ret, 0);
    msp_verify(msp);
    CU_ASSERT(num_events > n - 1);
    CU_ASSERT_EQUAL(1 + num_events,
            msp_get_num_recombination_events(msp) +
            msp_get_num_common_ancestor_events(msp))

    gsl_rng_set(rng, seed);
    ret = msp_reset(msp);
    num_events = 0;
    while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
        msp_verify(msp);
        num_events++;
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(1 + num_events,
            msp_get_num_recombination_events(msp) +
            msp_get_num_common_ancestor_events(msp))

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
}

static void
test_bottleneck_simulation(void)
{
    int ret;
    uint32_t n = 100;
    uint32_t m = 100;
    long seed = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double t1 = 0.1;
    double t2 = 0.5;
    int t1_found = 0;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    gsl_rng_set(rng, seed);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    /* set all the block sizes to something small to provoke the memory
     * expansions. */
    ret = msp_set_avl_node_block_size(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_node_mapping_block_size(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_segment_block_size(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_coalescence_record_block_size(msp, 2);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(msp, m);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    /* Add a bottleneck that does nothing at time t1 */
    ret = msp_add_bottleneck(msp, t1, 0, 0);
    CU_ASSERT_EQUAL(ret, 0);
    /* Add a bottleneck that coalesces everything at t2 */
    ret = msp_add_bottleneck(msp, t2, 0, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    ret = msp_run(msp, t1, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 2);
    CU_ASSERT_FALSE(msp_is_completed(msp));
    msp_print_state(msp, _devnull);

    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(msp));
    CU_ASSERT_EQUAL(msp->time, t2);
    msp_verify(msp);

    msp_reset(msp);
    while ((ret = msp_run(msp, DBL_MAX, 1)) == 1) {
        msp_verify(msp);
        if (msp->time == 0.1) {
            t1_found = 1;
        }
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(t1_found);
    CU_ASSERT_EQUAL(msp->time, t2);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
}

static void
test_large_bottleneck_simulation(void)
{
    int ret;
    uint32_t j;
    uint32_t n = 10000;
    uint32_t m = 100;
    long seed = 10;
    sample_t *samples = malloc(n * sizeof(sample_t));
    msp_t *msp = malloc(sizeof(msp_t));
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    uint32_t num_bottlenecks = 10;
    bottleneck_desc_t bottlenecks[num_bottlenecks];
    double t;

    CU_ASSERT_FATAL(msp != NULL);
    CU_ASSERT_FATAL(samples != NULL);
    CU_ASSERT_FATAL(rng != NULL);

    t = 0.1;
    for (j = 0; j < num_bottlenecks; j++) {
        bottlenecks[j].time = t;
        bottlenecks[j].intensity = 0.1;
        t += 0.01;
    }
    /* Set the last bottleneck to be full intensity */
    bottlenecks[num_bottlenecks - 1].intensity = 1.0;

    gsl_rng_set(rng, seed);
    memset(samples, 0, n * sizeof(sample_t));
    ret = msp_alloc(msp, n, samples, rng);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_scaled_recombination_rate(msp, 1.0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = msp_set_num_loci(msp, m);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_bottlenecks; j++) {
        msp_add_bottleneck(msp, bottlenecks[j].time, 0,
                bottlenecks[j].intensity);
        CU_ASSERT_EQUAL(ret, 0);
    }
    ret = msp_initialise(msp);
    CU_ASSERT_EQUAL(ret, 0);

    for (j = 0; j < num_bottlenecks - 1; j++) {
        ret = msp_run(msp, bottlenecks[j].time, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 2);
        CU_ASSERT_FALSE(msp_is_completed(msp));
        CU_ASSERT_EQUAL(msp->time, bottlenecks[j].time);
        msp_verify(msp);
    }
    ret = msp_run(msp, DBL_MAX, ULONG_MAX);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(msp_is_completed(msp));
    CU_ASSERT_EQUAL(msp->time, bottlenecks[num_bottlenecks - 1].time);
    msp_verify(msp);

    /* Test out resets on partially completed simulations. */
    ret = msp_reset(msp);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_bottlenecks - 1; j++) {
        ret = msp_run(msp, bottlenecks[j].time, ULONG_MAX);
        CU_ASSERT_EQUAL(ret, 2);
    }
    ret = msp_reset(msp);
    msp_verify(msp);

    ret = msp_free(msp);
    CU_ASSERT_EQUAL(ret, 0);
    gsl_rng_free(rng);
    free(msp);
    free(samples);
}

static void
test_simplest_records(void)
{
    int ret;
    uint32_t children[] = {0, 1};
    coalescence_record_t record;
    tree_sequence_t ts;

    record.left = 0.0;
    record.right = 1.0;
    record.node = 2;
    record.num_children = 2;
    record.children = children;
    record.population_id = 0;
    record.time = 1.0;

    ret = tree_sequence_load_records(&ts, 1, &record);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 2);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 3);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    tree_sequence_free(&ts);
}

static void
test_simplest_nonbinary_records(void)
{
    int ret;
    uint32_t children[] = {0, 1, 2, 3};
    coalescence_record_t record;
    tree_sequence_t ts;

    record.left = 0.0;
    record.right = 1.0;
    record.node = 4;
    record.num_children = 4;
    record.children = children;
    record.population_id = 0;
    record.time = 1.0;

    ret = tree_sequence_load_records(&ts, 1, &record);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 5);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    tree_sequence_free(&ts);
}

static void
test_simplest_bad_records(void)
{
    int ret = 0;
    uint32_t children[] = {0, 1};
    coalescence_record_t records[1];
    size_t num_records = 1;
    uint32_t j;
    tree_sequence_t ts;

    records[0].left = 0.0;
    records[0].right = 1.0;
    records[0].node = 2;
    records[0].num_children = 2;
    records[0].children = children;
    records[0].population_id = 0;
    records[0].time = 1.0;

    /* Make sure we have a good set of records */
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    tree_sequence_free(&ts);

    /* An empty sequence should be an error */
    ret = tree_sequence_load_records(&ts, 0, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);

    /* Bad sequence length */
    records[0].right = 0.0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].right = 1.0;

    /* Equal nodes in the children */
    records[0].children[0] = 1;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].children[0] = 0;

    /* children node >= parent */
    records[0].children[0] = 2;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].children[0] = 0;

    /* Unsorted nodes in the children */
    records[0].children[0] = 1;
    records[0].children[1] = 0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].children[0] = 0;
    records[0].children[1] = 1;

    /* Null parent */
    records[0].node = MSP_NULL_NODE;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].node = 2;

    /* Null child */
    records[0].children[0] = MSP_NULL_NODE;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].children[0] = 0;

    /* 0 or 1 children */
    for (j = 0; j < 2; j++) {
        records[0].num_children = j;
        ret = tree_sequence_load_records(&ts, num_records, records);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
        tree_sequence_free(&ts);
    }
    records[0].num_children = 2;

    /* Make sure we've preserved a good tree sequence */
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    tree_sequence_free(&ts);
}

static void
test_single_tree_good_records(void)
{
    int ret = 0;
    const char * text_records[] = {
        "0 1 4 0,1 1.0 0",
        "0 1 5 2,3 2.0 0",
        "0 1 6 4,5 3.0 0"
    };
    coalescence_record_t *records = NULL;
    size_t num_records = 3;
    tree_sequence_t ts;

    parse_text_records(num_records, text_records, &records);
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 7);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

static void
test_single_nonbinary_tree_good_records(void)
{
    int ret = 0;
    const char * text_records[] = {
        "0 1 7 0,1,2,3 1.0 0",
        "0 1 8 4,5     1.0 0",
        "0 1 9 6,7,8   1.0 0",
    };
    coalescence_record_t *records = NULL;
    size_t num_records = 3;
    tree_sequence_t ts;

    parse_text_records(num_records, text_records, &records);
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 7);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 1.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 10);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

static void
test_single_tree_bad_records(void)
{
    int ret = 0;
    const char * text_records[] = {
        "0 1 4 0,1 1.0 0",
        "0 1 5 2,3 2.0 0",
        "0 1 6 4,5 3.0 0"
    };
    coalescence_record_t *records = NULL;
    size_t num_records = 3;
    tree_sequence_t ts;

    parse_text_records(num_records, text_records, &records);

    /* Not sorted in time order */
    records[2].time = 0.5;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[2].time = 3;

    /* Left value greater than sequence length */
    records[2].left = 2.0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[2].left = 0.0;

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    tree_sequence_free(&ts);

    free_local_records(num_records, records);
}

static void
test_single_tree_good_mutations(void)
{
    int ret = 0;
    const char * text_records[] = {
        "0 6 4 0,1 1.0 0",
        "0 6 5 2,3 2.0 0",
        "0 6 6 4,5 3.0 0"
    };
    coalescence_record_t *records = NULL;
    mutation_t *mutations = NULL;
    mutation_t *other_mutations;
    size_t num_mutations = 6;
    size_t num_records = 3;
    size_t j;
    tree_sequence_t ts;

    parse_text_records(num_records, text_records, &records);

    mutations = malloc(num_mutations * sizeof(mutation_t));
    CU_ASSERT_FATAL(mutations != NULL);
    for (j = 0; j < num_mutations; j++) {
        mutations[j].position = (double) j;
        mutations[j].node = (uint32_t) j;
    }
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_sample_size(&ts), 4);
    CU_ASSERT_EQUAL(tree_sequence_get_sequence_length(&ts), 6.0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), 7);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);

    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), num_mutations);
    ret = tree_sequence_get_mutations(&ts, &other_mutations);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_mutations; j++) {
        CU_ASSERT_EQUAL(other_mutations[j].index, j);
        CU_ASSERT_EQUAL(mutations[j].position, other_mutations[j].position);
        CU_ASSERT_EQUAL(mutations[j].node, other_mutations[j].node);
    }
    free(mutations);
    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

static void
test_single_tree_bad_mutations(void)
{
    int ret = 0;
    const char * text_records[] = {
        "0 1 4 0,1 1.0 0",
        "0 1 5 2,3 2.0 0",
        "0 1 6 4,5 3.0 0"
    };
    coalescence_record_t *records = NULL;
    mutation_t mutations[] = {{0, 0}, {0, 1}};
    size_t num_mutations = 2;
    size_t num_records = 3;
    tree_sequence_t ts;

    parse_text_records(num_records, text_records, &records);

    /* negative coordinate */
    mutations[0].position = -1.0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_MUTATION);
    tree_sequence_free(&ts);
    mutations[0].position = 0.0;

    /* coordinate > sequence length */
    mutations[0].position = 1.1;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_MUTATION);
    tree_sequence_free(&ts);
    mutations[0].position = 0.0;

    /* node = NULL */
    mutations[0].node = MSP_NULL_NODE;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_MUTATION);
    tree_sequence_free(&ts);
    mutations[0].node = 0;

    /* node >= num_nodes */
    mutations[0].node = 7;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_MUTATION);
    tree_sequence_free(&ts);
    mutations[0].node = 0;

    /* Check to make sure we've maintained legal mutations */
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), num_mutations);
    tree_sequence_free(&ts);

    free_local_records(num_records, records);
}

static void
test_single_tree_iter(void)
{
    int ret = 0;
    const char * text_records[] = {
        "0 1 4 0,1 1.0 0",
        "0 1 5 2,3 2.0 0",
        "0 1 6 4,5 3.0 0"
    };
    coalescence_record_t *records = NULL;
    uint32_t parents[] = {4, 4, 5, 5, 6, 6, MSP_NULL_NODE};
    size_t num_records = 3;
    uint32_t num_nodes = 7;
    uint32_t u, v, num_leaves, w;
    tree_sequence_t ts;
    sparse_tree_t tree;

    parse_text_records(num_records, text_records, &records);

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    sparse_tree_print_state(&tree, _devnull);

    for (u = 0; u < num_nodes; u++) {
        ret = sparse_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
    }
    ret = sparse_tree_get_num_leaves(&tree, 0, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 1);
    ret = sparse_tree_get_num_leaves(&tree, 4, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 2);
    ret = sparse_tree_get_num_leaves(&tree, 6, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 4);
    ret = sparse_tree_get_mrca(&tree, 0, 1, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 4);
    ret = sparse_tree_get_mrca(&tree, 0, 2, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 6);

    ret = sparse_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

static void
test_single_nonbinary_tree_iter(void)
{
    int ret = 0;
    const char * text_records[] = {
        "0 1 7 0,1,2,3 1.0 0",
        "0 1 8 4,5     1.0 0",
        "0 1 9 6,7,8   1.0 0",
    };
    coalescence_record_t *records = NULL;
    uint32_t parents[] = {7, 7, 7, 7, 8, 8, 9, 9, 9, MSP_NULL_NODE};
    size_t num_records = 3;
    uint32_t num_nodes = 10;
    uint32_t sample_size = 7;
    uint32_t u, v, num_leaves, num_children, w, *children;
    tree_sequence_t ts;
    sparse_tree_t tree;

    parse_text_records(num_records, text_records, &records);

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    sparse_tree_print_state(&tree, _devnull);

    for (u = 0; u < num_nodes; u++) {
        ret = sparse_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
    }
    for (u = 0; u < sample_size; u++) {
        ret = sparse_tree_get_num_leaves(&tree, u, &num_leaves);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(num_leaves, 1);
        ret = sparse_tree_get_children(&tree, u, &num_children, &children);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(children, NULL);
        CU_ASSERT_EQUAL(num_children, 0);
    }

    u = 7;
    ret = sparse_tree_get_num_leaves(&tree, u, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 4);
    ret = sparse_tree_get_children(&tree, u, &num_children, &children);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_children, 4);
    CU_ASSERT_EQUAL(children[0], 0);
    CU_ASSERT_EQUAL(children[1], 1);
    CU_ASSERT_EQUAL(children[2], 2);
    CU_ASSERT_EQUAL(children[3], 3);

    u = 8;
    ret = sparse_tree_get_num_leaves(&tree, u, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 2);
    ret = sparse_tree_get_children(&tree, u, &num_children, &children);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_children, 2);
    CU_ASSERT_EQUAL(children[0], 4);
    CU_ASSERT_EQUAL(children[1], 5);

    u = 9;
    ret = sparse_tree_get_num_leaves(&tree, u, &num_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_leaves, 7);
    ret = sparse_tree_get_children(&tree, u, &num_children, &children);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(num_children, 3);
    CU_ASSERT_EQUAL(children[0], 6);
    CU_ASSERT_EQUAL(children[1], 7);
    CU_ASSERT_EQUAL(children[2], 8);

    ret = sparse_tree_get_root(&tree, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 9);

    ret = sparse_tree_get_mrca(&tree, 0, 1, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 7);
    ret = sparse_tree_get_mrca(&tree, 0, 4, &w);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(w, 9);

    ret = sparse_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

static void
test_single_tree_iter_times(void)
{
    int ret = 0;
    const char * text_records[] = {
        "0 1 4 0,1 1.0 0",
        "0 1 5 2,3 4.0 0",
        "0 1 6 4,5 5.0 0"
    };
    /* 0 and 1 sampled at time 0 and 2 and 3 later. */
    sample_t samples[] = {
        {0, 0.0},
        {0, 0.0},
        {0, 2.0},
        {0, 3.0}};
    size_t num_records = 3;
    size_t sample_size = sizeof(samples) / sizeof(sample_t);
    uint32_t parents[] = {4, 4, 5, 5, 6, 6, MSP_NULL_NODE};
    double times[] = {0.0, 0.0, 2.0, 3.0, 1.0, 4.0, 5.0};
    double t;
    uint32_t u, v;
    uint32_t num_nodes = 7;
    tree_sequence_t ts;
    sparse_tree_t tree;
    coalescence_record_t *records;

    parse_text_records(num_records, text_records, &records);

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_set_samples(&ts, sample_size, samples);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, 1);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), 1);
    sparse_tree_print_state(&tree, _devnull);

    for (u = 0; u < num_nodes; u++) {
        ret = sparse_tree_get_parent(&tree, u, &v);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(v, parents[u]);
        ret = sparse_tree_get_time(&tree, u, &t);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(t, times[u]);
    }
    ret = sparse_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

static void
test_single_tree_hapgen(void)
{
    int ret = 0;
    const char * text_records[] = {
        "0 1 4 0,1 1.0 0",
        "0 1 5 2,3 4.0 0",
        "0 1 6 4,5 5.0 0"
    };
    mutation_t mutations[] = {{0, 0}, {0, 1}};
    size_t sample_size = 4;
    size_t num_mutations = 2;
    char *haplotype;
    size_t j;
    size_t num_records = 3;
    tree_sequence_t ts;
    coalescence_record_t *records;
    hapgen_t hapgen;

    parse_text_records(num_records, text_records, &records);
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);
    for (j = 0; j < sample_size; j++) {
        ret = hapgen_get_haplotype(&hapgen, j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_STRING_EQUAL(haplotype, "");
    }
    ret = hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = hapgen_alloc(&hapgen, &ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);

    ret = hapgen_get_haplotype(&hapgen, 0, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "10");
    ret = hapgen_get_haplotype(&hapgen, 1, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "01");
    ret = hapgen_get_haplotype(&hapgen, 2, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "00");
    ret = hapgen_get_haplotype(&hapgen, 3, &haplotype);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_STRING_EQUAL(haplotype, "00");

    ret = hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

static void
test_single_tree_mutgen(void)
{
    int ret = 0;
    const char * text_records[] = {
        "0 1 4 0,1 1.0 0",
        "0 1 5 2,3 4.0 0",
        "0 1 6 4,5 5.0 0"
    };
    mutation_t *mutations, *before, *after;

    size_t j, num_mutations;
    size_t num_records = 3;
    tree_sequence_t ts;
    coalescence_record_t *records;
    mutgen_t mutgen;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

    parse_text_records(num_records, text_records, &records);
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(rng != NULL);

    ret = mutgen_alloc(&mutgen, &ts, 0.0, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(mutgen_get_num_mutations(&mutgen), 0);
    ret = mutgen_generate(&mutgen);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(mutgen_get_num_mutations(&mutgen), 0);
    mutgen_print_state(&mutgen, _devnull);
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    gsl_rng_set(rng, 1);
    ret = mutgen_alloc(&mutgen, &ts, 10.0, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(mutgen_get_num_mutations(&mutgen), 0);
    ret = mutgen_generate(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    mutgen_print_state(&mutgen, _devnull);
    num_mutations = mutgen_get_num_mutations(&mutgen);
    CU_ASSERT_TRUE(num_mutations > 0);
    ret = mutgen_get_mutations(&mutgen, &mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < num_mutations; j++) {
        CU_ASSERT_TRUE(mutations[j].position >= 0.0);
        CU_ASSERT_TRUE(mutations[j].position <= 1.0);
        CU_ASSERT_TRUE(mutations[j].node < 6);
    }
    before = malloc(num_mutations * sizeof(mutation_t));
    CU_ASSERT_FATAL(before != NULL);
    memcpy(before, mutations, num_mutations * sizeof(mutation_t));
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Test the reallocing behavior by setting a very small
     * block size.
     */
    gsl_rng_set(rng, 1);
    ret = mutgen_alloc(&mutgen, &ts, 10.0, rng);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(mutgen_get_num_mutations(&mutgen), 0);
    ret = mutgen_set_mutation_block_size(&mutgen, 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = mutgen_set_mutation_block_size(&mutgen, 1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = mutgen_generate(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(mutgen_get_num_mutations(&mutgen), num_mutations);
    ret = mutgen_get_mutations(&mutgen, &mutations);
    after = malloc(num_mutations * sizeof(mutation_t));
    CU_ASSERT_FATAL(after != NULL);
    memcpy(after, mutations, num_mutations * sizeof(mutation_t));
    ret = mutgen_free(&mutgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* before and after should be identical */
    for (j = 0; j < num_mutations; j++) {
        CU_ASSERT_EQUAL(before[j].position, after[j].position);
        CU_ASSERT_EQUAL(before[j].node, after[j].node);
    }
    free(before);
    free(after);
    tree_sequence_free(&ts);
    free_local_records(num_records, records);
    gsl_rng_free(rng);
}

static void
verify_trees(size_t num_records, coalescence_record_t *records,
        size_t num_trees, size_t num_nodes, uint32_t* parents,
        size_t num_mutations, mutation_t *mutations)
{
    int ret;
    uint32_t u, v, j, k, mutation_index;
    tree_sequence_t ts;
    sparse_tree_t tree;
    mutation_t *tree_mutations;
    size_t num_tree_mutations;

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_nodes(&ts), num_nodes);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(&ts), num_trees);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), num_mutations);
    ret = sparse_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    mutation_index = 0;
    j = 0;
    for (ret = sparse_tree_first(&tree); ret == 1;
            ret = sparse_tree_next(&tree)) {
        CU_ASSERT_EQUAL(j, tree.index);
        sparse_tree_print_state(&tree, _devnull);
        for (u = 0; u < num_nodes; u++) {
            ret = sparse_tree_get_parent(&tree, u, &v);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(v, parents[j * num_nodes + u]);
        }
        ret = sparse_tree_get_mutations(&tree, &num_tree_mutations,
                &tree_mutations);
        CU_ASSERT_EQUAL(ret, 0);
        for (k = 0; k < num_tree_mutations; k++) {
            CU_ASSERT_EQUAL(tree_mutations[k].index, mutation_index);
            CU_ASSERT_EQUAL(
                tree_mutations[k].position, mutations[mutation_index].position);
            CU_ASSERT_EQUAL(
                tree_mutations[k].node, mutations[mutation_index].node);
            mutation_index++;
        }
        j++;
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(mutation_index, num_mutations);

    ret = sparse_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);

    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
}

static void
verify_trees_consistent(tree_sequence_t *ts)
{
    int ret, found;
    size_t sample_size, num_trees;
    uint32_t u, v, j, root, k, num_children, *children;
    sparse_tree_t tree;

    sample_size = tree_sequence_get_sample_size(ts);
    ret = sparse_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    num_trees = 0;
    for (ret = sparse_tree_first(&tree); ret == 1;
            ret = sparse_tree_next(&tree)) {
        ret = sparse_tree_get_root(&tree, &root);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(tree.index, num_trees);
        num_trees++;
        for (j = 0; j < sample_size; j++) {
            v = j;
            while (v != MSP_NULL_NODE) {
                u = v;
                ret = sparse_tree_get_parent(&tree, u, &v);
                CU_ASSERT_EQUAL(ret, 0);
                if (v != MSP_NULL_NODE) {
                    ret = sparse_tree_get_children(&tree, v,
                            &num_children, &children);
                    CU_ASSERT_EQUAL(ret, 0);
                    CU_ASSERT(num_children >= 2);
                    found = 0;
                    for (k = 0; k < num_children; k++) {
                        if (children[k] == u) {
                            found = 1;
                        }
                    }
                    CU_ASSERT(found);
                }
            }
            CU_ASSERT_EQUAL(u, root);
        }
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_trees(ts), num_trees);

    sparse_tree_free(&tree);
}

static void
verify_tree_iter_fails(size_t num_records, coalescence_record_t *records,
        size_t num_mutations, mutation_t *mutations,
        uint32_t tree_index, int error_code)
{
    int ret;
    uint32_t index = 0;
    tree_sequence_t ts;
    sparse_tree_t tree;

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    for (ret = sparse_tree_first(&tree); ret == 1;
            ret = sparse_tree_next(&tree)) {
        index++;
    }
    CU_ASSERT_EQUAL(index, tree_index);
    CU_ASSERT_EQUAL(ret, error_code);

    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
}

static void
verify_diff_iter_fails(size_t num_records, coalescence_record_t *records,
        uint32_t tree_index, int error_code)
{
    int ret;
    uint32_t index = 0;
    tree_sequence_t ts;
    tree_diff_iterator_t iter;
    node_record_t *records_out, *records_in;
    double length;

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_diff_iterator_alloc(&iter, &ts);
    CU_ASSERT_EQUAL(ret, 0);

    while ((ret = tree_diff_iterator_next(
            &iter, &length, &records_out, &records_in)) == 1) {
        index++;
    }
    CU_ASSERT_EQUAL(index, tree_index);
    CU_ASSERT_EQUAL(ret, error_code);

    tree_diff_iterator_free(&iter);
    tree_sequence_free(&ts);
}

static void
test_sparse_tree_errors(void)
{
    int ret;
    const char * text_records[] = {
        "2 10 4 2,3 0.071 0",
        "0 2  5 1,3 0.090 0",
        "2 10 5 1,4 0.090 0",
        "0 7  6 0,5 0.170 0",
        "7 10 7 0,5 0.202 0",
        "0 2  8 2,6 0.253 0"
    };
    size_t j;
    size_t num_records = 6;
    uint32_t num_nodes = 9;
    uint32_t u;
    coalescence_record_t *records;
    tree_sequence_t ts, other_ts;
    sparse_tree_t t, other_t;
    uint32_t bad_nodes[] = {num_nodes, num_nodes + 1, UINT32_MAX};
    uint32_t tracked_leaves[] = {0, 0, 0};

    parse_text_records(num_records, text_records, &records);
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_load_records(&other_ts, num_records, records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = sparse_tree_alloc(&t, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = sparse_tree_alloc(&t, &ts, MSP_LEAF_COUNTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 1);

    /* Out-of-bounds queries */
    for (j = 0; j < sizeof(bad_nodes) / sizeof(uint32_t); j++) {
        u = bad_nodes[j];
        ret = sparse_tree_get_children(&t, u, NULL, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_parent(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_time(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_mrca(&t, u, 0, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_mrca(&t, 0, u, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_num_leaves(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_num_tracked_leaves(&t, u, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = sparse_tree_get_leaf_list(&t, u, NULL, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    }

    tracked_leaves[0] = 0;
    tracked_leaves[1] = tree_sequence_get_sample_size(&ts);
    ret = sparse_tree_set_tracked_leaves(&t, 2, tracked_leaves);
    CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    tracked_leaves[1] = 0;
    ret = sparse_tree_set_tracked_leaves(&t, 2, tracked_leaves);
    CU_ASSERT_EQUAL(ret, MSP_ERR_DUPLICATE_TRACKED_LEAF);
    ret = sparse_tree_set_tracked_leaves_from_leaf_list(&t, NULL, NULL);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    ret = sparse_tree_alloc(&other_t, &other_ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_copy(&t, &t);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    ret = sparse_tree_copy(&t, &other_t);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);

    /* TODO run checks for the various unsupported operations with flags */

    sparse_tree_free(&t);
    sparse_tree_free(&other_t);
    tree_sequence_free(&other_ts);
    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

static void
test_tree_sequence_iter(void)
{
    const char * text_records[] = {
        "2 10 4 2,3 0.071 0",
        "0 2  5 1,3 0.090 0",
        "2 10 5 1,4 0.090 0",
        "0 7  6 0,5 0.170 0",
        "7 10 7 0,5 0.202 0",
        "0 2  8 2,6 0.253 0"
    };
    /* We make one mutation for each tree */
    mutation_t mutations[] = {{1, 2}, {4.5, 0}, {8.5, 5}};
    uint32_t parents[] = {
        6, 5, 8, 5, MSP_NULL_NODE, 6, 8, MSP_NULL_NODE, MSP_NULL_NODE,
        6, 5, 4, 4, 5, 6, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
        7, 5, 4, 4, 5, 7, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
    };
    size_t num_records = 6;
    uint32_t num_nodes = 9;
    uint32_t num_trees = 3;
    uint32_t num_mutations = 3;
    coalescence_record_t *records;

    parse_text_records(num_records, text_records, &records);
    verify_trees(num_records, records, num_trees, num_nodes, parents,
            num_mutations, mutations);
    free_local_records(num_records, records);
}

static void
test_nonbinary_tree_sequence_iter(void)
{
    const char * text_records[] = {
        "0	100	8	0,1,2,3	0.01	0",
        "0	100	9	6,8	    0.068   0",
        "0	17	10	4,5,7	0.2	    0",
        "17	100	10	4,7	    0.2	    0",
        "17	100	11	5,9	    0.279   0",
        "0	17	12	9,10	0.405   0",
        "17	100	12	10,11	0.405   0",
    };
    /* We make one mutation for each tree */
    mutation_t mutations[] = {{1, 2}, {18, 11}};
    uint32_t parents[] = {
        8, 8, 8, 8, 10, 10, 9, 10, 9, 12, 12, MSP_NULL_NODE, MSP_NULL_NODE,
        8, 8, 8, 8, 10, 11, 9, 10, 9, 11, 12, 12, MSP_NULL_NODE,
    };
    uint32_t num_nodes = 13;
    uint32_t num_trees = 2;
    uint32_t num_mutations = 2;
    size_t num_records = sizeof(text_records) / sizeof(char *);
    coalescence_record_t *records;

    parse_text_records(num_records, text_records, &records);
    verify_trees(num_records, records, num_trees, num_nodes, parents,
            num_mutations, mutations);
    free_local_records(num_records, records);
}

static void
test_left_to_right_tree_sequence_iter(void)
{
    const char * text_records[] = {
        "2 10 7 2,3 0.071 0",
        "0 2  4 1,3 0.090 0",
        "2 10 4 1,7 0.090 0",
        "0 7  5 0,4 0.170 0",
        "7 10 8 0,4 0.202 0",
        "0 2  6 2,5 0.253 0"
    };
    /* We make one mutation for each tree */
    mutation_t mutations[] = {{1, 2}, {4.5, 0}, {8.5, 4}};
    uint32_t parents[] = {
        5, 4, 6, 4, 5, 6, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
        5, 4, 7, 7, 5, MSP_NULL_NODE, MSP_NULL_NODE, 4, MSP_NULL_NODE,
        8, 4, 7, 7, 8, MSP_NULL_NODE, MSP_NULL_NODE, 4, MSP_NULL_NODE,
    };
    size_t num_records = 6;
    uint32_t num_nodes = 9;
    uint32_t num_trees = 3;
    uint32_t num_mutations = 3;
    coalescence_record_t *records;

    parse_text_records(num_records, text_records, &records);
    verify_trees(num_records, records, num_trees, num_nodes, parents,
            num_mutations, mutations);
    free_local_records(num_records, records);
}

static void
test_tree_sequence_set_mutations(void)
{
    const char * text_records[] = {
        "2 10 4 2,3 0.071 0",
        "0 2  5 1,3 0.090 0",
        "2 10 5 1,4 0.090 0",
        "0 7  6 0,5 0.170 0",
        "7 10 7 0,5 0.202 0",
        "0 2  8 2,6 0.253 0"
    };
    /* We make one mutation for each tree */
    mutation_t mutations[] = {{1, 2}, {4.5, 0}, {8.5, 5}};
    tree_sequence_t ts;
    size_t j;
    const size_t num_trees = 3;
    const size_t num_records = 6;
    uint32_t num_mutations = 3;
    coalescence_record_t *records;
    sparse_tree_t trees[num_trees];
    int ret;

    parse_text_records(num_records, text_records, &records);

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tree_sequence_set_mutations(&ts, 1, mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 1);
    ret = tree_sequence_set_mutations(&ts, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 0);
    ret = tree_sequence_set_mutations(&ts, num_mutations, mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), num_mutations);

    for (j = 0; j < num_trees; j++) {
        ret = sparse_tree_alloc(&trees[j], &ts, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(ts.refcount, j + 1);
        ret = tree_sequence_set_mutations(&ts, 1, mutations);
        CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_REFCOUNT_NONZERO);

    }
    for (j = 0; j < num_trees; j++) {
        ret = tree_sequence_set_mutations(&ts, 1, mutations);
        CU_ASSERT_EQUAL_FATAL(ret, MSP_ERR_REFCOUNT_NONZERO);
        ret = sparse_tree_free(&trees[j]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(ts.refcount, num_trees - j - 1);
    }
    ret = tree_sequence_set_mutations(&ts, 1, mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tree_sequence_get_num_mutations(&ts), 1);

    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

typedef struct {
    uint32_t tree_index;
    uint32_t node;
    uint32_t count;
} leaf_count_test_t;

static void
verify_leaf_counts(tree_sequence_t *ts, size_t num_tests,
        leaf_count_test_t *tests)
{
    int ret;
    uint32_t j, num_leaves, n, k;
    uint32_t *tracked_leaves = NULL;
    sparse_tree_t tree;
    leaf_list_node_t *u, *head, *tail;

    n = tree_sequence_get_sample_size(ts);

    /* First run without the MSP_LEAF_COUNTS feature */
    ret = sparse_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = sparse_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = sparse_tree_get_num_leaves(&tree, tests[j].node, &num_leaves);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_leaves);
        /* all operations depending on tracked leaves should fail. */
        ret = sparse_tree_get_num_tracked_leaves(&tree, 0, &num_leaves);
        CU_ASSERT_EQUAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);
        ret = sparse_tree_get_leaf_list(&tree, 0, NULL, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);
    }
    sparse_tree_free(&tree);

    /* Now run with MSP_LEAF_COUNTS but with no leaves tracked. */
    ret = sparse_tree_alloc(&tree, ts, MSP_LEAF_COUNTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = sparse_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = sparse_tree_get_num_leaves(&tree, tests[j].node, &num_leaves);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_leaves);
        /* all operations depending on tracked leaves should fail. */
        ret = sparse_tree_get_num_tracked_leaves(&tree, 0, &num_leaves);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(num_leaves, 0);
        /* Getting leaf lists should still fail, as it's not enabled. */
        ret = sparse_tree_get_leaf_list(&tree, 0, NULL, NULL);
        CU_ASSERT_EQUAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);
    }
    sparse_tree_free(&tree);

    /* Run with MSP_LEAF_LISTS, but without MSP_LEAF_COUNTS */
    ret = sparse_tree_alloc(&tree, ts, MSP_LEAF_LISTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = sparse_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = sparse_tree_get_num_leaves(&tree, tests[j].node, &num_leaves);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_leaves);
        /* all operations depending on tracked leaves should fail. */
        ret = sparse_tree_get_num_tracked_leaves(&tree, 0, &num_leaves);
        CU_ASSERT_EQUAL(ret, MSP_ERR_UNSUPPORTED_OPERATION);
        /* Getting leaf lists should still fail, as it's not enabled. */
        ret = sparse_tree_get_leaf_list(&tree, tests[j].node, &head, &tail);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        u = head;
        k = 0;
        while (1) {
            k++;
            if (u == tail) {
                break;
            }
            u = u->next;
        }
        CU_ASSERT_EQUAL(tests[j].count, k);
    }
    sparse_tree_free(&tree);

    /* Now use MSP_LEAF_COUNTS|MSP_LEAF_LISTS */
    tracked_leaves = malloc(n * sizeof(uint32_t));
    for (j = 0; j < n; j++) {
        tracked_leaves[j] = j;
    }
    ret = sparse_tree_alloc(&tree, ts, MSP_LEAF_COUNTS|MSP_LEAF_LISTS);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_set_tracked_leaves(&tree, n, tracked_leaves);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 0; j < num_tests; j++) {
        while (tree.index < tests[j].tree_index) {
            ret = sparse_tree_next(&tree);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        ret = sparse_tree_get_num_leaves(&tree, tests[j].node, &num_leaves);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_leaves);
        /* We're tracking all leaves, so the count should be the same */
        ret = sparse_tree_get_num_tracked_leaves(&tree, tests[j].node,
                &num_leaves);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(tests[j].count, num_leaves);
        ret = sparse_tree_get_leaf_list(&tree, tests[j].node, &head, &tail);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        u = head;
        k = 0;
        while (1) {
            k++;
            if (u == tail) {
                break;
            }
            u = u->next;
        }
        CU_ASSERT_EQUAL(tests[j].count, k);
    }
    sparse_tree_free(&tree);
    free(tracked_leaves);
}


static void
verify_leaf_sets_for_tree(sparse_tree_t *tree)
{
    int ret, stack_top, j;
    uint32_t u, v, n, num_nodes, num_leaves;
    uint32_t *stack, *leaves;
    leaf_list_node_t *z, *head, *tail;
    tree_sequence_t *ts = tree->tree_sequence;

    n = tree_sequence_get_sample_size(ts);
    num_nodes = tree_sequence_get_num_nodes(ts);
    stack = malloc(n * sizeof(uint32_t));
    leaves = malloc(n * sizeof(uint32_t));
    CU_ASSERT_FATAL(stack != NULL);
    CU_ASSERT_FATAL(leaves != NULL);
    for (u = 0; u < num_nodes; u++) {
        if (tree->num_children[u] == 0 && u >= n) {
            ret = sparse_tree_get_leaf_list(tree, u, &head, &tail);
            CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        } else {
            stack_top = 0;
            num_leaves = 0;
            stack[stack_top] = u;
            while (stack_top >= 0) {
                v = stack[stack_top];
                stack_top--;
                if (v < n) {
                    leaves[num_leaves] = v;
                    num_leaves++;
                }
                for (j = tree->num_children[v] - 1; j >= 0; j--) {
                    stack_top++;
                    stack[stack_top] = tree->children[v][j];
                }
            }
            ret = sparse_tree_get_num_leaves(tree, u, &v);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(num_leaves, v);
            ret = sparse_tree_get_leaf_list(tree, u, &head, &tail);
            CU_ASSERT_EQUAL(ret, 0);
            z = head;
            j = 0;
            while (1) {
                CU_ASSERT_TRUE_FATAL(j < num_leaves);
                CU_ASSERT_EQUAL(z->node, leaves[j]);
                j++;
                if (z == tail) {
                    break;
                }
                z = z->next;
            }
            CU_ASSERT_EQUAL(j, num_leaves);
        }
    }
    free(stack);
    free(leaves);
}

static void
verify_leaf_sets(tree_sequence_t *ts)
{
    int ret;
    sparse_tree_t t;

    ret = sparse_tree_alloc(&t, ts, MSP_LEAF_COUNTS|MSP_LEAF_LISTS);
    CU_ASSERT_EQUAL(ret, 0);

    for (ret = sparse_tree_first(&t); ret == 1; ret = sparse_tree_next(&t)) {
        verify_leaf_sets_for_tree(&t);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (ret = sparse_tree_last(&t); ret == 1; ret = sparse_tree_prev(&t)) {
        verify_leaf_sets_for_tree(&t);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    sparse_tree_free(&t);
}

static sparse_tree_t *
get_tree_list(tree_sequence_t *ts)
{
    int ret;
    sparse_tree_t t, *trees;
    size_t num_trees;

    num_trees = tree_sequence_get_num_trees(ts);
    ret = sparse_tree_alloc(&t, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    trees = malloc(num_trees * sizeof(sparse_tree_t));
    CU_ASSERT_FATAL(trees != NULL);
    for (ret = sparse_tree_first(&t); ret == 1; ret = sparse_tree_next(&t)) {
        CU_ASSERT_FATAL(t.index < num_trees);
        ret = sparse_tree_alloc(&trees[t.index], ts, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = sparse_tree_copy(&trees[t.index], &t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = sparse_tree_equal(&trees[t.index], &t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        /* Make sure the left and right coordinates are also OK */
        CU_ASSERT_DOUBLE_EQUAL(trees[t.index].left, t.left, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(trees[t.index].right, t.right, 1e-6);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    return trees;
}

static void
verify_tree_equals(tree_sequence_t *ts)
{
    int ret;
    sparse_tree_t *trees, t;
    size_t j, k;
    tree_sequence_t *other_ts = get_example_tree_sequence(
            10, 0, 100, 100.0, 1.0, 1.0, 0, NULL);
    int flags[] = {0, MSP_LEAF_LISTS, MSP_LEAF_COUNTS,
        MSP_LEAF_LISTS | MSP_LEAF_COUNTS};

    trees = get_tree_list(ts);
    ret = sparse_tree_alloc(&t, other_ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < tree_sequence_get_num_trees(ts); j++) {
        ret = sparse_tree_equal(&t, &trees[j]);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
        for (k = 0; k < tree_sequence_get_num_trees(ts); k++) {
            ret = sparse_tree_equal(&trees[j], &trees[k]);
            if (j == k) {
                CU_ASSERT_EQUAL_FATAL(ret, 0);
            } else {
                CU_ASSERT_EQUAL_FATAL(ret, 1);
            }
        }
    }
    ret = sparse_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < sizeof(flags) / sizeof(int); j++) {
        ret = sparse_tree_alloc(&t, ts, flags[j]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        for (ret = sparse_tree_first(&t); ret == 1;
                ret = sparse_tree_next(&t)) {
            for (k = 0; k < tree_sequence_get_num_trees(ts); k++) {
                ret = sparse_tree_equal(&t, &trees[k]);
                if (t.index == k) {
                    CU_ASSERT_EQUAL_FATAL(ret, 0);
                } else {
                    CU_ASSERT_EQUAL_FATAL(ret, 1);
                }
            }
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = sparse_tree_free(&t);
        CU_ASSERT_EQUAL(ret, 0);
    }
    for (j = 0; j < tree_sequence_get_num_trees(ts); j++) {
        ret = sparse_tree_free(&trees[j]);
    }
    free(trees);
    tree_sequence_free(other_ts);
    free(other_ts);
}

static void
verify_tree_next_prev(tree_sequence_t *ts)
{
    int ret;
    sparse_tree_t *trees, t;
    size_t j;
    size_t num_trees = tree_sequence_get_num_trees(ts);

    trees = get_tree_list(ts);
    ret = sparse_tree_alloc(&t, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Single forward pass */
    j = 0;
    for (ret = sparse_tree_first(&t); ret == 1; ret = sparse_tree_next(&t)) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);

    /* Single reverse pass */
    j = num_trees;
    for (ret = sparse_tree_last(&t); ret == 1; ret = sparse_tree_prev(&t)) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);

    /* Full forward, then reverse */
    j = 0;
    for (ret = sparse_tree_first(&t); ret == 1; ret = sparse_tree_next(&t)) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);
    j--;
    while ((ret = sparse_tree_prev(&t)) == 1) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);
    CU_ASSERT_EQUAL_FATAL(t.index, 0);
    /* Calling prev should return 0 and have no effect. */
    for (j = 0; j < 10; j++) {
        ret = sparse_tree_prev(&t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, 0);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Full reverse then forward */
    j = num_trees;
    for (ret = sparse_tree_last(&t); ret == 1; ret = sparse_tree_prev(&t)) {
        CU_ASSERT_EQUAL_FATAL(j - 1, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j--;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, 0);
    j++;
    while ((ret = sparse_tree_next(&t)) == 1) {
        CU_ASSERT_EQUAL_FATAL(j, t.index);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(j, num_trees);
    CU_ASSERT_EQUAL_FATAL(t.index, num_trees - 1);
    /* Calling next should return 0 and have no effect. */
    for (j = 0; j < 10; j++) {
        ret = sparse_tree_next(&t);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(t.index, num_trees - 1);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Do a zigzagging traversal */
    ret = sparse_tree_first(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    for (j = 1; j < GSL_MIN(10, num_trees / 2); j++) {
        while (t.index < num_trees - j) {
            ret = sparse_tree_next(&t);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        CU_ASSERT_EQUAL_FATAL(t.index, num_trees - j);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        while (t.index > j) {
            ret = sparse_tree_prev(&t);
            CU_ASSERT_EQUAL_FATAL(ret, 1);
        }
        CU_ASSERT_EQUAL_FATAL(t.index, j);
        ret = sparse_tree_equal(&t, &trees[t.index]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    /* Free the trees. */
    ret = sparse_tree_free(&t);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < tree_sequence_get_num_trees(ts); j++) {
        ret = sparse_tree_free(&trees[j]);
    }
    free(trees);
}

static void
test_leaf_sets(void)
{
    int ret;
    const char * text_records[] = {
        "2 10 4 2,3 0.071 0",
        "0 2  5 1,3 0.090 0",
        "2 10 5 1,4 0.090 0",
        "0 7  6 0,5 0.170 0",
        "7 10 7 0,5 0.202 0",
        "0 2  8 2,6 0.253 0"
    };
    leaf_count_test_t tests[] = {
        {0, 0, 1}, {0, 5, 2}, {0, 6, 3},
        {1, 4, 2}, {1, 5, 3}, {1, 6, 4}};
    size_t num_records = 6;
    uint32_t num_tests = 6;
    tree_sequence_t ts;
    coalescence_record_t *records;

    parse_text_records(num_records, text_records, &records);
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    verify_leaf_counts(&ts, num_tests, tests);
    verify_leaf_sets(&ts);

    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

static void
test_nonbinary_leaf_sets(void)
{
    int ret;
    const char * text_records[] = {
        "0	100	8	0,1,2,3	0.01	0",
        "0	100	9	6,8	    0.068   0",
        "0	17	10	4,5,7	0.2	    0",
        "17	100	10	4,7	    0.2	    0",
        "17	100	11	5,9	    0.279   0",
        "0	17	12	9,10	0.405   0",
        "17	100	12	10,11	0.405   0",
    };
    size_t num_records = sizeof(text_records) / sizeof(char *);
    leaf_count_test_t tests[] = {
        {0, 0, 1}, {0, 8, 4}, {0, 9, 5}, {0, 10, 3}, {0, 12, 8},
        {1, 5, 1}, {1, 8, 4}, {1, 9, 5}, {0, 10, 2}, {0, 11, 1}};
    uint32_t num_tests = 8;
    coalescence_record_t *records;
    tree_sequence_t ts;

    parse_text_records(num_records, text_records, &records);
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    verify_leaf_counts(&ts, num_tests, tests);
    verify_leaf_sets(&ts);

    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

static void
test_tree_sequence_bad_records(void)
{
    int ret = 0;

    const char * text_records[] = {
        "2 10 4 2,3 0.071 0",
        "0 2  5 1,3 0.090 0",
        "2 10 5 1,4 0.090 0",
        "0 7  6 0,5 0.170 0",
        "7 10 7 0,5 0.202 0",
        "0 2  8 2,6 0.253 0"
    };
    size_t num_records = 6;
    tree_sequence_t ts;
    coalescence_record_t *records;

    parse_text_records(num_records, text_records, &records);

    /* Not sorted in time order */
    records[2].time = 0.5;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[2].time = 0.090;

    /* Left value greater than right */
    records[0].left = 10.0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].left = 2.0;

    /* Children equal */
    records[3].children[1] = 0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[3].children[1] = 5;

    /* Children not sorted */
    records[3].children[0] = 5;
    records[3].children[1] = 0;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[3].children[0] = 0;
    records[3].children[1] = 5;

    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    tree_sequence_free(&ts);

    free_local_records(num_records, records);
}

static void
test_single_tree_iter_failure(void)
{
    int ret = 0;
    const char * text_records[] = {
        "0 1 4 0,1 1.0 0",
        "0 1 5 2,3 2.0 0",
        "0 1 6 4,5 3.0 0"
    };
    size_t num_records = 3;
    tree_sequence_t ts;
    sparse_tree_t tree;
    tree_diff_iterator_t diff_iter;
    coalescence_record_t *records;
    node_record_t *records_out, *records_in;
    double length;

    parse_text_records(num_records, text_records, &records);

    /* change the left coordinate of one record so we can't build the tree */
    records[0].left = 0.5;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_alloc(&tree, &ts, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);

    ret = tree_diff_iterator_alloc(&diff_iter, &ts);
    CU_ASSERT_EQUAL(ret, 0);

    ret = tree_diff_iterator_next(&diff_iter, &length, &records_out,
            &records_in);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);

    tree_diff_iterator_free(&diff_iter);
    sparse_tree_free(&tree);
    tree_sequence_free(&ts);
    free_local_records(num_records, records);
}

static void
test_tree_sequence_iter_failure(void)
{
    const char * text_records[] = {
        "2 10 4 2,3 0.071 0",
        "0 2  5 1,3 0.090 0",
        "2 10 5 1,4 0.090 0",
        "0 7  6 0,5 0.170 0",
        "7 10 7 0,5 0.202 0",
        "0 2  8 2,6 0.253 0"
    };
    uint32_t parents[] = {
        6, 5, 8, 5, MSP_NULL_NODE, 6, 8, MSP_NULL_NODE, MSP_NULL_NODE,
        6, 5, 4, 4, 5, 6, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
        7, 5, 4, 4, 5, 7, MSP_NULL_NODE, MSP_NULL_NODE, MSP_NULL_NODE,
    };
    size_t num_records = 6;
    uint32_t num_nodes = 9;
    uint32_t num_trees = 3;
    coalescence_record_t *records;
    tree_sequence_t ts;
    int ret;

    parse_text_records(num_records, text_records, &records);

    /* The first tree is missing a record */
    records[5].left = 1;
    verify_tree_iter_fails(num_records, records, 0, NULL, 0,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    verify_diff_iter_fails(num_records, records, 0,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    records[5].left = 0;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Make a gap between adjacent records */
    records[1].right = 1;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[1].right = 2;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Make a gap in the middle of the sequence */
    records[0].left = 7;
    records[2].left = 7;
    records[3].right = 2;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[0].left = 2;
    records[2].left = 2;
    records[3].right = 7;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Make a gap before the last tree */
    records[4].left = 8;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[4].left = 7;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Add an extra record to the first tree */
    records[4].left = 2;
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_COALESCENCE_RECORDS);
    tree_sequence_free(&ts);
    records[4].left = 7;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Add an extra record to the second tree */
    records[0].left = 0;
    verify_tree_iter_fails(num_records, records, 0, NULL, 0,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    verify_diff_iter_fails(num_records, records, 0,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    records[0].left = 2;
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    /* Remove the last record */
    verify_tree_iter_fails(num_records - 1, records, 0, NULL, 0,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    verify_diff_iter_fails(num_records - 1, records, 0,
            MSP_ERR_BAD_COALESCENCE_RECORDS);
    verify_trees(num_records, records, num_trees, num_nodes, parents, 0, NULL);

    free_local_records(num_records, records);
}

static void
verify_tree_diffs(tree_sequence_t *ts)
{
    int ret;
    tree_diff_iterator_t iter;
    sparse_tree_t tree;
    node_record_t *record, *records_out, *records_in;
    size_t num_nodes = tree_sequence_get_num_nodes(ts);
    size_t j, k, num_in, num_out, num_trees;
    double length, t, x;
    sample_t sample;
    uint32_t u;
    uint32_t *pi = malloc(num_nodes * sizeof(uint32_t));
    double *tau = malloc(num_nodes * sizeof(double));
    int first_tree;

    CU_ASSERT_FATAL(pi != NULL);
    CU_ASSERT_FATAL(tau != NULL);
    for (j = 0; j < num_nodes; j++) {
        pi[j] = MSP_NULL_NODE;
        tau[j] = 0.0;
    }
    ret = tree_diff_iterator_alloc(&iter, ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_alloc(&tree, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_first(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    tree_diff_iterator_print_state(&iter, _devnull);
    for (j = 0; j < tree_sequence_get_sample_size(ts); j++) {
        ret = tree_sequence_get_sample(ts, j, &sample);
        CU_ASSERT_EQUAL(ret, 0);
        tau[j] = sample.time;
    }

    first_tree = 1;
    x = 0.0;
    num_trees = 0;
    while ((ret = tree_diff_iterator_next(
                &iter, &length, &records_out, &records_in)) == 1) {
        tree_diff_iterator_print_state(&iter, _devnull);
        num_trees++;
        record = records_out;
        num_out = 0;
        while (record != NULL) {
            for (k = 0; k < record->num_children; k++) {
                pi[record->children[k]] = MSP_NULL_NODE;
            }
            tau[record->node] = 0;
            num_out += record->num_children - 1;
            record = record->next;
        }
        record = records_in;
        num_in = 0;
        while (record != NULL) {
            for (k = 0; k < record->num_children; k++) {
                pi[record->children[k]] = record->node;
            }
            tau[record->node] = record->time;
            num_in += record->num_children - 1;
            record = record->next;
        }
        if (first_tree) {
            CU_ASSERT_EQUAL(num_in, tree_sequence_get_sample_size(ts) - 1);
        } else {
            CU_ASSERT_EQUAL(num_in, num_out);
        }
        /* Now check against the sparse tree iterator. */
        for (j = 0; j < num_nodes; j++) {
            ret = sparse_tree_get_parent(&tree, j, &u);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(pi[j], u);
            ret = sparse_tree_get_time(&tree, j, &t);
            CU_ASSERT_EQUAL(ret, 0);
            CU_ASSERT_EQUAL(tau[j], t);
        }
        CU_ASSERT_EQUAL(tree.left, x);
        x += length;
        CU_ASSERT_EQUAL(tree.right, x);
        first_tree = 0;
        ret = sparse_tree_next(&tree);
        if (num_trees < tree_sequence_get_num_trees(ts)) {
            CU_ASSERT_EQUAL(ret, 1);
        } else {
            CU_ASSERT_EQUAL(ret, 0);
        }
    }
    if (num_trees != tree_sequence_get_num_trees(ts)) {
        printf("ERROR: %d %d\n", (int) num_trees, (int) tree_sequence_get_num_trees(ts));
    }
    CU_ASSERT_EQUAL(num_trees, tree_sequence_get_num_trees(ts));
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_next(&tree);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_diff_iterator_free(&iter);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = sparse_tree_free(&tree);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    free(pi);
    free(tau);
}

static void
test_tree_sequence_diff_iter(void)
{
    int ret;
    const char * text_records[] = {
        "2 10 4 2,3 0.071 0",
        "0 2  5 1,3 0.090 0",
        "2 10 5 1,4 0.090 0",
        "0 7  6 0,5 0.170 0",
        "7 10 7 0,5 0.202 0",
        "0 2  8 2,6 0.253 0"
    };
    size_t num_records = 6;
    coalescence_record_t *records;
    tree_sequence_t ts;

    parse_text_records(num_records, text_records, &records);
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    verify_tree_diffs(&ts);

    ret = tree_sequence_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
    free_local_records(num_records, records);
}

static void
test_nonbinary_tree_sequence_diff_iter(void)
{
    int ret;
    const char * text_records[] = {
        "0	100	8	0,1,2,3	0.01	0",
        "0	100	9	6,8	    0.068   0",
        "0	17	10	4,5,7	0.2	    0",
        "17	100	10	4,7	    0.2	    0",
        "17	100	11	5,9	    0.279   0",
        "0	17	12	9,10	0.405   0",
        "17	100	12	10,11	0.405   0",
    };
    size_t num_records = sizeof(text_records) / sizeof(char *);
    coalescence_record_t *records;
    tree_sequence_t ts;

    parse_text_records(num_records, text_records, &records);
    ret = tree_sequence_load_records(&ts, num_records, records);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    verify_tree_diffs(&ts);

    ret = tree_sequence_free(&ts);
    CU_ASSERT_EQUAL(ret, 0);
    free_local_records(num_records, records);
}

static void
test_diff_iter_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_tree_diffs(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_tree_iter_from_examples(void)
{

    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_trees_consistent(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_leaf_sets_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_leaf_sets(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_tree_equals_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_tree_equals(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
test_next_prev_from_examples(void)
{

    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_tree_next_prev(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_hapgen(tree_sequence_t *ts)
{
    int ret;
    hapgen_t hapgen;
    char *haplotype;
    size_t sample_size = tree_sequence_get_sample_size(ts);
    size_t num_mutations = tree_sequence_get_num_mutations(ts);
    size_t j;

    ret = hapgen_alloc(&hapgen, ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    hapgen_print_state(&hapgen, _devnull);

    for (j = 0; j < sample_size; j++) {
        ret = hapgen_get_haplotype(&hapgen, j, &haplotype);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(strlen(haplotype), num_mutations);
    }
    for (j = sample_size; j < sample_size + 10; j++) {
        ret = hapgen_get_haplotype(&hapgen, j, &haplotype);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    }
    ret = hapgen_free(&hapgen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
}

static void
test_hapgen_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_hapgen(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_ld(tree_sequence_t *ts)
{
    int ret;
    size_t num_mutations = tree_sequence_get_num_mutations(ts);
    mutation_t *mutations;
    ld_calc_t ld_calc;
    double *r2, *r2_prime, x;
    size_t j, num_r2_values;
    double eps = 1e-6;

    r2 = malloc(num_mutations * sizeof(double));
    r2_prime = malloc(num_mutations * sizeof(double));
    CU_ASSERT_FATAL(r2 != NULL);

    ret = ld_calc_alloc(&ld_calc, ts);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ld_calc_print_state(&ld_calc, _devnull);

    ret = tree_sequence_get_mutations(ts, &mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < num_mutations; j++) {
        ret = ld_calc_get_r2(&ld_calc, j, j, &x);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_DOUBLE_EQUAL_FATAL(x, 1.0, eps);
    }

    if (num_mutations > 0) {
        /* Some checks in the forward direction */
        ret = ld_calc_get_r2_array(&ld_calc, 0, MSP_DIR_FORWARD,
                num_mutations, DBL_MAX, r2, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, num_mutations - 1);
        ld_calc_print_state(&ld_calc, _devnull);

        ret = ld_calc_get_r2_array(&ld_calc, num_mutations - 2, MSP_DIR_FORWARD,
                num_mutations, DBL_MAX, r2_prime, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);
        ld_calc_print_state(&ld_calc, _devnull);

        ret = ld_calc_get_r2_array(&ld_calc, 0, MSP_DIR_FORWARD,
                num_mutations, DBL_MAX, r2_prime, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, num_mutations - 1);
        ld_calc_print_state(&ld_calc, _devnull);

        for (j = 0; j < num_r2_values; j++) {
            CU_ASSERT_EQUAL_FATAL(r2[j], r2_prime[j]);
            ret = ld_calc_get_r2(&ld_calc, 0, j + 1, &x);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_DOUBLE_EQUAL_FATAL(r2[j], x, eps);
        }

        /* Some checks in the reverse direction */
        ret = ld_calc_get_r2_array(&ld_calc, num_mutations - 1,
                MSP_DIR_REVERSE, num_mutations, DBL_MAX,
                r2, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, num_mutations - 1);
        ld_calc_print_state(&ld_calc, _devnull);

        ret = ld_calc_get_r2_array(&ld_calc, 1, MSP_DIR_REVERSE,
                num_mutations, DBL_MAX, r2_prime, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);
        ld_calc_print_state(&ld_calc, _devnull);

        ret = ld_calc_get_r2_array(&ld_calc, num_mutations - 1,
                MSP_DIR_REVERSE, num_mutations, DBL_MAX,
                r2_prime, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, num_mutations - 1);
        ld_calc_print_state(&ld_calc, _devnull);

        for (j = 0; j < num_r2_values; j++) {
            CU_ASSERT_EQUAL_FATAL(r2[j], r2_prime[j]);
            ret = ld_calc_get_r2(&ld_calc, num_mutations - 1,
                    num_mutations - j - 2, &x);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            CU_ASSERT_DOUBLE_EQUAL_FATAL(r2[j], x, eps);
        }

        /* Check some error conditions */
        ret = ld_calc_get_r2_array(&ld_calc, 0, 0, num_mutations, DBL_MAX,
            r2, &num_r2_values);
        CU_ASSERT_EQUAL(ret, MSP_ERR_BAD_PARAM_VALUE);
    }

    if (num_mutations > 3) {
        /* Check for some basic distance calculations */
        j = num_mutations / 2;
        x = mutations[j + 1].position - mutations[j].position;
        ret = ld_calc_get_r2_array(&ld_calc, j, MSP_DIR_FORWARD, num_mutations,
                x, r2, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);

        x = mutations[j].position - mutations[j - 1].position;
        ret = ld_calc_get_r2_array(&ld_calc, j, MSP_DIR_REVERSE, num_mutations,
                x, r2, &num_r2_values);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(num_r2_values, 1);
    }

    /* Check some error conditions */
    for (j = num_mutations; j < num_mutations + 2; j++) {
        ret = ld_calc_get_r2_array(&ld_calc, j, MSP_DIR_FORWARD,
                num_mutations, DBL_MAX, r2, &num_r2_values);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = ld_calc_get_r2(&ld_calc, j, 0, r2);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
        ret = ld_calc_get_r2(&ld_calc, 0, j, r2);
        CU_ASSERT_EQUAL(ret, MSP_ERR_OUT_OF_BOUNDS);
    }

    ld_calc_free(&ld_calc);
    free(r2);
    free(r2_prime);
}

static void
test_ld_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_ld(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_vargen(tree_sequence_t *ts)
{
    int ret;
    vargen_t vargen;
    mutation_t *mut;
    size_t sample_size = tree_sequence_get_sample_size(ts);
    size_t num_mutations = tree_sequence_get_num_mutations(ts);
    char *genotypes = malloc((sample_size + 1) * sizeof(char));
    size_t j, k;

    CU_ASSERT_FATAL(genotypes != NULL);
    ret = vargen_alloc(&vargen, ts, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    vargen_print_state(&vargen, _devnull);
    j = 0;
    while ((ret = vargen_next(&vargen, &mut, genotypes)) == 1) {
        CU_ASSERT_EQUAL(mut->index, j);
        for (k = 0; k < sample_size; k++) {
            CU_ASSERT(genotypes[k] == 0 || genotypes[k] == 1);
        }
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(j, num_mutations);
    CU_ASSERT_EQUAL_FATAL(vargen_next(&vargen, &mut, genotypes), 0);
    ret = vargen_free(&vargen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = vargen_alloc(&vargen, ts, MSP_GENOTYPES_AS_CHAR);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    vargen_print_state(&vargen, _devnull);
    j = 0;
    genotypes[sample_size] = '\0';
    while ((ret = vargen_next(&vargen, &mut, genotypes)) == 1) {
        CU_ASSERT_EQUAL(mut->index, j);
        for (k = 0; k < sample_size; k++) {
            CU_ASSERT(genotypes[k] == '0' || genotypes[k] == '1');
        }
        CU_ASSERT_EQUAL_FATAL(genotypes[sample_size], '\0');
        j++;
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(j, num_mutations);
    CU_ASSERT_EQUAL_FATAL(vargen_next(&vargen, &mut, genotypes), 0);
    ret = vargen_free(&vargen);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    free(genotypes);
}

static void
test_vargen_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(1);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_vargen(examples[j]);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_newick(tree_sequence_t *ts, int should_fail)
{
    newick_converter_t nc;
    double length;
    char *tree;
    int ret;

    ret = newick_converter_alloc(&nc, ts, 1, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    while ((ret = newick_converter_next(&nc, &length, &tree)) == 1) {
        CU_ASSERT(length > 0);
        newick_converter_print_state(&nc, _devnull);
        CU_ASSERT_FATAL(tree != NULL);
        CU_ASSERT(strlen(tree) > 0);
    }
    if (should_fail) {
        CU_ASSERT_EQUAL(ret, MSP_ERR_NONBINARY_NEWICK);
    } else {
        CU_ASSERT_EQUAL(ret, 0);
    }
    newick_converter_free(&nc);
}

static void
test_newick_from_examples(void)
{
    tree_sequence_t **examples = get_example_tree_sequences(0);
    uint32_t j;

    CU_ASSERT_FATAL(examples != NULL);
    for (j = 0; examples[j] != NULL; j++) {
        verify_newick(examples[j], 0);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);

    examples = get_example_nonbinary_tree_sequences();
    for (j = 0; examples[j] != NULL; j++) {
        verify_newick(examples[j], 1);
        tree_sequence_free(examples[j]);
        free(examples[j]);
    }
    free(examples);
}

static void
verify_tree_sequences_equal(tree_sequence_t *ts1, tree_sequence_t *ts2,
        int check_provenance_strings)
{
    int ret, err1, err2;
    size_t j, nps1, nps2;
    sample_t sample1, sample2;
    coalescence_record_t *r1, *r2;
    char **ps1, **ps2;
    size_t num_mutations = tree_sequence_get_num_mutations(ts1);
    mutation_t *mutations_1, *mutations_2;
    sparse_tree_t t1, t2;

    CU_ASSERT_EQUAL(
        tree_sequence_get_sample_size(ts1),
        tree_sequence_get_sample_size(ts2))
    CU_ASSERT_EQUAL(
        tree_sequence_get_sequence_length(ts1),
        tree_sequence_get_sequence_length(ts2))
    CU_ASSERT_EQUAL(
        tree_sequence_get_num_coalescence_records(ts1),
        tree_sequence_get_num_coalescence_records(ts2));
    CU_ASSERT_EQUAL_FATAL(
        tree_sequence_get_num_mutations(ts1),
        tree_sequence_get_num_mutations(ts2));
    CU_ASSERT_EQUAL(
        tree_sequence_get_num_nodes(ts1),
        tree_sequence_get_num_nodes(ts2));
    CU_ASSERT_EQUAL(
        tree_sequence_get_num_trees(ts1),
        tree_sequence_get_num_trees(ts2));

    for (j = 0; j < tree_sequence_get_num_coalescence_records(ts1); j++) {
        ret = tree_sequence_get_record(ts1, j, &r1, MSP_ORDER_TIME);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_get_record(ts2, j, &r2, MSP_ORDER_TIME);
        CU_ASSERT_EQUAL(ret, 0);
        verify_coalescence_records_equal(r1, r2, 1.0);
        ret = tree_sequence_get_record(ts1, j, &r1, MSP_ORDER_LEFT);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_get_record(ts2, j, &r2, MSP_ORDER_LEFT);
        CU_ASSERT_EQUAL(ret, 0);
        verify_coalescence_records_equal(r1, r2, 1.0);
        ret = tree_sequence_get_record(ts1, j, &r1, MSP_ORDER_RIGHT);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_get_record(ts2, j, &r2, MSP_ORDER_RIGHT);
        CU_ASSERT_EQUAL(ret, 0);
        verify_coalescence_records_equal(r1, r2, 1.0);
    }

    ret = tree_sequence_get_mutations(ts1, &mutations_1);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tree_sequence_get_mutations(ts2, &mutations_2);
    CU_ASSERT_EQUAL(ret, 0);
    for (j = 0; j < num_mutations; j++) {
        CU_ASSERT_EQUAL(mutations_1[j].index, j);
        CU_ASSERT_EQUAL(mutations_1[j].index, mutations_2[j].index);
        CU_ASSERT_EQUAL(mutations_1[j].position, mutations_2[j].position);
        CU_ASSERT_EQUAL(mutations_1[j].node, mutations_2[j].node);
    }

    for (j = 0; j < tree_sequence_get_sample_size(ts1); j++) {
        ret = tree_sequence_get_sample(ts1, (uint32_t) j, &sample1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_get_sample(ts2, (uint32_t) j, &sample2);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL(sample1.population_id, sample2.population_id);
        CU_ASSERT_EQUAL(sample1.time, sample2.time);
    }

    if (check_provenance_strings) {
        ret = tree_sequence_get_provenance_strings(ts1, &nps1, &ps1);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_get_provenance_strings(ts2, &nps2, &ps2);
        CU_ASSERT_EQUAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(nps1, nps2);
        for (j = 0; j < nps1; j++) {
            CU_ASSERT_FATAL(ps1[j] != NULL);
            CU_ASSERT_FATAL(ps2[j] != NULL);
            CU_ASSERT_STRING_EQUAL(ps1[j], ps2[j]);
        }
    }
    ret = sparse_tree_alloc(&t1, ts1, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_alloc(&t2, ts2, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = sparse_tree_first(&t1);
    CU_ASSERT_EQUAL(ret, 1);
    ret = sparse_tree_first(&t2);
    CU_ASSERT_EQUAL(ret, 1);
    while (1) {
        err1 = sparse_tree_next(&t1);
        err2 = sparse_tree_next(&t2);
        CU_ASSERT_EQUAL_FATAL(err1, err2);
        if (err1 != 1) {
            break;
        }
    }
    sparse_tree_free(&t1);
    sparse_tree_free(&t2);
}

static void
test_save_hdf5(void)
{
    int ret;
    size_t j, k;
    tree_sequence_t **examples = get_example_tree_sequences(1);
    tree_sequence_t ts2;
    tree_sequence_t *ts1;
    int dump_flags[] = {0, MSP_ZLIB_COMPRESSION};

    CU_ASSERT_FATAL(examples != NULL);

    for (j = 0; examples[j] != NULL; j++) {
        ts1 = examples[j];
        for (k = 0; k < sizeof(dump_flags) / sizeof(int); k++) {
            ret = tree_sequence_dump(ts1, _tmp_file_name, dump_flags[k]);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            ret = tree_sequence_load(&ts2, _tmp_file_name, 0);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
            verify_tree_sequences_equal(ts1, &ts2, 1);
            tree_sequence_print_state(&ts2, _devnull);
            verify_hapgen(&ts2);
            verify_vargen(&ts2);
            tree_sequence_free(&ts2);
        }
        tree_sequence_free(ts1);
        free(ts1);
    }
    free(examples);
}

static void
test_save_records_hdf5(void)
{
    int ret;
    size_t j, k, num_records, num_mutations;
    uint32_t sample_size;
    coalescence_record_t *r, *records;
    sample_t *samples;
    mutation_t *mutations;
    tree_sequence_t *ts1, ts2, ts3, **examples;

    examples = get_example_tree_sequences(1);
    for (k = 0; examples[k] != NULL; k++) {
        ts1 = examples[k];
        CU_ASSERT_FATAL(ts1 != NULL);
        sample_size = tree_sequence_get_sample_size(ts1);
        num_records = tree_sequence_get_num_coalescence_records(ts1);
        records = malloc(num_records * sizeof(coalescence_record_t));
        CU_ASSERT_FATAL(records != NULL);
        for (j = 0; j < num_records; j++) {
            ret = tree_sequence_get_record(ts1, j, &r, MSP_ORDER_TIME);
            CU_ASSERT_EQUAL(ret, 0);
            copy_record(&records[j], r);
        }
        samples = malloc(sample_size * sizeof(sample_t));
        CU_ASSERT_FATAL(samples != NULL);
        for (j = 0; j < sample_size; j++) {
            ret = tree_sequence_get_sample(ts1, (uint32_t) j, &samples[j]);
            CU_ASSERT_EQUAL(ret, 0);
        }
        num_mutations = tree_sequence_get_num_mutations(ts1);
        ret = tree_sequence_get_mutations(ts1, &mutations);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_load_records(&ts2, num_records, records);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_set_samples(&ts2, sample_size, samples);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_set_mutations(&ts2, num_mutations, mutations);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_dump(&ts2, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tree_sequence_load(&ts3, _tmp_file_name, 0);
        CU_ASSERT_EQUAL(ret, 0);
        verify_tree_sequences_equal(ts1, &ts3, 0);
        tree_sequence_print_state(&ts2, _devnull);

        tree_sequence_free(&ts2);
        tree_sequence_free(&ts3);
        tree_sequence_free(ts1);
        free(ts1);
        free(samples);
        free_local_records(num_records, records);
    }
    free(examples);
}

static void
test_records_equivalent(void)
{
    int ret;
    tree_sequence_t *ts1 = get_example_tree_sequence(10, 0, 100, 10.0, 1.0,
            1.0, 0, NULL);
    tree_sequence_t ts2;
    coalescence_record_t *records, *r1, *r2;
    size_t j, num_records;

    CU_ASSERT_FATAL(ts1 != NULL);
    num_records = tree_sequence_get_num_coalescence_records(ts1);
    records = malloc(num_records * sizeof(coalescence_record_t));
    CU_ASSERT_FATAL(records != NULL);
    for (j = 0; j < num_records; j++) {
        ret = tree_sequence_get_record(ts1, j, &r1, MSP_ORDER_TIME);
        CU_ASSERT_EQUAL(ret, 0);
        copy_record(&records[j], r1);
    }
    ret = tree_sequence_load_records(&ts2, num_records, records);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
        tree_sequence_get_sample_size(ts1),
        tree_sequence_get_sample_size(&ts2));
    CU_ASSERT_EQUAL(
        tree_sequence_get_sequence_length(ts1),
        tree_sequence_get_sequence_length(&ts2));
    for (j = 0; j < num_records; j++) {
        ret = tree_sequence_get_record(ts1, j, &r1, MSP_ORDER_TIME);
        CU_ASSERT_EQUAL(ret, 0);
        ret = tree_sequence_get_record(&ts2, j, &r2, MSP_ORDER_TIME);
        CU_ASSERT_EQUAL(ret, 0);
        verify_coalescence_records_equal(r1, r2, 1.0);
    }
    tree_sequence_free(&ts2);
    tree_sequence_free(ts1);
    free(ts1);
    free_local_records(num_records, records);
}

static void
test_strerror(void)
{
    int ret;
    int j;
    const char *msg;
    int max_error_code = 1024; /* totally arbitrary */
    FILE *f;
    tree_sequence_t ts;

    for (j = 0; j < max_error_code; j++) {
        msg = msp_strerror(-j);
        CU_ASSERT_FATAL(msg != NULL);
        if (-j == MSP_ERR_HDF5) {
            /* There is no HDF5 error, so... */
            CU_ASSERT_EQUAL(strlen(msg), 0);
        } else {
            CU_ASSERT(strlen(msg) > 0);
        }
    }
    /* Provoke an HDF5 error */
    ret = tree_sequence_load(&ts, "/file/does/not/exist", 0);
    CU_ASSERT_EQUAL(ret, MSP_ERR_HDF5);
    msg = msp_strerror(ret);
    CU_ASSERT_FATAL(msg != NULL);
    CU_ASSERT(strlen(msg) > 0);
    /* Provoke an IO error */
    f = fopen("/file/does/not/exist", "r");
    CU_ASSERT_EQUAL(f, NULL);
    msg = msp_strerror(MSP_ERR_IO);
    CU_ASSERT_FATAL(msg != NULL);
    CU_ASSERT_STRING_EQUAL(msg, strerror(errno));
}

static int
msprime_suite_init(void)
{
    int fd;
    static char template[] = "/tmp/msp_c_test_XXXXXX";

    _tmp_file_name = NULL;
    _devnull = NULL;

    _tmp_file_name = malloc(sizeof(template));
    if (_tmp_file_name == NULL) {
        return CUE_NOMEMORY;
    }
    strcpy(_tmp_file_name, template);
    fd = mkstemp(_tmp_file_name);
    if (fd == -1) {
        return CUE_SINIT_FAILED;
    }
    close(fd);
    _devnull = fopen("/dev/null", "w");
    if (_devnull == NULL) {
        return CUE_SINIT_FAILED;
    }
    /* Silence HDF5 errors */
    if (H5Eset_auto(H5E_DEFAULT, NULL, NULL) < 0) {
        return CUE_SINIT_FAILED;
    }
    return CUE_SUCCESS;
}

static int
msprime_suite_cleanup(void)
{
    if (_tmp_file_name != NULL) {
        unlink(_tmp_file_name);
        free(_tmp_file_name);
    }
    if (_devnull != NULL) {
        fclose(_devnull);
    }
    return CUE_SUCCESS;
}

static void
handle_cunit_error()
{
    fprintf(stderr, "CUnit error occured: %d: %s\n",
            CU_get_error(), CU_get_error_msg());
    exit(EXIT_FAILURE);
}

int
main(void)
{
    CU_TestInfo tests[] = {
        {"Fenwick tree", test_fenwick},
        {"VCF", test_vcf},
        {"Simple recombination map", test_simple_recomb_map},
        {"Recombination map errors", test_recomb_map_errors},
        {"Recombination map examples", test_recomb_map_examples},
        {"Simplest records", test_simplest_records},
        {"Simplest nonbinary records", test_simplest_nonbinary_records},
        {"Simplest bad records", test_simplest_bad_records},
        {"Single tree good records", test_single_tree_good_records},
        {"Single nonbinary tree good records",
            test_single_nonbinary_tree_good_records},
        {"Single tree bad records", test_single_tree_bad_records},
        {"Single tree good mutations", test_single_tree_good_mutations},
        {"Single tree bad mutations", test_single_tree_bad_mutations},
        {"Single tree iterator", test_single_tree_iter},
        {"Single nonbinary tree iterator", test_single_nonbinary_tree_iter},
        {"Single tree iterator times", test_single_tree_iter_times},
        {"Single tree hapgen", test_single_tree_hapgen},
        {"Single tree mutgen", test_single_tree_mutgen},
        {"Sparse tree errors", test_sparse_tree_errors},
        {"Tree sequence iterator", test_tree_sequence_iter},
        {"Leaf sets", test_leaf_sets},
        {"Nonbinary leaf sets", test_nonbinary_leaf_sets},
        {"Tree nonbinary sequence iterator", test_nonbinary_tree_sequence_iter},
        {"Left-to-right tree sequence iterator",
            test_left_to_right_tree_sequence_iter},
        {"Tree sequence bad records", test_tree_sequence_bad_records},
        {"Set mutations", test_tree_sequence_set_mutations},
        {"Single tree iterator failure", test_single_tree_iter_failure},
        {"Tree sequence iterator failure", test_tree_sequence_iter_failure},
        {"Tree sequence diff iter", test_tree_sequence_diff_iter},
        {"Nonbinary Tree sequence diff iter",
            test_nonbinary_tree_sequence_diff_iter},
        {"diff iter from examples", test_diff_iter_from_examples},
        {"tree iter from examples", test_tree_iter_from_examples},
        {"tree equals from examples", test_tree_equals_from_examples},
        {"tree next and prev from examples", test_next_prev_from_examples},
        {"leaf sets from examples", test_leaf_sets_from_examples},
        {"Test hapgen from examples", test_hapgen_from_examples},
        {"Test vargen from examples", test_vargen_from_examples},
        {"Test newick from examples", test_newick_from_examples},
        {"Test ld from examples", test_ld_from_examples},
        {"Test records equivalent after import", test_records_equivalent},
        {"Test saving to HDF5", test_save_hdf5},
        {"Test saving records to HDF5", test_save_records_hdf5},
        {"Historical samples two populations",
            test_single_locus_two_populations},
        {"Historical samples", test_single_locus_historical_sample},
        {"Simulator getters/setters", test_simulator_getters_setters},
        {"Demographic events", test_simulator_demographic_events},
        {"Single locus simulation", test_single_locus_simulation},
        {"Simulation memory limit", test_simulation_memory_limit},
        {"Multi locus simulation", test_multi_locus_simulation},
        {"Bottleneck simulation", test_bottleneck_simulation},
        {"Large bottleneck simulation", test_large_bottleneck_simulation},
        {"Test error messages", test_strerror},
        CU_TEST_INFO_NULL,
    };

    /* We use initialisers here as the struct definitions change between
     * versions of CUnit */
    CU_SuiteInfo suites[] = {
        {
            .pName = "msprime",
            .pInitFunc = msprime_suite_init,
            .pCleanupFunc = msprime_suite_cleanup,
            .pTests = tests
        },
        CU_SUITE_INFO_NULL,
    };
    /* initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry()) {
        handle_cunit_error();
    }
    if (CUE_SUCCESS != CU_register_suites(suites)) {
        handle_cunit_error();
    }
    /* Run all tests using the CUnit Basic interface */
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return EXIT_SUCCESS;
}
