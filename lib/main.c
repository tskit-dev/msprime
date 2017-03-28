/*
** Copyright (C) 2015-2017 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

#include <regex.h>
#include <libconfig.h>
#include <gsl/gsl_math.h>
#include "argtable3.h"

#include "msprime.h"
#include "err.h"

/* This file defines a crude CLI for msprime. It is intended for development
 * use only.
 */

typedef struct {
    int alphabet;
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
        ret_samples[j].population_id = (population_id_t) config_setting_get_int(t);
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
read_model_config(msp_t *msp, config_t *config)
{
    int ret = 0;
    config_setting_t *setting = config_lookup(config, "model");
    config_setting_t *s;
    const char *name;
    double psi, alpha, truncation_point;

    if (setting == NULL) {
        fatal_error("model is a required parameter");
    }
    if (config_setting_is_group(setting) == CONFIG_FALSE) {
        fatal_error("model not a group");
    }
    s = config_setting_get_member(setting, "name");
    if (s == NULL) {
        fatal_error("model name not specified");
    }
    name = config_setting_get_string(s);
    if (strcmp(name, "hudson") == 0) {
        ret = msp_set_simulation_model_non_parametric(msp, MSP_MODEL_HUDSON);
    } else if (strcmp(name, "smc") == 0) {
        ret = msp_set_simulation_model_non_parametric(msp, MSP_MODEL_SMC);
    } else if (strcmp(name, "smc_prime") == 0) {
        ret = msp_set_simulation_model_non_parametric(msp, MSP_MODEL_SMC_PRIME);
    } else if (strcmp(name, "dirac") == 0) {
        s = config_setting_get_member(setting, "psi");
        if (s == NULL) {
            fatal_error("dirac model psi not specified");
        }
        psi = config_setting_get_float(s);
        ret = msp_set_simulation_model_dirac(msp, psi);
    } else if (strcmp(name, "beta") == 0) {
        s = config_setting_get_member(setting, "alpha");
        if (s == NULL) {
            fatal_error("beta model alpha not specified");
        }
        alpha = config_setting_get_float(s);
        s = config_setting_get_member(setting, "truncation_point");
        if (s == NULL) {
            fatal_error("beta model truncation_point not specified");
        }
        truncation_point = config_setting_get_float(s);
        ret = msp_set_simulation_model_beta(msp, alpha, truncation_point);
    } else {
        fatal_error("Unknown simulation model '%s'", name);
    }
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
           intensity, strength;
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
        } else if (strcmp(type, "simple_bottleneck") == 0) {
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
            ret = msp_add_simple_bottleneck(msp, time, population_id, intensity);
        } else if (strcmp(type, "instantaneous_bottleneck") == 0) {
            t = config_setting_get_member(s, "strength");
            if (t == NULL) {
                fatal_error("strength not specified");
            }
            strength = config_setting_get_float(t);
            t = config_setting_get_member(s, "population_id");
            if (t == NULL) {
                fatal_error("population_id not specified");
            }
            population_id = config_setting_get_int(t);
            ret = msp_add_instantaneous_bottleneck(msp, time, population_id,
                    strength);
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
    config_setting_t *t;

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
    if (config_lookup_float(config,
            "mutation_rate", &mutation_params->mutation_rate)
            == CONFIG_FALSE) {
        fatal_error("mutation_rate is a required parameter");
    }
    if (config_lookup_int(config, "mutation_alphabet", &int_tmp)
            == CONFIG_FALSE) {
        fatal_error("mutation_alphabet is a required parameter.");
    }
    mutation_params->alphabet = int_tmp;
    if (config_lookup_int(config, "avl_node_block_size", &int_tmp)
            == CONFIG_FALSE) {
        fatal_error("avl_node_block_size is a required parameter");
    }
    ret = msp_set_avl_node_block_size(msp, (size_t) int_tmp);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    if (config_lookup_int(config, "segment_block_size", &int_tmp)
            == CONFIG_FALSE) {
        fatal_error("segment_block_size is a required parameter");
    }
    ret = msp_set_segment_block_size(msp, (size_t) int_tmp);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    if (config_lookup_int(config, "node_mapping_block_size", &int_tmp)
            == CONFIG_FALSE) {
        fatal_error("node_mapping_block_size is a required parameter");
    }
    ret = msp_set_node_mapping_block_size(msp, (size_t) int_tmp);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    if (config_lookup_int(config, "coalescence_record_block_size", &int_tmp)
            == CONFIG_FALSE) {
        fatal_error("coalescence_record_block_size is a required parameter");
    }
    ret = msp_set_coalescence_record_block_size(msp, (size_t) int_tmp);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    if (config_lookup_int(config, "max_memory", &int_tmp)
            == CONFIG_FALSE) {
        fatal_error("max_memory is a required parameter");
    }
    ret = msp_set_max_memory(msp, (size_t) int_tmp * 1024 * 1024);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    if (config_lookup_int(config, "store_migrations", &int_tmp) == CONFIG_FALSE) {
        fatal_error("store_migrations is a required parameter");
    }
    ret = msp_set_store_migrations(msp, (bool) int_tmp);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    t = config_lookup(config, "model");
    if (t == NULL) {
        fatal_error("model not specified");
    }
    ret = read_model_config(msp, config);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
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
    site_t *site;
    char *genotypes = malloc(tree_sequence_get_sample_size(ts) * sizeof(char));

    if (genotypes == NULL) {
        fatal_error("no memory");
    }
    printf("variants (%d) \n", (int) ts->sites.num_records);
    ret = vargen_alloc(&vg, ts, 0);
    if (ret != 0) {
        fatal_library_error(ret, "vargen_alloc");
    }
    j = 0;
    while ((ret = vargen_next(&vg, &site, genotypes)) == 1) {
        assert(site->mutations_length == 1);
        printf("%d\t%f\t%s\t%s\t", j, site->position, site->ancestral_state,
                site->mutations[0].derived_state);
        for (k = 0; k < tree_sequence_get_sample_size(ts); k++) {
            printf("%d", genotypes[k]);
        }
        printf("\n");
        j++;
    }
    if (ret != 0) {
        fatal_library_error(ret, "vargen_next");
    }
    if (j != ts->sites.num_records) {
        printf("ERROR!! missing variants %d %d\n", j, (int) ts->sites.num_records);
    }

    while ((ret = vargen_next(&vg, &site, genotypes)) == 1) {
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
        ret = hapgen_get_haplotype(&hg, (node_id_t) j, &haplotype);
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
    size_t num_sites = tree_sequence_get_num_sites(ts);
    site_t sA, sB;
    double *r2 = malloc(num_sites * sizeof(double));
    size_t j, k, num_r2_values;
    ld_calc_t ld_calc;

    if (r2 == NULL) {
        fatal_error("no memory");
    }
    ret = ld_calc_alloc(&ld_calc, ts);
    printf("alloc: ret = %d\n", ret);
    if (ret != 0) {
        fatal_library_error(ret, "ld_calc_alloc");
    }
    ld_calc_print_state(&ld_calc, stdout);
    for (j = 0; j < num_sites; j++) {
        ret = ld_calc_get_r2_array(&ld_calc, j, MSP_DIR_FORWARD, num_sites,
                DBL_MAX, r2, &num_r2_values);
        if (ret != 0) {
            fatal_library_error(ret, "ld_calc_get_r2_array");
        }
        for (k = 0; k < num_r2_values; k++) {
            ret = tree_sequence_get_site(ts, (site_id_t) j, &sA);
            if (ret != 0) {
                fatal_library_error(ret, "get_site");
            }
            ret = tree_sequence_get_site(ts, (site_id_t) (j + k + 1), &sB);
            if (ret != 0) {
                fatal_library_error(ret, "get_site");
            }
            printf("%d\t%f\t%d\t%f\t%.3f\n",
                (int) sA.id, sA.position, (int) sB.id, sB.position, r2[k]);
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
    size_t sample_size = tree_sequence_get_sample_size(ts) / 2;
    node_id_t *sample = malloc(sample_size * sizeof(node_id_t));
    double pi;

    if (sample == NULL) {
        fatal_error("no memory");
    }
    for (j = 0; j < sample_size; j++) {
        sample[j] = (node_id_t) j;
    }
    ret = tree_sequence_get_pairwise_diversity(ts, sample, sample_size, &pi);
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
print_tree_sequence(tree_sequence_t *ts, int verbose)
{
    int ret = 0;
    sparse_tree_t tree;

    tree_sequence_print_state(ts, stdout);
    if (verbose > 0) {
        printf("========================\n");
        printf("trees\n");
        printf("========================\n");
        ret = sparse_tree_alloc(&tree, ts, MSP_LEAF_COUNTS);
        if (ret != 0) {
            fatal_error("ERROR: %d: %s\n", ret, msp_strerror(ret));
        }
        for (ret = sparse_tree_first(&tree); ret == 1; ret = sparse_tree_next(&tree)) {
            printf("-------------------------\n");
            printf("New tree: %d: %f (%d)\n", (int) tree.index,
                    tree.right - tree.left, (int) tree.num_nodes);
            printf("-------------------------\n");
            sparse_tree_print_state(&tree, stdout);
        }
        if (ret < 0) {
            fatal_error("ERROR: %d: %s\n", ret, msp_strerror(ret));
        }
        sparse_tree_free(&tree);
    }
}

static void
run_simulate(const char *conf_file, const char *output_file, int verbose, int num_replicates)
{
    int ret = -1;
    int j;
    mutation_params_t mutation_params;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    msp_t *msp = calloc(1, sizeof(msp_t));
    tree_sequence_t *tree_seq = calloc(1, sizeof(tree_sequence_t));
    recomb_map_t *recomb_map = calloc(1, sizeof(recomb_map_t));
    mutgen_t *mutgen = calloc(1, sizeof(mutgen_t));
    const char *provenance[] = {"main.simulate"};
    node_table_t *nodes = malloc(sizeof(node_table_t));
    edgeset_table_t *edgesets = malloc(sizeof(edgeset_table_t));
    site_table_t *sites = malloc(sizeof(site_table_t));
    mutation_table_t *mutations = malloc(sizeof(mutation_table_t));
    migration_table_t *migrations = malloc(sizeof(migration_table_t));


    if (rng == NULL || msp == NULL || tree_seq == NULL || recomb_map == NULL
            || mutgen == NULL || nodes == NULL || edgesets == NULL
            || sites == NULL || mutations == NULL || migrations == NULL) {
        goto out;
    }
    ret = get_configuration(rng, msp, &mutation_params, recomb_map, conf_file);
    if (ret != 0) {
        goto out;
    }
    ret = edgeset_table_alloc(edgesets, 10, 10);
    if (ret != 0) {
        goto out;
    }
    ret = node_table_alloc(nodes, 10, 10);
    if (ret != 0) {
        goto out;
    }
    ret = site_table_alloc(sites, 1, 1);
    if (ret != 0) {
        goto out;
    }
    ret = mutation_table_alloc(mutations, 10, 10);
    if (ret != 0) {
        goto out;
    }
    ret = migration_table_alloc(migrations, 10);
    if (ret != 0) {
        goto out;
    }
    ret = mutgen_alloc(mutgen, mutation_params.mutation_rate, rng,
            mutation_params.alphabet);
    if (ret != 0) {
        goto out;
    }
    ret = msp_initialise(msp);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_initialise(tree_seq);
    if (ret != 0) {
        goto out;
    }
    for (j = 0; j < num_replicates; j++) {
        if (verbose >= 1) {
            printf("=====================\n");
            printf("replicate %d\n", j);
            printf("=====================\n");
        }
        ret = msp_reset(msp);
        if (ret != 0) {
            goto out;
        }
        ret = msp_run(msp, DBL_MAX, UINT32_MAX);
        if (ret < 0) {
            goto out;
        }
        msp_verify(msp);
        if (verbose >= 1) {
            ret = msp_print_state(msp, stdout);
        }
        if (ret != 0) {
            goto out;
        }
        /* Create the tree_sequence from the state of the simulator.
         * We want to use coalescent time here, so use an Ne of 1/4
         * to cancel scaling factor. */
        ret = msp_populate_tables(msp, 0.25, recomb_map, nodes, edgesets, migrations);
        if (ret != 0) {
            goto out;
        }
        ret = mutgen_generate_tables_tmp(mutgen, nodes, edgesets);
        if (ret != 0) {
            goto out;
        }
        ret = mutgen_populate_tables(mutgen, sites, mutations);
        if (ret != 0) {
            goto out;
        }
        ret = tree_sequence_load_tables_tmp(tree_seq, nodes, edgesets, migrations,
                sites, mutations, 1, (char **) &provenance);
        if (ret != 0) {
            goto out;
        }
        if (output_file != NULL) {
            ret = tree_sequence_dump(tree_seq, output_file, 0);
            if (ret != 0) {
                goto out;
            }
        }
        if (verbose >= 1) {
            node_table_print_state(nodes, stdout);
            edgeset_table_print_state(edgesets, stdout);
            site_table_print_state(sites, stdout);
            mutation_table_print_state(mutations, stdout);
            migration_table_print_state(migrations, stdout);
            printf("-----------------\n");
            mutgen_print_state(mutgen, stdout);
            printf("-----------------\n");
            tree_sequence_print_state(tree_seq, stdout);
        }
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
    if (edgesets != NULL) {
        edgeset_table_free(edgesets);
        free(edgesets);
    }
    if (nodes != NULL) {
        node_table_free(nodes);
        free(nodes);
    }
    if (mutations != NULL) {
        mutation_table_free(mutations);
        free(mutations);
    }
    if (sites != NULL) {
        site_table_free(sites);
        free(sites);
    }
    if (migrations != NULL) {
        migration_table_free(migrations);
        free(migrations);
    }
    if (ret != 0) {
        printf("error occured:%d:%s\n", ret, msp_strerror(ret));
    }
}

static void
load_tree_sequence(tree_sequence_t *ts, const char *filename)
{
    int ret = tree_sequence_initialise(ts);
    if (ret != 0) {
        fatal_library_error(ret, "Init error");
    }
    ret = tree_sequence_load(ts, filename, 0);
    if (ret != 0) {
        fatal_library_error(ret, "Load error");
    }
}

static void
run_ld(const char *filename, int verbose)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_ld_matrix(&ts);
    tree_sequence_free(&ts);
}

static void
run_haplotypes(const char *filename, int verbose)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_haplotypes(&ts);
    tree_sequence_free(&ts);
}

static void
run_variants(const char *filename, int verbose)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_variants(&ts);
    tree_sequence_free(&ts);
}

static void
run_vcf(const char *filename, int verbose)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_vcf(&ts, 1);
    tree_sequence_free(&ts);
}

static void
run_print(const char *filename, int verbose)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_tree_sequence(&ts, verbose);
    tree_sequence_free(&ts);
}

static void
run_newick(const char *filename, int verbose)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_newick_trees(&ts);
    tree_sequence_free(&ts);
}

static void
run_stats(const char *filename, int verbose)
{
    tree_sequence_t ts;

    load_tree_sequence(&ts, filename);
    print_stats(&ts);
    tree_sequence_free(&ts);
}

static void
run_simplify(const char *input_filename, const char *output_filename, int verbose)
{
    tree_sequence_t ts, subset;
    size_t j, num_samples;
    node_id_t *samples;
    int flags = 0;
    int ret;

    load_tree_sequence(&ts, input_filename);
    num_samples = tree_sequence_get_sample_size(&ts);
    samples = malloc(num_samples * sizeof(node_id_t));
    if (samples == NULL) {
        fatal_error("out of memory");
    }
    for (j = 0; j < num_samples; j++) {
        samples[j] = (node_id_t) j;
    }
    ret = tree_sequence_initialise(&subset);
    if (ret != 0) {
        fatal_library_error(ret, "init error");
    }
    ret = tree_sequence_simplify(&ts, samples, num_samples, flags, &subset);
    if (ret != 0) {
        fatal_library_error(ret, "Subset error");
    }
    ret = tree_sequence_dump(&subset, output_filename, 0);
    if (ret != 0) {
        fatal_library_error(ret, "Write error");
    }
    /* tree_sequence_print_state(&subset, stdout); */
    tree_sequence_free(&ts);
    tree_sequence_free(&subset);
    free(samples);
}

int
main(int argc, char** argv)
{
    /* SYNTAX 1: simulate [-v] <config-file> -o <output-file> */
    struct arg_rex *cmd1 = arg_rex1(NULL, NULL, "simulate", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose1 = arg_lit0("v", "verbose", NULL);
    struct arg_int *replicates1 = arg_int0("r", "replicates", "<num-replicates>",
            "number of replicates to run");
    struct arg_file *infiles1 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_file *output1 = arg_file0("o", "output", "output-file",
            "Output HDF5 file");
    struct arg_end *end1 = arg_end(20);
    void* argtable1[] = {cmd1, verbose1, infiles1, output1, replicates1, end1};
    int nerrors1;

    /* SYNTAX 2: ld [-v] <input-file> */
    struct arg_rex *cmd2 = arg_rex1(NULL, NULL, "ld", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose2 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles2 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end2 = arg_end(20);
    void* argtable2[] = {cmd2, verbose2, infiles2, end2};
    int nerrors2;

    /* SYNTAX 3: haplotypes [-v] <input-file> */
    struct arg_rex *cmd3 = arg_rex1(NULL, NULL, "haplotypes", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose3 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles3 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end3 = arg_end(20);
    void* argtable3[] = {cmd3, verbose3, infiles3, end3};
    int nerrors3;

    /* SYNTAX 4: variants [-v] <input-file> */
    struct arg_rex *cmd4 = arg_rex1(NULL, NULL, "variants", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose4 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles4 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end4 = arg_end(20);
    void* argtable4[] = {cmd4, verbose4, infiles4, end4};
    int nerrors4;

    /* SYNTAX 5: vcf [-v] <input-file> */
    struct arg_rex *cmd5 = arg_rex1(NULL, NULL, "vcf", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose5 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles5 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end5 = arg_end(20);
    void* argtable5[] = {cmd5, verbose5, infiles5, end5};
    int nerrors5;

    /* SYNTAX 6: print  [-v] <input-file> */
    struct arg_rex *cmd6 = arg_rex1(NULL, NULL, "print", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose6 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles6 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end6 = arg_end(20);
    void* argtable6[] = {cmd6, verbose6, infiles6, end6};
    int nerrors6;

    /* SYNTAX 7: newick [-v] <input-file> */
    struct arg_rex *cmd7 = arg_rex1(NULL, NULL, "newick", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose7 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles7 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end7 = arg_end(20);
    void* argtable7[] = {cmd7, verbose7, infiles7, end7};
    int nerrors7;

    /* SYNTAX 8: stats [-v] <input-file> */
    struct arg_rex *cmd8 = arg_rex1(NULL, NULL, "stats", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose8 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles8 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end8 = arg_end(20);
    void* argtable8[] = {cmd8, verbose8, infiles8, end8};
    int nerrors8;

    /* SYNTAX 9: simplify [-v] <input-file> */
    struct arg_rex *cmd9 = arg_rex1(NULL, NULL, "simplify", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose9 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles9 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_file *outfiles9 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end9 = arg_end(20);
    void* argtable9[] = {cmd9, verbose9, infiles9, outfiles9, end9};
    int nerrors9;

    int exitcode = EXIT_SUCCESS;
    const char *progname = "main";

    /* Set defaults */
    replicates1->ival[0] = 1;
    output1->filename[0] = NULL;

    nerrors1 = arg_parse(argc, argv, argtable1);
    nerrors2 = arg_parse(argc, argv, argtable2);
    nerrors3 = arg_parse(argc, argv, argtable3);
    nerrors4 = arg_parse(argc, argv, argtable4);
    nerrors5 = arg_parse(argc, argv, argtable5);
    nerrors6 = arg_parse(argc, argv, argtable6);
    nerrors7 = arg_parse(argc, argv, argtable7);
    nerrors8 = arg_parse(argc, argv, argtable8);
    nerrors9 = arg_parse(argc, argv, argtable9);

    if (nerrors1 == 0) {
        run_simulate(infiles1->filename[0], output1->filename[0], verbose1->count,
                replicates1->ival[0]);
    } else if (nerrors2 == 0) {
        run_ld(infiles2->filename[0], verbose2->count);
    } else if (nerrors3 == 0) {
        run_haplotypes(infiles3->filename[0], verbose3->count);
    } else if (nerrors4 == 0) {
        run_variants(infiles4->filename[0], verbose4->count);
    } else if (nerrors5 == 0) {
        run_vcf(infiles5->filename[0], verbose5->count);
    } else if (nerrors6 == 0) {
        run_print(infiles6->filename[0], verbose6->count);
    } else if (nerrors7 == 0) {
        run_newick(infiles7->filename[0], verbose7->count);
    } else if (nerrors8 == 0) {
        run_stats(infiles8->filename[0], verbose8->count);
    } else if (nerrors9 == 0) {
        run_simplify(infiles9->filename[0], outfiles9->filename[0], verbose9->count);
    } else {
        /* We get here if the command line matched none of the possible syntaxes */
        if (cmd1->count > 0) {
            arg_print_errors(stdout, end1, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable1, "\n");
        } else if (cmd2->count > 0) {
            arg_print_errors(stdout, end2, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable2, "\n");
        } else if (cmd3->count > 0) {
            arg_print_errors(stdout, end3, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable3, "\n");
        } else if (cmd4->count > 0) {
            arg_print_errors(stdout, end4, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable4, "\n");
        } else if (cmd5->count > 0) {
            arg_print_errors(stdout, end5, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable5, "\n");
        } else if (cmd6->count > 0) {
            arg_print_errors(stdout, end6, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable6, "\n");
        } else if (cmd7->count > 0) {
            arg_print_errors(stdout, end7, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable7, "\n");
        } else if (cmd8->count > 0) {
            arg_print_errors(stdout, end8, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable8, "\n");
        } else if (cmd9->count > 0) {
            arg_print_errors(stdout, end9, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable9, "\n");
        } else {
            /* no correct cmd literals were given, so we cant presume which syntax was intended */
            printf("%s: missing command.\n",progname);
            printf("usage 1: %s ", progname);  arg_print_syntax(stdout, argtable1, "\n");
            printf("usage 2: %s ", progname);  arg_print_syntax(stdout, argtable2, "\n");
            printf("usage 3: %s ", progname);  arg_print_syntax(stdout, argtable3, "\n");
            printf("usage 4: %s ", progname);  arg_print_syntax(stdout, argtable4, "\n");
            printf("usage 5: %s ", progname);  arg_print_syntax(stdout, argtable5, "\n");
            printf("usage 6: %s ", progname);  arg_print_syntax(stdout, argtable6, "\n");
            printf("usage 7: %s ", progname);  arg_print_syntax(stdout, argtable7, "\n");
            printf("usage 8: %s ", progname);  arg_print_syntax(stdout, argtable8, "\n");
            printf("usage 9: %s ", progname);  arg_print_syntax(stdout, argtable9, "\n");
        }
    }

    arg_freetable(argtable1, sizeof(argtable1) / sizeof(argtable1[0]));
    arg_freetable(argtable2, sizeof(argtable2) / sizeof(argtable2[0]));
    arg_freetable(argtable3, sizeof(argtable3) / sizeof(argtable3[0]));
    arg_freetable(argtable4, sizeof(argtable4) / sizeof(argtable4[0]));
    arg_freetable(argtable5, sizeof(argtable5) / sizeof(argtable5[0]));
    arg_freetable(argtable6, sizeof(argtable6) / sizeof(argtable6[0]));
    arg_freetable(argtable7, sizeof(argtable7) / sizeof(argtable7[0]));
    arg_freetable(argtable8, sizeof(argtable8) / sizeof(argtable8[0]));
    arg_freetable(argtable9, sizeof(argtable9) / sizeof(argtable9[0]));

    return exitcode;
}
