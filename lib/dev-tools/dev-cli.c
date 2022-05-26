/*
** Copyright (C) 2015-2019 University of Oxford
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
#include "util.h"

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
fatal_tskit_error(int err, int line)
{
    fprintf(stderr, "error line %d::", line);
    fprintf(stderr, ":%d:'%s'\n", err, tsk_strerror(err));
    exit(EXIT_FAILURE);
}

static void
fatal_msprime_error(int err, int line)
{
    fprintf(stderr, "error line %d::", line);
    fprintf(stderr, ":%d:'%s'\n", err, msp_strerror(err));
    exit(EXIT_FAILURE);
}

static bool
read_pedigree(tsk_table_collection_t *tables, config_t *config)
{
    int ret = 0;
    size_t i, j, k, size, num_inds;
    tsk_flags_t flags;
    tsk_id_t ind_id;
    size_t ploidy = 2;
    tsk_id_t *parents = NULL;
    double *time = NULL;
    tsk_flags_t *is_sample = NULL;
    config_setting_t *s, *row;
    config_setting_t *setting;

    setting = config_lookup(config, "time");
    if (config_setting_is_array(setting) == CONFIG_FALSE) {
        fatal_error("Pedigree times must be an array");
    }
    size = (size_t) config_setting_length(setting);
    num_inds = size;

    if (num_inds == 0) {
        /* If we have zero individuals, skip the rest of this method. */
        return false;
    }

    parents = malloc(num_inds * ploidy * sizeof(tsk_id_t));
    time = malloc(num_inds * sizeof(double));
    is_sample = malloc(num_inds * sizeof(tsk_flags_t));
    if (parents == NULL || time == NULL || is_sample == NULL) {
        fatal_error("Out of memory");
    }

    setting = config_lookup(config, "parents");
    if (config_setting_is_list(setting) == CONFIG_FALSE) {
        fatal_error("Pedigree parents must be an list");
    }
    for (j = 0; j < num_inds; j++) {
        row = config_setting_get_elem(setting, (unsigned int) j);
        if (row == NULL) {
            fatal_error("error reading pedigree parents %d", j);
        }
        if (config_setting_is_array(row) == CONFIG_FALSE) {
            fatal_error("pedigree parents row %d not an array", j);
        }
        if (config_setting_length(row) != ploidy) {
            fatal_error("pedigree parents must have 'ploidy' columns");
        }
        for (i = 0; i < ploidy; i++) {
            s = config_setting_get_elem(row, (unsigned int) i);
            if (!config_setting_is_number(s)) {
                fatal_error("pedigree parent entries must be numbers");
            }
            parents[j * ploidy + i] = (tsk_id_t) config_setting_get_int(s);
        }
    }

    setting = config_lookup(config, "time");
    if (config_setting_is_array(setting) == CONFIG_FALSE) {
        fatal_error("Pedigree times must be an array");
    }
    for (j = 0; j < num_inds; j++) {
        s = config_setting_get_elem(setting, (unsigned int) j);
        if (s == NULL) {
            fatal_error("error reading pedigree individual %d", j);
        }
        time[j] = (double) config_setting_get_float(s);
    }

    setting = config_lookup(config, "is_sample");
    if (config_setting_is_array(setting) == CONFIG_FALSE) {
        fatal_error("Pedigree sample flags must be an array");
    }
    for (j = 0; j < num_inds; j++) {
        s = config_setting_get_elem(setting, (unsigned int) j);
        if (s == NULL) {
            fatal_error("error reading pedigree individual %d", j);
        }
        is_sample[j] = (tsk_flags_t) config_setting_get_int(s);
        assert(is_sample[j] == 0 || is_sample[j] == 1);
    }

    /* Update the tables to store this information */
    for (j = 0; j < num_inds; j++) {
        ret = tsk_individual_table_add_row(&tables->individuals, 0, NULL, 0,
            /* encode the parents in the metadata */
            NULL, 0, (void *) (parents + j * ploidy), ploidy * sizeof(tsk_id_t));
        if (ret < 0) {
            fatal_tskit_error(ret, __LINE__);
        }
        ind_id = ret;
        for (k = 0; k < ploidy; k++) {
            flags = is_sample[j] ? TSK_NODE_IS_SAMPLE : 0;
            ret = tsk_node_table_add_row(
                &tables->nodes, flags, time[j], 0, ind_id, NULL, 0);
            if (ret < 0) {
                fatal_tskit_error(ret, __LINE__);
            }
        }
    }

    free(parents);
    free(time);
    free(is_sample);
    return true;
}
static void
read_initial_tables(config_t *config, tsk_table_collection_t *tables)
{
    int ret;
    size_t j, n, num_populations;
    double time;
    tsk_id_t population;
    tsk_flags_t flags = TSK_NODE_IS_SAMPLE;
    config_setting_t *s, *t;
    config_setting_t *setting = config_lookup(config, "samples");
    bool pedigree_specified;

    if (setting == NULL) {
        fatal_error("samples is a required parameter");
    }
    if (config_setting_is_list(setting) == CONFIG_FALSE) {
        fatal_error("samples must be a list");
    }
    pedigree_specified = read_pedigree(tables, config);

    n = (size_t) config_setting_length(setting);
    if (pedigree_specified && n > 0) {
        fatal_error("Cannot specify samples and pedigree");
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
        population = (population_id_t) config_setting_get_int(t);
        t = config_setting_get_member(s, "time");
        if (t == NULL) {
            fatal_error("population not specified");
        }
        time = config_setting_get_float(t);
        ret = tsk_node_table_add_row(
            &tables->nodes, flags, time, population, TSK_NULL, NULL, 0);
        if (ret < 0) {
            fatal_tskit_error(ret, __LINE__);
        }
    }
    setting = config_lookup(config, "population_configuration");
    if (setting == NULL) {
        fatal_error("population_configuration is a required parameter");
    }
    if (config_setting_is_list(setting) == CONFIG_FALSE) {
        fatal_error("population_configuration must be a list");
    }
    num_populations = config_setting_length(setting);

    for (j = 0; j < num_populations; j++) {
        ret = tsk_population_table_add_row(&tables->populations, NULL, 0);
        if (ret < 0) {
            fatal_tskit_error(ret, __LINE__);
        }
    }
}

static void
read_sweep_genic_selection_model_config(msp_t *msp, config_setting_t *model_setting)
{
    int ret = 0;
    config_setting_t *s;
    double position, start_frequency, end_frequency, alpha, dt;

    s = config_setting_get_member(model_setting, "position");
    if (s == NULL) {
        fatal_error("sweep genic selection model position not specified");
    }
    position = config_setting_get_float(s);

    s = config_setting_get_member(model_setting, "start_frequency");
    if (s == NULL) {
        fatal_error("sweep genic selection model start_frequency not specified");
    }
    start_frequency = config_setting_get_float(s);

    s = config_setting_get_member(model_setting, "end_frequency");
    if (s == NULL) {
        fatal_error("sweep genic selection model end_frequency not specified");
    }
    end_frequency = config_setting_get_float(s);

    s = config_setting_get_member(model_setting, "dt");
    if (s == NULL) {
        fatal_error("sweep genic selection model dt not specified");
    }
    dt = config_setting_get_float(s);

    s = config_setting_get_member(model_setting, "alpha");
    if (s == NULL) {
        fatal_error("sweep genic selection model alpha not specified");
    }
    alpha = config_setting_get_float(s);

    ret = msp_set_simulation_model_sweep_genic_selection(
        msp, position, start_frequency, end_frequency, alpha, dt);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
}

static void
read_model_config(msp_t *msp, config_t *config)
{
    int ret = 0;
    config_setting_t *setting = config_lookup(config, "model");
    config_setting_t *s;
    const char *name;
    double psi, c, alpha, truncation_point;

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

    /* printf("Checking model specification\n"); */
    if (strcmp(name, "hudson") == 0) {
        ret = msp_set_simulation_model_hudson(msp);
    } else if (strcmp(name, "smc") == 0) {
        ret = msp_set_simulation_model_smc(msp);
    } else if (strcmp(name, "smc_prime") == 0) {
        ret = msp_set_simulation_model_smc_prime(msp);
    } else if (strcmp(name, "dtwf") == 0) {
        ret = msp_set_simulation_model_dtwf(msp);
    } else if (strcmp(name, "fixed_pedigree") == 0) {
        ret = msp_set_simulation_model_fixed_pedigree(msp);
    } else if (strcmp(name, "dirac") == 0) {
        s = config_setting_get_member(setting, "psi");
        if (s == NULL) {
            fatal_error("dirac model psi not specified");
        }
        psi = config_setting_get_float(s);
        s = config_setting_get_member(setting, "c");
        if (s == NULL) {
            fatal_error("dirac model c not specified");
        }
        c = config_setting_get_float(s);
        ret = msp_set_simulation_model_dirac(msp, psi, c);
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
    } else if (strcmp(name, "sweep_genic_selection") == 0) {
        read_sweep_genic_selection_model_config(msp, setting);
    } else {
        fatal_error("Unknown simulation model '%s'", name);
    }
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
}

static void
read_population_configuration(msp_t *msp, config_t *config)
{
    int ret = 0;
    int j;
    double growth_rate, initial_size;
    int num_populations, num_labels;
    config_setting_t *s, *t;
    config_setting_t *setting = config_lookup(config, "population_configuration");

    if (setting == NULL) {
        fatal_error("population_configuration is a required parameter");
    }
    if (config_setting_is_list(setting) == CONFIG_FALSE) {
        fatal_error("population_configuration must be a list");
    }
    num_populations = config_setting_length(setting);
    if (config_lookup_int(config, "num_labels", &num_labels) == CONFIG_FALSE) {
        fatal_error("num_labels is a required parameter");
    }
    ret = msp_set_num_labels(msp, num_labels);
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
        ret = msp_set_population_configuration(msp, j, initial_size, growth_rate, true);
        if (ret != 0) {
            fatal_msprime_error(ret, __LINE__);
        }
    }
}

static void
read_demographic_events(msp_t *msp, config_t *config)
{
    int ret = 0;
    int j;
    const char *type;
    double time, growth_rate, initial_size, migration_rate, proportion, intensity,
        strength;
    int num_demographic_events, population_id, source, dest;
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
            ret = msp_add_population_parameters_change(
                msp, time, population_id, initial_size, growth_rate);
        } else if (strcmp(type, "migration_rate_change") == 0) {
            t = config_setting_get_member(s, "migration_rate");
            if (t == NULL) {
                fatal_error("migration_rate not specified");
            }
            migration_rate = config_setting_get_float(t);
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
            ret = msp_add_migration_rate_change(msp, time, source, dest, migration_rate);
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
            ret = msp_add_instantaneous_bottleneck(msp, time, population_id, strength);
        } else {
            fatal_error("unknown demographic event type '%s'", type);
        }
        if (ret != 0) {
            fatal_msprime_error(ret, __LINE__);
        }
    }
}

static void
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
        fatal_error("Out of memory");
    }
    for (j = 0; j < size; j++) {
        s = config_setting_get_elem(setting, (unsigned int) j);
        if (s == NULL) {
            fatal_error("error reading migration_matrix[%d]", j);
        }
        migration_matrix[j] = (double) config_setting_get_float(s);
    }
    ret = msp_set_migration_matrix(msp, size, migration_matrix);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    free(migration_matrix);
}

static void
read_recomb_map(rate_map_t *recomb_map, config_t *config)
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
        fatal_error("Out of memory");
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
            fatal_error("recombination_map entries must be [coord, rate] pairs");
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
    ret = rate_map_alloc(recomb_map, size - 1, coordinates, rates);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    free(rates);
    free(coordinates);
}

static void
get_configuration(gsl_rng *rng, msp_t *msp, tsk_table_collection_t *tables,
    mutation_params_t *mutation_params, rate_map_t *recomb_map, const char *filename)
{
    int ret = 0;
    int err;
    int int_tmp;
    double start_time, gene_conversion_rate, gene_conversion_tract_length;
    const char *from_ts_path;
    config_t *config = malloc(sizeof(config_t));
    config_setting_t *t;

    if (config == NULL) {
        fatal_error("no memory");
    }
    config_init(config);
    err = config_read_file(config, filename);
    if (err == CONFIG_FALSE) {
        fatal_error("configuration error:%s at line %d in file %s\n",
            config_error_text(config), config_error_line(config), filename);
    }
    if (config_lookup_int(config, "random_seed", &int_tmp) == CONFIG_FALSE) {
        fatal_error("random_seed is a required parameter");
    }
    gsl_rng_set(rng, (unsigned long) int_tmp);
    if (config_lookup_string(config, "from", &from_ts_path) == CONFIG_TRUE) {
        err = tsk_table_collection_load(tables, from_ts_path, TSK_NO_INIT);
        if (err != 0) {
            fatal_tskit_error(err, __LINE__);
        }
    } else {
        read_initial_tables(config, tables);
    }

    read_recomb_map(recomb_map, config);

    tables->sequence_length = rate_map_get_sequence_length(recomb_map);

    ret = msp_alloc(msp, tables, rng);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    ret = msp_set_recombination_map(
        msp, recomb_map->size, recomb_map->position, recomb_map->rate);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }

    if (config_lookup_float(config, "gene_conversion_rate", &gene_conversion_rate)
        == CONFIG_FALSE) {
        fatal_error("gene_conversion_rate is a required parameter");
    }
    if (config_lookup_float(
            config, "gene_conversion_tract_length", &gene_conversion_tract_length)
        == CONFIG_FALSE) {
        fatal_error("gene_conversion_tract_length is a required parameter");
    }
    ret = msp_set_gene_conversion_rate(msp, gene_conversion_rate);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    ret = msp_set_gene_conversion_tract_length(msp, gene_conversion_tract_length);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    if (config_lookup_float(config, "mutation_rate", &mutation_params->mutation_rate)
        == CONFIG_FALSE) {
        fatal_error("mutation_rate is a required parameter");
    }
    if (config_lookup_int(config, "mutation_alphabet", &int_tmp) == CONFIG_FALSE) {
        fatal_error("mutation_alphabet is a required parameter.");
    }
    mutation_params->alphabet = int_tmp;
    if (config_lookup_int(config, "avl_node_block_size", &int_tmp) == CONFIG_FALSE) {
        fatal_error("avl_node_block_size is a required parameter");
    }
    ret = msp_set_avl_node_block_size(msp, (size_t) int_tmp);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    if (config_lookup_int(config, "segment_block_size", &int_tmp) == CONFIG_FALSE) {
        fatal_error("segment_block_size is a required parameter");
    }
    ret = msp_set_segment_block_size(msp, (size_t) int_tmp);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    if (config_lookup_int(config, "node_mapping_block_size", &int_tmp) == CONFIG_FALSE) {
        fatal_error("node_mapping_block_size is a required parameter");
    }
    ret = msp_set_node_mapping_block_size(msp, (size_t) int_tmp);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    if (config_lookup_int(config, "store_migrations", &int_tmp) == CONFIG_FALSE) {
        fatal_error("store_migrations is a required parameter");
    }
    ret = msp_set_store_migrations(msp, (bool) int_tmp);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    if (config_lookup_int(config, "store_full_arg", &int_tmp) == CONFIG_FALSE) {
        fatal_error("store_ful_arg is a required parameter");
    }
    if (config_lookup_float(config, "start_time", &start_time) == CONFIG_FALSE) {
        fatal_error("start_time is a required parameter");
    }
    ret = msp_set_start_time(msp, start_time);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    ret = msp_set_store_full_arg(msp, (bool) int_tmp);
    if (ret != 0) {
        fatal_error(msp_strerror(ret));
    }
    t = config_lookup(config, "model");
    if (t == NULL) {
        fatal_error("model not specified");
    }
    read_model_config(msp, config);
    read_population_configuration(msp, config);
    read_migration_matrix(msp, config);
    read_demographic_events(msp, config);
    config_destroy(config);
    free(config);
}

static void
record_provenance(tsk_provenance_table_t *provenance)
{
    time_t timer;
    size_t timestamp_size = 64;
    char buffer[timestamp_size];
    struct tm *tm_info;
    const char *provenance_str = "{\"program\"=\"main\"}";
    int ret;

    time(&timer);
    tm_info = localtime(&timer);
    strftime(buffer, timestamp_size, "%Y-%m-%dT%H:%M:%S", tm_info);

    ret = tsk_provenance_table_add_row(
        provenance, buffer, strlen(buffer), provenance_str, strlen(provenance_str));
    if (ret < 1) {
        fatal_error("Error recording provenance");
    }
}

static void
run_simulate(
    const char *conf_file, const char *output_file, int verbose, int num_replicates)
{
    int ret = -1;
    int j;
    mutation_params_t mutation_params;
    msp_t msp;
    rate_map_t recomb_map;
    mutation_model_t mut_model;
    mutgen_t mutgen;
    tsk_table_collection_t tables;
    tsk_treeseq_t tree_seq;

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    if (rng == NULL) {
        fatal_error("No memory");
    }
    ret = tsk_table_collection_init(&tables, 0);
    if (ret != 0) {
        fatal_tskit_error(ret, __LINE__);
    }
    get_configuration(rng, &msp, &tables, &mutation_params, &recomb_map, conf_file);

    ret = msp_initialise(&msp);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }

    ret = matrix_mutation_model_factory(&mut_model, mutation_params.alphabet);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    ret = mutgen_alloc(&mutgen, rng, &tables, &mut_model, 0);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    ret = mutgen_set_rate(&mutgen, mutation_params.mutation_rate);
    if (ret != 0) {
        fatal_msprime_error(ret, __LINE__);
    }
    record_provenance(&tables.provenances);

    msp_print_state(&msp, stdout);
    msp_verify(&msp, MSP_VERIFY_BREAKPOINTS);
    for (j = 0; j < num_replicates; j++) {
        if (verbose >= 1) {
            printf("=====================\n");
            printf("replicate %d\n", j);
            printf("=====================\n");
        }
        ret = msp_reset(&msp);
        if (ret != 0) {
            fatal_msprime_error(ret, __LINE__);
        }
        msp_verify(&msp, 0);
        ret = msp_run(&msp, DBL_MAX, UINT32_MAX);
        if (ret < 0) {
            fatal_msprime_error(ret, __LINE__);
        }
        if (verbose >= 1) {
            msp_print_state(&msp, stdout);
        }
        msp_verify(&msp, 0);
        ret = msp_finalise_tables(&msp);

        if (ret != 0) {
            fatal_msprime_error(ret, __LINE__);
        }
        ret = mutgen_generate(&mutgen, 0);
        if (ret != 0) {
            fatal_msprime_error(ret, __LINE__);
        }
        ret = tsk_treeseq_init(&tree_seq, &tables, TSK_TS_INIT_BUILD_INDEXES);
        if (ret != 0) {
            fatal_tskit_error(ret, __LINE__);
        }
        if (output_file != NULL) {
            ret = tsk_treeseq_dump(&tree_seq, output_file, 0);
            if (ret != 0) {
                fatal_tskit_error(ret, __LINE__);
            }
        }
        if (verbose >= 1) {
            tsk_table_collection_print_state(&tables, stdout);
            printf("-----------------\n");
            mutgen_print_state(&mutgen, stdout);
            printf("-----------------\n");
            tsk_treeseq_print_state(&tree_seq, stdout);
        }
        tsk_treeseq_free(&tree_seq);
    }

    msp_free(&msp);
    rate_map_free(&recomb_map);
    mutation_model_free(&mut_model);
    mutgen_free(&mutgen);
    gsl_rng_free(rng);
    tsk_table_collection_free(&tables);
}

int
main(int argc, char **argv)
{
    /* SYNTAX 1: simulate [-v] <config-file> -o <output-file> */
    struct arg_rex *cmd1 = arg_rex1(NULL, NULL, "simulate", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose1 = arg_lit0("v", "verbose", NULL);
    struct arg_int *replicates1
        = arg_int0("r", "replicates", "<num-replicates>", "number of replicates to run");
    struct arg_file *infiles1 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_file *output1
        = arg_file0("o", "output", "output-file", "Output trees file");
    struct arg_end *end1 = arg_end(20);
    void *argtable1[] = { cmd1, verbose1, infiles1, output1, replicates1, end1 };
    int nerrors1;

    int exitcode = EXIT_SUCCESS;
    const char *progname = "main";

    /* Set defaults */
    replicates1->ival[0] = 1;
    output1->filename[0] = NULL;

    nerrors1 = arg_parse(argc, argv, argtable1);

    if (nerrors1 == 0) {
        run_simulate(infiles1->filename[0], output1->filename[0], verbose1->count,
            replicates1->ival[0]);
    } else {
        /* We get here if the command line matched none of the possible syntaxes */
        if (cmd1->count > 0) {
            arg_print_errors(stdout, end1, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable1, "\n");
        } else {
            /* no correct cmd literals were given, so we cant presume which syntax was
             * intended */
            printf("%s: missing command.\n", progname);
            printf("usage 1: %s ", progname);
            arg_print_syntax(stdout, argtable1, "\n");
        }
    }

    arg_freetable(argtable1, sizeof(argtable1) / sizeof(argtable1[0]));
    return exitcode;
}
