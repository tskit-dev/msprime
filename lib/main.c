/*
** Copyright (C) 2014 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

#include <libconfig.h>

#include "msprime.h"

#include "util.h"

static void
msp_read_population_models(msp_t *self, config_t *config)
{
    int j;
    const char *type;
    double time, param;
    unsigned int num_population_models;
    config_setting_t *s, *t;
    config_setting_t *setting = config_lookup(config, "population_models");

    if (setting == NULL) {
        fatal_error("population_models is a required parameter");
    }
    if (config_setting_is_list(setting) == CONFIG_FALSE) {
        fatal_error("population_models must be a list");
    }
    num_population_models = config_setting_length(setting);
    for (j = 0; j < num_population_models; j++) {
        s = config_setting_get_elem(setting, j);
        if (s == NULL) {
            fatal_error("error reading population_models[%d]", j);
        }
        if (config_setting_is_group(s) == CONFIG_FALSE) {
            fatal_error("population_models[%d] not a group", j);
        }
        t = config_setting_get_member(s, "time");
        if (t == NULL) {
            fatal_error("time not specified");
        }
        time = config_setting_get_float(t);
        if (time < 0.0) {
            fatal_error("population_model time must be > 0");
        }
        t = config_setting_get_member(s, "param");
        if (t == NULL) {
            fatal_error("param not specified");
        }
        param = config_setting_get_float(t);
        t = config_setting_get_member(s, "type");
        if (t == NULL) {
            fatal_error("type not specified");
        }
        type = config_setting_get_string(t);
        if (strcmp(type, "constant") == 0) {
            msp_add_constant_population_model(self, time, param);
        } else if (strcmp(type, "exponential") == 0) {
            msp_add_exponential_population_model(self, time, param);
        } else {
            fatal_error("unknown population_model type '%s'", type);
        }
    }
}

static void
msp_read_config(msp_t *self, const char *filename)
{
    int err;
    int tmp;
    size_t s;
    const char *str;
    config_t *config = xmalloc(sizeof(config_t));

    config_init(config);
    err = config_read_file(config, filename);
    if (err == CONFIG_FALSE) {
        fatal_error("configuration error:%s at line %d in file %s\n",
                config_error_text(config), config_error_line(config),
                filename);
    }
    if (config_lookup_int(config, "verbosity", &tmp) == CONFIG_FALSE) {
        fatal_error("verbosity is a required parameter");
    }
    self->verbosity = tmp;
    if (config_lookup_int(config, "approx", &tmp) == CONFIG_FALSE) {
        fatal_error("approx is a required parameter");
    }
    self->approx = tmp;
    if (config_lookup_int(config, "sample_size", &tmp) == CONFIG_FALSE) {
        fatal_error("sample_size is a required parameter");
    }
    self->sample_size = tmp;
    if (config_lookup_int(config, "num_loci", &tmp) == CONFIG_FALSE) {
        fatal_error("num_loci is a required parameter");
    }
    self->num_loci = tmp;
    if (config_lookup_int(config, "max_avl_nodes", &tmp) == CONFIG_FALSE) {
        fatal_error("max_avl_nodes is a required parameter");
    }
    self->max_avl_nodes = tmp;
    if (config_lookup_int(config, "max_segments", &tmp) == CONFIG_FALSE) {
        fatal_error("max_segments is a required parameter");
    }
    self->max_segments = tmp;
    if (config_lookup_int(config, "max_trees", &tmp) == CONFIG_FALSE) {
        fatal_error("max_trees is a required parameter");
    }
    self->max_trees = tmp;
    if (config_lookup_int(config, "max_coalescence_records", &tmp)
            == CONFIG_FALSE) {
        fatal_error("max_coalescence_records is a required parameter");
    }
    self->max_coalescence_records = tmp;
    if (config_lookup_float(config, "recombination_rate",
            &self->recombination_rate) == CONFIG_FALSE) {
        fatal_error("recombination_rate is a required parameter");
    }
    if (config_lookup_string(config, "coalescence_record_file", &str)
            == CONFIG_FALSE) {
        fatal_error("coalescence_record_file is a required parameter");
    }
    s = strlen(str);
    self->coalescence_record_file_pattern = xmalloc(s + 1);
    strcpy(self->coalescence_record_file_pattern, str);
    msp_read_population_models(self, config);
    config_destroy(config);
    free(config);
}

static void
run_simulate(char *conf_file, long seed, long num_replicates)
{
    long j;
    msp_t *self = xmalloc(sizeof(msp_t));

    self->random_seed = seed;
    /* Add default population model; TODO this is error prone*/
    self->population_models = NULL;
    msp_add_constant_population_model(self, -1.0, 1.0);
    msp_read_config(self, conf_file);
    msp_alloc(self);
    printf("t\tnum_trees\tre_events\tca_events\tcoalescence_records\tmax_segments"
            "\tmax_population_size\n");
    for (j = 0; j < num_replicates; j++) {
        msp_simulate(self);
        msp_reset(self);
    }
    msp_free(self);
    free(self);
}

/*
 * Tree file reader. This is a very simple reader for the tree files that
 * inputs the binary format and writes out the trees to stdout in
 * oriented forest form.
 */

static void
run_process(char *tree_file)
{
    size_t file_size, ret;
    coalescence_record_t *coalescence_records, *cr;
    unsigned int header[2], sample_size, num_loci, j, k;
    int *pi;
    float *tau;
    FILE *f = fopen(tree_file, "r");
    if (f == NULL) {
        fatal_error("cannot open %s: %s", tree_file, strerror(errno));
    }
    /* read the header */
    ret = fread(header, sizeof(header), 1, f);
    if (ret != 1) {
        fatal_error("error reading %s: %s", tree_file, strerror(errno));
    }
    if (header[0] != TREEFILE_MAGIC) {
        fatal_error("error reading %s: magic number mismatchs", tree_file);
    }
    sample_size = header[1];
    /* allocate the tree */
    pi = xmalloc(2 * sample_size * sizeof(int));
    tau = xmalloc(2 * sample_size * sizeof(float));
    /* Allocate the coalescence records. */
    fseek(f, 0, SEEK_END);
    file_size = ftell(f);
    file_size -= sizeof(header);
    coalescence_records = xmalloc(file_size);
    fseek(f, sizeof(header), SEEK_SET);
    ret = fread(coalescence_records, 1, file_size, f);
    if (ret != file_size) {
        fatal_error("error reading %s: %s", tree_file, strerror(errno));
    }
    num_loci = 0;
    for (cr = coalescence_records; cr->left != 0; cr++) {
        if (cr->right > num_loci) {
            num_loci = cr->right;
        }
    }
    /* Now go through each locus and print out the corresponding tree.
     * A more sophisticated implementation would only write out trees
     * at the breakpoints, which can be easily found by getting the
     * unique left values among the coalescence records
     */
    for (j = 1; j <= num_loci; j++) {
        memset(pi, 0, 2 * sample_size * sizeof(int));
        memset(tau, 0, 2 * sample_size * sizeof(float));
        for (cr = coalescence_records; cr->left != 0; cr++) {
            if (cr->left <= j && j <= cr->right) {
                pi[cr->children[0]] = cr->parent;
                pi[cr->children[1]] = cr->parent;
                tau[cr->parent] = cr->time;
            }
        }
        printf("%d\t", j);
        for (k = 1; k < 2 * sample_size; k++) {
            printf(" %3d", pi[k]);
        }
        printf(" :: ");
        for (k = 1; k < 2 * sample_size; k++) {
            printf(" %.3f", tau[k]);
        }
        printf("\n");
    }

    if (fclose(f) != 0) {
        fatal_error("cannot close %s: %s", tree_file, strerror(errno));
    }
    free(pi);
    free(tau);
    free(coalescence_records);
}



int
main(int argc, char** argv)
{
    char *cmd;
    long num_replicates, seed;
    if (argc < 2) {
        fatal_error("usage: %s simulate|process [ARGS]", argv[0]);
    }
    cmd = argv[1];
    if (strcmp("simulate", cmd) == 0) {
        if (argc != 5) {
            fatal_error("usage: %s simulate CONFIG_FILE SEED NUM_REPS", argv[0]);
        }
        if (parse_long(argv[3], &seed, 0, LONG_MAX) != 0) {
            fatal_error("cannot parse seed '%s'", argv[3]);
        }
        if (parse_long(argv[4], &num_replicates, 0, LONG_MAX) != 0) {
            fatal_error("cannot parse replicates '%s'", argv[4]);
        }
        run_simulate(argv[2], seed, num_replicates);
    } else if (strcmp("process", cmd) == 0) {
        if (argc != 3) {
            fatal_error("usage: %s process TREE_FILE", argv[0]);
        }
        run_process(argv[2]);
    } else {
        fatal_error("unrecognised command '%s'", cmd);
    }

    return EXIT_SUCCESS;
}
