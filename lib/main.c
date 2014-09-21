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
#include <stdarg.h>
#include <float.h>

#include <libconfig.h>

#include "msprime.h"

void
fatal_error(const char *msg, ...)
{
    va_list argp;
    fprintf(stderr, "sms:");
    va_start(argp, msg);
    vfprintf(stderr, msg, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}

/*
 * Parses the specified string into a long and assigns the value into
 * the specified pointer. Returns EINVAL if the string cannot be
 * converted to double or if min <= x <= max does not hold; returns 0 if
 * the value is converted successfully.
 */
int
parse_long(const char *str, long *value, const long min,
        const long max)
{
    int ret = 0;
    long x;
    char *tail;
    x = strtol(str, &tail, 10);
    if (tail[0] != '\0') {
        ret = EINVAL;
    } else if (min > x || max < x) {
        ret = EINVAL;
    } else {
        *value = x;
    }
    return ret;
}


static void
read_population_models(msp_t *msp, config_t *config)
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
            msp_add_constant_population_model(msp, time, param);
        } else if (strcmp(type, "exponential") == 0) {
            msp_add_exponential_population_model(msp, time, param);
        } else {
            fatal_error("unknown population_model type '%s'", type);
        }
    }
}

static void
read_config(msp_t *msp, const char *filename)
{
    int err;
    int tmp;
    size_t s;
    const char *str;
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
    if (config_lookup_int(config, "sample_size", &tmp) == CONFIG_FALSE) {
        fatal_error("sample_size is a required parameter");
    }
    msp->sample_size = tmp;
    if (config_lookup_int(config, "num_loci", &tmp) == CONFIG_FALSE) {
        fatal_error("num_loci is a required parameter");
    }
    msp->num_loci = tmp;
    if (config_lookup_int(config, "avl_node_block_size", &tmp) == CONFIG_FALSE) {
        fatal_error("avl_node_block_size is a required parameter");
    }
    msp->avl_node_block_size = tmp;
    if (config_lookup_int(config, "segment_block_size", &tmp) == CONFIG_FALSE) {
        fatal_error("segment_block_size is a required parameter");
    }
    msp->segment_block_size = tmp;
    if (config_lookup_int(config, "node_mapping_block_size", &tmp)
            == CONFIG_FALSE) {
        fatal_error("node_mapping_block_size is a required parameter");
    }
    msp->node_mapping_block_size = tmp;
    if (config_lookup_int(config, "max_memory", &tmp)
            == CONFIG_FALSE) {
        fatal_error("max_memory is a required parameter");
    }
    msp->max_memory = tmp * 1024 * 1024;
    if (config_lookup_float(config, "recombination_rate",
            &msp->recombination_rate) == CONFIG_FALSE) {
        fatal_error("recombination_rate is a required parameter");
    }
    if (config_lookup_string(config, "coalescence_record_file", &str)
            == CONFIG_FALSE) {
        fatal_error("coalescence_record_file is a required parameter");
    }
    s = strlen(str);
    msp->coalescence_record_filename = malloc(s + 1);
    if (msp->coalescence_record_filename == NULL) {
        fatal_error("no memory");
    }
    strcpy(msp->coalescence_record_filename, str);
    read_population_models(msp, config);
    config_destroy(config);
    free(config);
}

static void
run_simulate(char *conf_file, long seed)
{
    int ret = -1;
    msp_t *msp = calloc(1, sizeof(msp_t));

    if (msp == NULL) {
        goto out;
    }
    msp->random_seed = seed;
    ret = msp_add_constant_population_model(msp, -1.0, 1.0);
    if (ret != 0) {
        goto out;
    }
    read_config(msp, conf_file);
    ret = msp_alloc(msp);
    if (ret != 0) {
        goto out;
    }
    ret = msp_initialise(msp);
    if (ret != 0) {
        goto out;
    }
    do {
        ret = msp_run(msp, DBL_MAX, ULONG_MAX);
        if (ret < 0) {
            goto out;
        }
        if (msp_print_state(msp) != 0) {
            goto out;
        }
    } while (ret > 0);
out:
    if (msp != NULL) {
        free(msp->coalescence_record_filename);
        msp_free(msp);
        free(msp);
    }
    if (ret != 0) {
        printf("error occured:%d:%s\n", ret, msp_strerror(ret));
    }
}


int
main(int argc, char** argv)
{
    long seed;
    if (argc != 3) {
        fatal_error("usage: %s CONFIG_FILE SEED", argv[0]);
    }
    if (parse_long(argv[2], &seed, 0, LONG_MAX) != 0) {
        fatal_error("cannot parse seed '%s'", argv[3]);
    }
    run_simulate(argv[1], seed);
    return EXIT_SUCCESS;
}
