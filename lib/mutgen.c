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
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "err.h"
#include "msprime.h"

static void
mutgen_check_state(mutgen_t *self)
{
    /* TODO some checks! */
}

void
mutgen_print_state(mutgen_t *self)
{
    printf("Mutgen state\n");
    printf("\tparamters = %s\n", self->parameters);
    printf("\trandom_seed = %d\n", (int) self->random_seed);
    printf("\tmutation_rate = %f\n", (double) self->mutation_rate);
    mutgen_check_state(self);
}

static int
mutgen_encode_parameters(mutgen_t *self)
{
    int ret = -1;
    const char *pattern = "{"
        "\"random_seed\":%lu,"
        "\"scaled_mutation_rate\":" MSP_LOSSLESS_DBL "}";
    int written;
    size_t size = 1 + (size_t) snprintf(NULL, 0, pattern,
            self->random_seed,
            self->mutation_rate);
    char *str = malloc(size);

    if (str == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    written = snprintf(str, size, pattern,
            self->random_seed,
            self->mutation_rate);
    if (written < 0) {
        ret = MSP_ERR_IO;
        goto out;
    }
    assert(written == (int) size - 1);
    self->parameters = str;
    ret = 0;
out:
    return ret;
}

int
mutgen_alloc(mutgen_t *self, tree_sequence_t *tree_sequence,
        recomb_map_t *recomb_map, double mutation_rate,
        unsigned long random_seed)
{
    int ret = MSP_ERR_NO_MEMORY;

    assert(tree_sequence != NULL);
    assert(recomb_map != NULL);
    memset(self, 0, sizeof(mutgen_t));
    self->random_seed = random_seed;
    self->mutation_rate = mutation_rate;
    self->tree_sequence = tree_sequence;
    self->recomb_map = recomb_map;
    self->rng = gsl_rng_alloc(gsl_rng_default);
    if (self->rng == NULL) {
        goto out;
    }
    gsl_rng_set(self->rng, random_seed);
    ret = mutgen_encode_parameters(self);
out:
    return ret;
}

int
mutgen_free(mutgen_t *self)
{

    if (self->rng != NULL) {
        gsl_rng_free(self->rng);
    }
    if (self->parameters != NULL) {
        free(self->parameters);
    }
    return 0;
}

#if 0
/* TODO mutations generation should be spun out into a separate class.
 * We should really just do this in Python and send the resulting
 * mutation objects here for storage. It's not an expensive operation.
 */
int
tree_sequence_generate_mutations(tree_sequence_t *self,
        recomb_map_t *recomb_map, double mutation_rate,
        unsigned long random_seed)
{
    int ret = -1;
    coalescence_record_t cr;
    uint32_t j, k, l, child;
    gsl_rng *rng = NULL;
    double *times = NULL;
    mutation_t *mutations = NULL;
    unsigned int branch_mutations;
    size_t num_mutations;
    double branch_length, distance, mu, position;
    size_t buffer_size;
    size_t block_size = 2 << 10; /* alloc in blocks of 1M */
    void *p;
    char *parameters = NULL;
    char *environment = NULL;

    buffer_size = block_size;
    num_mutations = 0;
    rng = gsl_rng_alloc(gsl_rng_default);
    if (rng == NULL) {
        goto out;
    }
    gsl_rng_set(rng, random_seed);
    times = calloc(self->num_nodes + 1, sizeof(double));
    mutations = malloc(buffer_size * sizeof(mutation_t));
    if (times == NULL || mutations == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < self->num_records; j++) {
        ret = tree_sequence_get_record(self, j, &cr, MSP_ORDER_TIME);
        if (ret != 0) {
            goto out;
        }
        times[cr.node] = cr.time;
        distance = cr.right - cr.left;
        for (k = 0; k < 2; k++) {
            child = cr.children[k];
            branch_length = cr.time - times[child];
            mu = branch_length * distance * mutation_rate;
            branch_mutations = gsl_ran_poisson(rng, mu);
            for (l = 0; l < branch_mutations; l++) {
                position = gsl_ran_flat(rng, cr.left, cr.right);
                if (num_mutations >= buffer_size) {
                    buffer_size += block_size;
                    p = realloc(mutations, buffer_size * sizeof(mutation_t));
                    if (p == NULL) {
                        ret = MSP_ERR_NO_MEMORY;
                        goto out;
                    }
                    mutations = p;
                }
                mutations[num_mutations].node = child;
                mutations[num_mutations].position = position;
                num_mutations++;
            }
        }
    }
    ret = encode_mutation_parameters(mutation_rate, random_seed, &parameters);
    if (ret != 0) {
        goto out;
    }
    ret = encode_environment(&environment);
    if (ret != 0) {
        goto out;
    }
    ret = tree_sequence_set_mutations(self, num_mutations, mutations,
            parameters, environment);
    if (ret != 0) {
        goto out;
    }
out:
    if (times != NULL) {
        free(times);
    }
    if (mutations != NULL) {
        free(mutations);
    }
    if (rng != NULL) {
        gsl_rng_free(rng);
    }
    if (parameters != NULL) {
        free(parameters);
    }
    if (environment != NULL) {
        free(environment);
    }
    return ret;
}
#endif
