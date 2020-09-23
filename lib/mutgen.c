/*
** Copyright (C) 2015-2020 University of Oxford
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
#include <float.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_minmax.h>

#include "msprime.h"

static int
cmp_site(const void *a, const void *b)
{
    const site_t *ia = (const site_t *) a;
    const site_t *ib = (const site_t *) b;
    return (ia->position > ib->position) - (ia->position < ib->position);
}

static int
cmp_mutationp(const void *a, const void *b)
{
    int out;
    /* the extra cast is to avoid a gcc bug in -Werror=cast-qual:
     * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=81631 */
    const mutation_t *ia = *(const mutation_t **) (uintptr_t) a;
    const mutation_t *ib = *(const mutation_t **) (uintptr_t) b;
    out = (ia->new || ib->new) ? 0 : (ia->id - ib->id) - (ib->id - ia->id);
    if (out == 0) {
        out = (ib->time > ia->time) - (ib->time < ia->time);
    }
    return out;
}

static void
insert_mutation(site_t *site, mutation_t *new)
{
    mutation_t *u;
    u = site->mutations;
    new->next = u;
    site->mutations = new;
    site->mutations_length++;
}

static int MSP_WARN_UNUSED
sort_mutations(site_t *site)
{
    int ret = 0;
    size_t k;
    mutation_t *m;
    size_t num_mutations = site->mutations_length;
    mutation_t **p = NULL;
    if (num_mutations > 0) {
        p = malloc(num_mutations * sizeof(*p));
        if (p == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        k = 0;
        for (m = site->mutations; m != NULL; m = m->next) {
            p[k] = m;
            k++;
        }
        assert(k == num_mutations);
        qsort(p, (size_t) num_mutations, sizeof(mutation_t *), &cmp_mutationp);
        site->mutations = p[0];
        for (k = 0; k < num_mutations; k++) {
            m = p[k];
            if (k == num_mutations - 1) {
                m->next = NULL;
            } else {
                m->next = p[k + 1];
            }
        }
    }
out:
    msp_safe_free(p);
    return ret;
}

static int MSP_WARN_UNUSED
copy_string(tsk_blkalloc_t *allocator, char *source, tsk_size_t length, char **destp,
    tsk_size_t *dest_length)
{
    int ret = 0;
    char *buff;
    *dest_length = length;
    *destp = NULL;

    if (length > 0) {
        buff = tsk_blkalloc_get(allocator, length);
        if (buff == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        memcpy(buff, source, length);
        *destp = buff;
    }
out:
    return ret;
}

/**************************
 * Mutation matrix model */

static void
mutation_matrix_print_state(mutation_model_t *self, FILE *out)
{
    size_t j, k;
    double *mutation_row;
    mutation_matrix_t params = self->params.mutation_matrix;

    fprintf(out, "mutation_matrix :: num_alleles = %d\n", (int) params.num_alleles);
    fprintf(out, "\nroot_distribution =");
    for (j = 0; j < params.num_alleles; j++) {
        fprintf(out, " '%.*s'(len=%d p=%0.4f),", (int) params.allele_length[j],
            params.alleles[j], (int) params.allele_length[j],
            params.root_distribution[j]);
    }
    fprintf(out, "\n\t------------------------------\n");
    for (j = 0; j < params.num_alleles; j++) {
        mutation_row = params.transition_matrix + j * params.num_alleles;
        fprintf(out, "\t");
        for (k = 0; k < params.num_alleles; k++) {
            fprintf(out, " %0.4f", mutation_row[k]);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
}

static bool MSP_WARN_UNUSED
valid_probabilities(size_t num_alleles, double *distribution)
{
    double prob = 0;
    double min_prob = 0;
    size_t j;
    /* Check probabilities sum to one */
    for (j = 0; j < num_alleles; j++) {
        prob += distribution[j];
        min_prob = GSL_MIN(min_prob, distribution[j]);
    }
    return doubles_almost_equal(prob, 1.0, 1e-12) && min_prob >= 0.0;
}

static int MSP_WARN_UNUSED
mutation_matrix_copy_alleles(
    mutation_matrix_t *self, char **alleles, size_t *allele_length)
{
    int ret = 0;
    size_t j;
    tsk_size_t len;

    for (j = 0; j < self->num_alleles; j++) {
        len = (tsk_size_t) allele_length[j];
        self->allele_length[j] = len;
        self->alleles[j] = malloc(len);
        if (self->alleles[j] == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        memcpy(self->alleles[j], alleles[j], len);
    }
out:
    return ret;
}

static tsk_id_t
mutation_matrix_allele_index(mutation_matrix_t *self, const char *allele, size_t length)
{
    tsk_id_t ret = -1;
    tsk_size_t j;

    for (j = 0; j < self->num_alleles; j++) {
        if (length == self->allele_length[j]
            && memcmp(allele, self->alleles[j], length) == 0) {
            ret = (tsk_id_t) j;
            break;
        }
    }
    return ret;
}

static int
mutation_matrix_choose_root_state(mutation_model_t *self, gsl_rng *rng, site_t *site)
{
    int ret = 0;
    mutation_matrix_t params = self->params.mutation_matrix;
    double u = gsl_ran_flat(rng, 0.0, 1.0);
    size_t j = probability_list_select(u, params.num_alleles, params.root_distribution);
    assert(j < params.num_alleles);
    site->ancestral_state = params.alleles[j];
    site->ancestral_state_length = params.allele_length[j];
    return ret;
}

/* Return 1 if the new derived allele matches parent allele
 * and therefore no transition occurs
 */
static int
mutation_matrix_transition(mutation_model_t *self, gsl_rng *rng,
    const char *parent_allele, tsk_size_t parent_allele_length,
    const char *MSP_UNUSED(parent_metadata),
    tsk_size_t MSP_UNUSED(parent_metadata_length), mutation_t *mutation)
{
    int ret = 0;
    mutation_matrix_t params = self->params.mutation_matrix;
    double u = gsl_ran_flat(rng, 0.0, 1.0);
    double *probs;
    tsk_id_t j, pi;

    pi = mutation_matrix_allele_index(&params, parent_allele, parent_allele_length);
    if (pi < 0) {
        /* only error if we are actually trying to mutate an unknown allele */
        ret = MSP_ERR_UNKNOWN_ALLELE;
        goto out;
    }
    probs = params.transition_matrix + (tsk_size_t) pi * params.num_alleles;
    j = (tsk_id_t) probability_list_select(u, params.num_alleles, probs);
    ret = 1;
    if (j != pi) {
        /* Only return 0 in the case where we perform an actual transition */
        ret = 0;
        mutation->derived_state = params.alleles[j];
        mutation->derived_state_length = params.allele_length[j];
    }
out:
    return ret;
}

static int
mutation_matrix_free(mutation_model_t *self)
{
    mutation_matrix_t params = self->params.mutation_matrix;
    tsk_size_t j;

    if (params.alleles != NULL) {
        for (j = 0; j < params.num_alleles; j++) {
            msp_safe_free(params.alleles[j]);
        }
    }
    msp_safe_free(params.alleles);
    msp_safe_free(params.allele_length);
    msp_safe_free(params.root_distribution);
    msp_safe_free(params.transition_matrix);
    return 0;
}

/***********************
 * SLiM mutation model */

/* Typedefs from MutationMetadataRec in slim_sim.h:
 * line 125 in v3.4, git hash b2c2b634199f35e53c4e7e513bd26c91c6d99fd9
 *  int32_t mutation_type_id_; // 4 bytes (int32_t): the id of the mutation type the
 *                             // mutation belongs to
 *  float selection_coeff_;    // 4 bytes (float): the selection coefficient
 *  int32_t subpop_index_; // 4 bytes (int32_t): the id of the subpopulation in which the
 *                         // mutation arose
 *  int32_t origin_generation_; // 4 bytes (int32_t): the generation in which the
 *                              // mutation arose
 *  int8_t nucleotide_; // 1 byte (int8_t): the nucleotide for the mutation (0='A',
 *                      // 1='C', 2='G', 3='T'), or -1
 *
 * Note that these are defined there as a __packed__ struct, like
 * typedef struct __attribute__((__packed__))  but this is not available
 * in Windows compilers, so we're just copying the info in directly */

#define SLIM_MUTATION_METADATA_SIZE 17 // = 4 + 4 + 4 + 4 + 1

static void
copy_slim_mutation_metadata(slim_mutator_t *params, char *dest)
{
    size_t n;
    int32_t *mutation_type_id;
    float *selection_coeff;
    int32_t *subpop_index;
    int32_t *origin_generation;
    int8_t *nucleotide;

    n = 0;
    mutation_type_id = (int32_t *) (dest + n);
    *mutation_type_id = params->mutation_type_id;
    n += sizeof(int32_t);
    selection_coeff = (float *) (dest + n);
    *selection_coeff = 0.0;
    n += sizeof(float);
    subpop_index = (int32_t *) (dest + n);
    *subpop_index = TSK_NULL;
    n += sizeof(int32_t);
    // TODO: remove this when switch to mutation time
    origin_generation = (int32_t *) (dest + n);
    *origin_generation = 0.0;
    n += sizeof(int32_t);
    nucleotide = (int8_t *) (dest + n);
    *nucleotide = -1;
}

static void
slim_mutator_print_state(mutation_model_t *self, FILE *out)
{
    slim_mutator_t params = self->params.slim_mutator;
    fprintf(out, "SLiM mutation model :: mutation type ID = %d\n",
        (int) params.mutation_type_id);
    fprintf(out, "                       next mutation ID = %d\n",
        (int) params.next_mutation_id);
}

static int MSP_WARN_UNUSED
slim_mutator_check_validity(slim_mutator_t *self)
{
    int ret = 0;

    if (self->next_mutation_id < 0) {
        ret = MSP_ERR_BAD_SLIM_PARAMETERS;
        goto out;
    }

    if (self->mutation_type_id < 0) {
        ret = MSP_ERR_BAD_SLIM_PARAMETERS;
        goto out;
    }

out:
    return ret;
}

static int
slim_mutator_choose_root_state(
    mutation_model_t *MSP_UNUSED(self), gsl_rng *MSP_UNUSED(rng), site_t *site)
{
    site->ancestral_state = NULL;
    site->ancestral_state_length = 0;
    return 0;
}

static int
slim_mutator_transition(mutation_model_t *self, gsl_rng *MSP_UNUSED(rng),
    const char *parent_allele, tsk_size_t parent_allele_length,
    const char *parent_metadata, tsk_size_t parent_metadata_length, mutation_t *mutation)
{
    int ret = 0;
    slim_mutator_t *params = &self->params.slim_mutator;
    char *buff = NULL;
    int len;
    /* The maximum number of digits for a signed 64 bit integer (including
     * the leading "-") */
    const size_t max_digits = 20;
    /* We allow for a possible comma to separate the previous element
     * in the list, as well as a NULL byte added by snprintf. We don't bother
     * trying to alloc the exact number of bytes needed. */
    const size_t alloc_size = parent_allele_length + max_digits + 2;
    const char *sep = parent_allele_length == 0 ? "" : ",";

    /* Append to derived_state */
    buff = tsk_blkalloc_get(&params->allocator, alloc_size);
    if (buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    len = snprintf(buff, alloc_size, "%.*s%s%" PRId64, parent_allele_length,
        parent_allele, sep, params->next_mutation_id);
    if (len < 0) {
        /* Technically this can happen. Returning TSK_ERR_IO should result
         * in Python code checking errno. */
        ret = TSK_ERR_IO;
        goto out;
    }
    assert(len < (int) alloc_size);
    if (params->next_mutation_id == INT64_MAX) {
        ret = MSP_ERR_MUTATION_ID_OVERFLOW;
        goto out;
    }
    params->next_mutation_id++;
    mutation->derived_state = buff;
    mutation->derived_state_length = (tsk_size_t) len;

    /* Append to metadata */
    buff = tsk_blkalloc_get(
        &params->allocator, parent_metadata_length + SLIM_MUTATION_METADATA_SIZE);
    if (buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memcpy(buff, parent_metadata, parent_metadata_length);
    copy_slim_mutation_metadata(params, buff + parent_metadata_length);

    mutation->metadata = buff;
    mutation->metadata_length
        = (tsk_size_t)(parent_metadata_length + SLIM_MUTATION_METADATA_SIZE);
out:
    return ret;
}

static int
slim_mutator_free(mutation_model_t *self)
{
    slim_mutator_t params = self->params.slim_mutator;
    tsk_blkalloc_free(&params.allocator);
    return 0;
}

/***********************
 * Infinite alleles mutation model */

static void
infinite_alleles_print_state(mutation_model_t *self, FILE *out)
{
    infinite_alleles_t params = self->params.infinite_alleles;
    fprintf(out, "infinite alleles:\n");
    fprintf(out, "\tstart_allele:%" PRIu64 "\n", params.start_allele);
    fprintf(out, "\tnext_allele:%" PRIu64 "\n", params.next_allele);
}

/* The maximum number of digits for an unsigned 64 bit integer is 20
 * and one byte for the NULL terminator. */
#define MAX_UINT_BUFF_SIZE 21

static int
infinite_alleles_make_allele(
    mutation_model_t *self, char **dest, tsk_size_t *dest_length)
{
    int ret = 0;
    infinite_alleles_t *params = &self->params.infinite_alleles;
    char tmp_buff[MAX_UINT_BUFF_SIZE];
    char *buff = NULL;
    int num_digits;
    tsk_size_t len;

    num_digits = snprintf(tmp_buff, MAX_UINT_BUFF_SIZE, "%" PRIu64, params->next_allele);
    if (num_digits < 0) {
        /* Technically this can happen. Returning TSK_ERR_IO should result
         * in Python code checking errno. */
        ret = TSK_ERR_IO;
        goto out;
    }
    len = (tsk_size_t) num_digits;
    assert(len < MAX_UINT_BUFF_SIZE);

    buff = tsk_blkalloc_get(&params->allocator, len);
    if (buff == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memcpy(buff, tmp_buff, len);
    params->next_allele++;
    *dest = buff;
    *dest_length = len;
out:
    return ret;
}

static int
infinite_alleles_choose_root_state(
    mutation_model_t *self, gsl_rng *MSP_UNUSED(rng), site_t *site)
{
    return infinite_alleles_make_allele(
        self, &site->ancestral_state, &site->ancestral_state_length);
}

static int
infinite_alleles_transition(mutation_model_t *self, gsl_rng *MSP_UNUSED(rng),
    const char *MSP_UNUSED(parent_allele), tsk_size_t MSP_UNUSED(parent_allele_length),
    const char *MSP_UNUSED(parent_metadata),
    tsk_size_t MSP_UNUSED(parent_metadata_length), mutation_t *mutation)
{
    return infinite_alleles_make_allele(
        self, &mutation->derived_state, &mutation->derived_state_length);
}

static int
infinite_alleles_free(mutation_model_t *self)
{
    infinite_alleles_t params = self->params.infinite_alleles;
    tsk_blkalloc_free(&params.allocator);
    return 0;
}

/*****************************************
 * Mutation model public API
 *
 * This consists of a factory function for each type of model supported
 * (e.g., matrix_mutation_model_alloc), mutation_matrix_free,
 * mutation_model_choose_root_state and mutation_model_transition.
 * These delegate the actual implementations to function pointers that
 * are implemented by the actual models.
 */

int MSP_WARN_UNUSED
matrix_mutation_model_alloc(mutation_model_t *self, size_t num_alleles, char **alleles,
    size_t *allele_lengths, double *root_distribution, double *transition_matrix)
{
    int ret = 0;
    mutation_matrix_t *params = &self->params.mutation_matrix;
    size_t i;

    memset(self, 0, sizeof(*self));
    if (num_alleles < 2) {
        ret = MSP_ERR_INSUFFICIENT_ALLELES;
        goto out;
    }
    if (!valid_probabilities(num_alleles, root_distribution)) {
        ret = MSP_ERR_BAD_ROOT_PROBABILITIES;
        goto out;
    }
    for (i = 0; i < num_alleles; ++i) {
        double *probs = transition_matrix + i * num_alleles;
        if (!valid_probabilities(num_alleles, probs)) {
            ret = MSP_ERR_BAD_TRANSITION_MATRIX;
            goto out;
        }
    }
    params->num_alleles = num_alleles;
    params->alleles = calloc(num_alleles, sizeof(*params->alleles));
    params->allele_length = calloc(num_alleles, sizeof(*params->allele_length));
    params->root_distribution = malloc(num_alleles * sizeof(*params->root_distribution));
    params->transition_matrix
        = malloc(num_alleles * num_alleles * sizeof(*params->transition_matrix));
    if (params->alleles == NULL || params->allele_length == NULL
        || params->root_distribution == NULL || params->transition_matrix == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memcpy(params->root_distribution, root_distribution,
        num_alleles * sizeof(*root_distribution));
    memcpy(params->transition_matrix, transition_matrix,
        num_alleles * num_alleles * sizeof(*transition_matrix));
    ret = mutation_matrix_copy_alleles(params, alleles, allele_lengths);
    if (ret != 0) {
        goto out;
    }
    self->choose_root_state = &mutation_matrix_choose_root_state;
    self->transition = &mutation_matrix_transition;
    self->print_state = &mutation_matrix_print_state;
    self->free = &mutation_matrix_free;
out:
    return ret;
}

int MSP_WARN_UNUSED
slim_mutation_model_alloc(mutation_model_t *self, int32_t mutation_type_id,
    int64_t next_mutation_id, size_t block_size)
{
    int ret = 0;
    slim_mutator_t *params = &self->params.slim_mutator;

    memset(self, 0, sizeof(*self));

    self->choose_root_state = &slim_mutator_choose_root_state;
    self->transition = &slim_mutator_transition;
    self->print_state = &slim_mutator_print_state;
    self->free = &slim_mutator_free;
    if (block_size == 0) {
        /* 8K is a good default, but we need to have the
         * block_size argument here because the size of the allocations
         * we need to make are in principle unbounded since SLiM
         * copies the entire parent state for every mutation as it
         * goes down along the tree.
         */
        block_size = 8192;
    }
    ret = tsk_blkalloc_init(&params->allocator, block_size);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
    params->mutation_type_id = mutation_type_id;
    params->next_mutation_id = next_mutation_id;

    ret = slim_mutator_check_validity(params);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int MSP_WARN_UNUSED
infinite_alleles_mutation_model_alloc(
    mutation_model_t *self, uint64_t start_allele, tsk_flags_t MSP_UNUSED(options))
{
    int ret = 0;
    infinite_alleles_t *params = &self->params.infinite_alleles;

    memset(self, 0, sizeof(*self));

    self->choose_root_state = &infinite_alleles_choose_root_state;
    self->transition = &infinite_alleles_transition;
    self->print_state = &infinite_alleles_print_state;
    self->free = &infinite_alleles_free;
    ret = tsk_blkalloc_init(&params->allocator, 8192);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
    params->start_allele = start_allele;
    params->next_allele = start_allele;
out:
    return ret;
}

int
mutation_model_free(mutation_model_t *self)
{
    if (self->free != NULL) {
        self->free(self);
    }
    return 0;
}

static int MSP_WARN_UNUSED
mutation_model_choose_root_state(mutation_model_t *self, gsl_rng *rng, site_t *site)
{
    return self->choose_root_state(self, rng, site);
}

static int MSP_WARN_UNUSED
mutation_model_transition(mutation_model_t *self, gsl_rng *rng,
    const char *parent_allele, tsk_size_t parent_allele_length,
    const char *parent_metadata, tsk_size_t parent_metadata_length, mutation_t *mutation)
{
    return self->transition(self, rng, parent_allele, parent_allele_length,
        parent_metadata, parent_metadata_length, mutation);
}

void
mutation_model_print_state(mutation_model_t *self, FILE *out)
{
    self->print_state(self, out);
}

static const char *binary_alleles[] = { "0", "1" };
static double binary_model_root_distribution[] = { 1.0, 0.0 };
static double binary_model_transition_matrix[] = { 0.0, 1.0, 1.0, 0.0 };

#define ONE_THIRD (1.0 / 3.0)
static const char *acgt_alleles[] = { "A", "C", "G", "T" };
static double jukes_cantor_model_root_distribution[] = { 0.25, 0.25, 0.25, 0.25 };
static double jukes_cantor_model_transition_matrix[]
    = { 0.0, ONE_THIRD, ONE_THIRD, ONE_THIRD, ONE_THIRD, 0.0, ONE_THIRD, ONE_THIRD,
          ONE_THIRD, ONE_THIRD, 0.0, ONE_THIRD, ONE_THIRD, ONE_THIRD, ONE_THIRD, 0.0 };

/* Populates this mutation model instance with simple mutation model.
 *
 * 0 -> Symmetric 0/1 mutations
 * 1 -> Jukes Cantor ACGT mutations.
 *
 * This is only provided for testing purposes to make it easy to
 * obtain a populated model.
 */
int
matrix_mutation_model_factory(mutation_model_t *self, int model)
{
    int ret = MSP_ERR_GENERIC;
    size_t lengths[] = { 1, 1, 1, 1 };

    /* We need to cast to uintptr_t * first to work around the annoying pedantry
     * about discarding const qualifiers. */
    if (model == 0) {
        ret = matrix_mutation_model_alloc(self, 2,
            (char **) (uintptr_t *) binary_alleles, lengths,
            binary_model_root_distribution, binary_model_transition_matrix);
    } else if (model == 1) {
        ret = matrix_mutation_model_alloc(self, 4, (char **) (uintptr_t *) acgt_alleles,
            lengths, jukes_cantor_model_root_distribution,
            jukes_cantor_model_transition_matrix);
    }
    return ret;
}

/***********************
 * Mutation generator */

static void
mutgen_check_state(mutgen_t *self)
{
    size_t j;
    avl_node_t *a;
    site_t *s;
    mutation_t *m;

    for (a = self->sites.head; a != NULL; a = a->next) {
        s = (site_t *) a->item;
        m = s->mutations;
        for (j = 0; j < s->mutations_length; j++) {
            assert(m != NULL);
            assert(m->id >= -1);
            assert(m->node >= 0);
            if (j == s->mutations_length - 1) {
                assert(m->next == NULL);
            }
            m = m->next;
        }
        assert(m == NULL);
    }
}

void
mutgen_print_state(mutgen_t *self, FILE *out)
{
    avl_node_t *a;
    site_t *s;
    mutation_t *m;
    tsk_id_t parent_id;

    fprintf(out, "Mutgen state\n");
    fprintf(out, "\trate_map:\n");
    rate_map_print_state(&self->rate_map, out);
    fprintf(out, "\tstart_time = %f\n", self->start_time);
    fprintf(out, "\tend_time = %f\n", self->end_time);
    fprintf(out, "\tmodel:\n");
    mutation_model_print_state(self->model, out);
    tsk_blkalloc_print_state(&self->allocator, out);

    for (a = self->sites.head; a != NULL; a = a->next) {
        s = (site_t *) a->item;
        fprintf(out, "site:\t%f\t'%.*s'\t'%.*s'\t(%d)\t%d\n", s->position,
            (int) s->ancestral_state_length, s->ancestral_state,
            (int) s->metadata_length, s->metadata, s->new, (int) s->mutations_length);
        for (m = s->mutations; m != NULL; m = m->next) {
            parent_id = m->parent == NULL ? TSK_NULL : m->parent->id;
            fprintf(out, "\tmut:\t(%d)\t%f\t%d\t%d\t'%.*s'\t'%.*s'\t(%d)\t%d\n", m->id,
                m->time, m->node, parent_id, (int) m->derived_state_length,
                m->derived_state, (int) m->metadata_length, m->metadata, m->new,
                m->keep);
        }
    }
    mutgen_check_state(self);
}

int MSP_WARN_UNUSED
mutgen_alloc(mutgen_t *self, gsl_rng *rng, tsk_table_collection_t *tables,
    mutation_model_t *model, size_t block_size)
{
    int ret = 0;

    memset(self, 0, sizeof(mutgen_t));
    if (rng == NULL || tables == NULL || model == NULL) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->rng = rng;
    self->tables = tables;
    self->model = model;
    self->start_time = -DBL_MAX;
    self->end_time = DBL_MAX;
    self->block_size = block_size;

    if (tables->sequence_length <= 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    avl_init_tree(&self->sites, cmp_site, NULL);
    ret = mutgen_set_rate(self, 0);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int
mutgen_free(mutgen_t *self)
{
    tsk_blkalloc_free(&self->allocator);
    rate_map_free(&self->rate_map);
    return 0;
}

int MSP_WARN_UNUSED
mutgen_set_time_interval(mutgen_t *self, double start_time, double end_time)
{
    int ret = 0;

    if (end_time < start_time) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->start_time = start_time;
    self->end_time = end_time;
out:
    return ret;
}

int
mutgen_set_rate_map(mutgen_t *self, size_t size, double *position, double *rate)
{
    int ret = 0;

    rate_map_free(&self->rate_map);

    ret = rate_map_alloc(&self->rate_map, size, position, rate);
    if (ret != 0) {
        goto out;
    }
    if (rate_map_get_sequence_length(&self->rate_map) != self->tables->sequence_length) {
        ret = MSP_ERR_INCOMPATIBLE_MUTATION_MAP;
        goto out;
    }
out:
    return ret;
}

/* Short-cut for mutgen_set_recombination_map can be used in testing. */
int
mutgen_set_rate(mutgen_t *self, double rate)
{
    double position[] = { 0, self->tables->sequence_length };
    return mutgen_set_rate_map(self, 1, position, &rate);
}

static int MSP_WARN_UNUSED
mutgen_init_allocator(mutgen_t *self)
{
    int ret = -1;

    tsk_blkalloc_free(&self->allocator);
    if (self->block_size == 0) {
        /* Default */
        self->block_size = 8192;
    }
    /* This is the effective minimum */
    self->block_size = GSL_MAX(self->block_size, 128);
    /* Need to make sure we have enough space to store sites and mutations. We
     * allocate ancestral and derived states, as well as a list of mutations
     * for each site. This ensures that we can always allocate the required amount.
     * We need to add one because the assert trips when the block size is equal
     * to chunk size (probably wrongly).
     */
    self->block_size
        = GSL_MAX(self->block_size, 1 + self->tables->sites.ancestral_state_length);
    self->block_size
        = GSL_MAX(self->block_size, 1 + self->tables->sites.metadata_length);
    self->block_size
        = GSL_MAX(self->block_size, 1 + self->tables->mutations.derived_state_length);
    self->block_size
        = GSL_MAX(self->block_size, 1 + self->tables->mutations.metadata_length);
    self->block_size = GSL_MAX(
        self->block_size, (1 + self->tables->mutations.num_rows) * sizeof(mutation_t));
    ret = tsk_blkalloc_init(&self->allocator, self->block_size);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_add_new_site(mutgen_t *self, double position, site_t **new_site)
{
    int ret = 0;
    avl_node_t *avl_node;
    site_t *site;

    avl_node = tsk_blkalloc_get(&self->allocator, sizeof(*avl_node));
    site = tsk_blkalloc_get(&self->allocator, sizeof(*site));
    if (site == NULL || avl_node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memset(site, 0, sizeof(*site));
    site->position = position;
    site->new = true;

    avl_init_node(avl_node, site);
    avl_node = avl_insert_node(&self->sites, avl_node);
    if (avl_node == NULL) {
        ret = MSP_ERR_DUPLICATE_SITE_POSITION;
        goto out;
    }
    *new_site = site;
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_add_existing_site(mutgen_t *self, double position, char *ancestral_state,
    tsk_size_t ancestral_state_length, char *metadata, tsk_size_t metadata_length,
    site_t **new_site)
{
    int ret = 0;
    site_t *site;

    ret = mutgen_add_new_site(self, position, &site);
    if (ret != 0) {
        goto out;
    }
    site->new = false;

    /* We need to copy the ancestral state and metadata  */
    ret = copy_string(&self->allocator, ancestral_state, ancestral_state_length,
        &site->ancestral_state, &site->ancestral_state_length);
    if (ret != 0) {
        goto out;
    }
    ret = copy_string(&self->allocator, metadata, metadata_length, &site->metadata,
        &site->metadata_length);
    if (ret != 0) {
        goto out;
    }
    *new_site = site;
out:
    return ret;
}

static int
mutgen_add_mutation(
    mutgen_t *self, site_t *site, tsk_id_t node, double time, mutation_t **new_mutation)
{
    int ret = 0;

    mutation_t *mutation = tsk_blkalloc_get(&self->allocator, sizeof(*mutation));
    if (mutation == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    memset(mutation, 0, sizeof(*mutation));
    mutation->node = node;
    mutation->time = time;
    mutation->parent = NULL;
    mutation->new = true;
    insert_mutation(site, mutation);
    *new_mutation = mutation;
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_add_new_mutation(mutgen_t *self, site_t *site, tsk_id_t node, double time)
{
    mutation_t *not_used;
    return mutgen_add_mutation(self, site, node, time, &not_used);
}

static int MSP_WARN_UNUSED
mutgen_add_existing_mutation(mutgen_t *self, site_t *site, tsk_id_t id, tsk_id_t node,
    double time, char *derived_state, tsk_size_t derived_state_length, char *metadata,
    tsk_size_t metadata_length)
{
    int ret = 0;
    mutation_t *mutation;

    ret = mutgen_add_mutation(self, site, node, time, &mutation);
    if (ret != 0) {
        goto out;
    }
    mutation->id = id;
    mutation->new = false;

    /* Need to copy the derived state and metadata */
    ret = copy_string(&self->allocator, derived_state, derived_state_length,
        &mutation->derived_state, &mutation->derived_state_length);
    if (ret != 0) {
        goto out;
    }
    ret = copy_string(&self->allocator, metadata, metadata_length, &mutation->metadata,
        &mutation->metadata_length);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_initialise_sites(mutgen_t *self, bool check_kept_times)
{
    int ret = 0;
    tsk_site_table_t *sites = &self->tables->sites;
    tsk_mutation_table_t *mutations = &self->tables->mutations;
    tsk_id_t site_id;
    site_t *site;
    double time;
    const double end_time = self->end_time;
    char *state, *metadata;
    tsk_size_t j, length, metadata_length;

    j = 0;
    for (site_id = 0; site_id < (tsk_id_t) sites->num_rows; site_id++) {
        state = sites->ancestral_state + sites->ancestral_state_offset[site_id];
        length = sites->ancestral_state_offset[site_id + 1]
                 - sites->ancestral_state_offset[site_id];
        metadata = sites->metadata + sites->metadata_offset[site_id];
        metadata_length
            = sites->metadata_offset[site_id + 1] - sites->metadata_offset[site_id];

        ret = mutgen_add_existing_site(self, sites->position[site_id], state, length,
            metadata, metadata_length, &site);
        if (ret != 0) {
            goto out;
        }

        while (j < mutations->num_rows && mutations->site[j] == site_id) {
            assert(j < mutations->num_rows);
            time = mutations->time[j];
            if (tsk_is_unknown_time(time)) {
                ret = MSP_ERR_UNKNOWN_TIME_NOT_SUPPORTED;
                goto out;
            }
            // check if any kept mutations are younger than
            // the time period where new mutations can be added
            if (check_kept_times && time < end_time) {
                ret = MSP_ERR_MUTATION_GENERATION_OUT_OF_ORDER;
                goto out;
            }
            state = mutations->derived_state + mutations->derived_state_offset[j];
            length = mutations->derived_state_offset[j + 1]
                     - mutations->derived_state_offset[j];
            metadata = mutations->metadata + mutations->metadata_offset[j];
            metadata_length
                = mutations->metadata_offset[j + 1] - mutations->metadata_offset[j];
            ret = mutgen_add_existing_mutation(self, site, (int) j, mutations->node[j],
                time, state, length, metadata, metadata_length);
            if (ret != 0) {
                goto out;
            }
            j++;
        }
    }

out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_populate_tables(mutgen_t *self)
{
    int ret = 0;
    tsk_site_table_t *sites = &self->tables->sites;
    tsk_mutation_table_t *mutations = &self->tables->mutations;
    tsk_id_t site_id, mutation_id, parent_id;
    avl_node_t *a;
    site_t *site;
    mutation_t *m;
    size_t num_mutations;

    site_id = 0;
    for (a = self->sites.head; a != NULL; a = a->next) {
        site = (site_t *) a->item;
        num_mutations = 0;
        for (m = site->mutations; m != NULL; m = m->next) {
            if (m->keep) {
                if (m->parent == NULL) {
                    parent_id = TSK_NULL;
                } else {
                    parent_id = m->parent->id;
                    assert(parent_id != TSK_NULL);
                }
                mutation_id = tsk_mutation_table_add_row(mutations, site_id, m->node,
                    parent_id, m->time, m->derived_state, m->derived_state_length,
                    m->metadata, m->metadata_length);
                if (mutation_id < 0) {
                    ret = msp_set_tsk_error(mutation_id);
                    goto out;
                }
                assert(mutation_id > parent_id);
                m->id = mutation_id;
                num_mutations++;
            }
        }
        /* Omit any new sites that have no mutations */
        if ((!site->new) || num_mutations > 0) {
            ret = tsk_site_table_add_row(sites, site->position, site->ancestral_state,
                site->ancestral_state_length, site->metadata, site->metadata_length);
            if (ret < 0) {
                goto out;
            }
            site_id++;
        }
    }
    ret = 0;
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_place_mutations(mutgen_t *self, bool discrete_sites)
{
    /* The mutation model for discrete sites is that there is
     * a unit of "mutation mass" on each integer, so that
     * in a segment [left, right) there can be mutations at the
     * integers {ceil(left), ceil(left) + 1, ..., ceil(right) - 1},
     * and the total mutation "length" is ceil(right) - ceil(left). */
    int ret = 0;
    const double *map_position = self->rate_map.position;
    const double *map_rate = self->rate_map.rate;
    size_t branch_mutations, map_index;
    size_t j, k;
    const tsk_node_table_t nodes = self->tables->nodes;
    const tsk_edge_table_t edges = self->tables->edges;
    const double start_time = self->start_time;
    const double end_time = self->end_time;
    double left, right, site_left, site_right, edge_right;
    double time, mu, position;
    double branch_start, branch_end, branch_length;
    tsk_id_t parent, child;
    avl_node_t *avl_node;
    site_t *site;
    site_t search;

    for (j = 0; j < edges.num_rows; j++) {
        left = edges.left[j];
        edge_right = edges.right[j];
        parent = edges.parent[j];
        child = edges.child[j];
        assert(child >= 0 && child < (tsk_id_t) nodes.num_rows);
        branch_start = GSL_MAX(start_time, nodes.time[child]);
        branch_end = GSL_MIN(end_time, nodes.time[parent]);
        branch_length = branch_end - branch_start;

        map_index = rate_map_get_index(&self->rate_map, left);
        right = 0;
        while (right != edge_right) {
            right = GSL_MIN(edge_right, map_position[map_index + 1]);
            site_left = discrete_sites ? ceil(left) : left;
            site_right = discrete_sites ? ceil(right) : right;
            mu = branch_length * (site_right - site_left) * map_rate[map_index];
            branch_mutations = gsl_ran_poisson(self->rng, mu);
            for (k = 0; k < branch_mutations; k++) {
                /* Rejection sample positions until we get one we haven't seen before,
                 * unless we are doing discrete sites. Note that in principle this
                 * could lead to an infinite loop here, but in practise we'd need to
                 * use up all of the doubles before it could happen and so we'd
                 * certainly run out of memory first. */
                do {
                    position = gsl_ran_flat(self->rng, site_left, site_right);
                    if (discrete_sites) {
                        position = floor(position);
                    }
                    search.position = position;
                    avl_node = avl_search(&self->sites, &search);
                } while (avl_node != NULL && !discrete_sites);

                time = gsl_ran_flat(self->rng, branch_start, branch_end);
                assert(site_left <= position && position < site_right);
                assert(branch_start <= time && time < branch_end);
                if (avl_node != NULL) {
                    site = (site_t *) avl_node->item;
                } else {
                    ret = mutgen_add_new_site(self, position, &site);
                    if (ret != 0) {
                        goto out;
                    }
                }
                ret = mutgen_add_new_mutation(self, site, child, time);
                if (ret != 0) {
                    goto out;
                }
            }
            map_index++;
        }
    }
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_choose_alleles(mutgen_t *self, tsk_id_t *parent, mutation_t **bottom_mutation,
    tsk_size_t num_nodes, site_t *site)
{
    int ret = 0;
    const char *pa, *pm;
    tsk_size_t palen, pmlen;
    mutation_t *mut, *parent_mut;
    tsk_id_t u;

    ret = sort_mutations(site);
    if (ret != 0) {
        goto out;
    }
    if (site->new) {
        assert(site->ancestral_state == NULL);
        ret = mutation_model_choose_root_state(self->model, self->rng, site);
        if (ret != 0) {
            goto out;
        }
    }

    /* Create a mapping from mutations to nodes in bottom_mutation. If we see
     * more than one mutation at a node, the previously seen one must be the
     * parent of the current one since we assume they are in order. */
    for (mut = site->mutations; mut != NULL; mut = mut->next) {
        u = mut->node;
        assert((tsk_size_t) u < num_nodes);
        while (u != TSK_NULL && bottom_mutation[u] == NULL) {
            u = parent[u];
        }
        if (u == TSK_NULL) {
            pa = site->ancestral_state;
            palen = site->ancestral_state_length;
            pm = site->metadata;
            pmlen = site->metadata_length;
            assert(mut->parent == NULL);
        } else {
            parent_mut = bottom_mutation[u];
            mut->parent = parent_mut;
            assert(mut->time <= parent_mut->time);
            if (mut->new) {
                pa = parent_mut->derived_state;
                palen = parent_mut->derived_state_length;
                pm = parent_mut->metadata;
                pmlen = parent_mut->metadata_length;
            }
        }
        mut->keep = true;
        if (mut->new) {
            assert(mut->derived_state == NULL);
            ret = mutation_model_transition(
                self->model, self->rng, pa, palen, pm, pmlen, mut);
            if (ret < 0) {
                goto out;
            }
            /* A non-zero return value here indicates that we transitioned to the
             * same allele, and so the discard the mutation. */
            if (ret != 0) {
                mut->keep = false;
            }
        }
        if (mut->keep) {
            bottom_mutation[mut->node] = mut;
        }
    }
    /* Reset the mapping for the next site */
    for (mut = site->mutations; mut != NULL; mut = mut->next) {
        bottom_mutation[mut->node] = NULL;
    }
    ret = 0;
out:
    return ret;
}

static int MSP_WARN_UNUSED
mutgen_apply_mutations(mutgen_t *self)
{
    int ret = 0;
    const tsk_id_t *I, *O;
    const tsk_edge_table_t edges = self->tables->edges;
    const tsk_node_table_t nodes = self->tables->nodes;
    const tsk_id_t M = (tsk_id_t) edges.num_rows;
    tsk_id_t tj, tk;
    tsk_id_t *parent = NULL;
    mutation_t **bottom_mutation = NULL;
    double left, right;
    const double sequence_length = self->tables->sequence_length;
    avl_node_t *avl_node;
    site_t *site;

    parent = malloc(nodes.num_rows * sizeof(*parent));
    bottom_mutation = malloc(nodes.num_rows * sizeof(*bottom_mutation));
    if (parent == NULL || bottom_mutation == NULL) {
        ret = TSK_ERR_NO_MEMORY;
        goto out;
    }
    memset(parent, 0xff, nodes.num_rows * sizeof(*parent));
    memset(bottom_mutation, 0, nodes.num_rows * sizeof(*bottom_mutation));

    if (!tsk_table_collection_has_index(self->tables, 0)) {
        ret = tsk_table_collection_build_index(self->tables, 0);
        if (ret != 0) {
            goto out;
        }
    }

    I = self->tables->indexes.edge_insertion_order;
    O = self->tables->indexes.edge_removal_order;
    tj = 0;
    tk = 0;
    left = 0;
    avl_node = self->sites.head;
    while (tj < M || left < sequence_length) {
        while (tk < M && edges.right[O[tk]] == left) {
            parent[edges.child[O[tk]]] = TSK_NULL;
            tk++;
        }
        while (tj < M && edges.left[I[tj]] == left) {
            parent[edges.child[I[tj]]] = edges.parent[I[tj]];
            tj++;
        }
        right = sequence_length;
        if (tj < M) {
            right = TSK_MIN(right, edges.left[I[tj]]);
        }
        if (tk < M) {
            right = TSK_MIN(right, edges.right[O[tk]]);
        }

        /* Tree is now ready. We look at each site on this tree in turn */
        while (avl_node != NULL) {
            site = (site_t *) avl_node->item;
            if (site->position >= right) {
                break;
            }
            ret = mutgen_choose_alleles(
                self, parent, bottom_mutation, nodes.num_rows, site);
            if (ret != 0) {
                goto out;
            }
            avl_node = avl_node->next;
        }
        /* Move on to the next tree */
        left = right;
    }

out:
    msp_safe_free(parent);
    msp_safe_free(bottom_mutation);
    return ret;
}

int MSP_WARN_UNUSED
mutgen_generate(mutgen_t *self, int flags)
{
    int ret = 0;
    bool discrete_sites = flags & MSP_DISCRETE_SITES;
    bool kept_mutations_before_end_time = flags & MSP_KEPT_MUTATIONS_BEFORE_END_TIME;

    avl_clear_tree(&self->sites);

    ret = mutgen_init_allocator(self);
    if (ret != 0) {
        goto out;
    }
    ret = tsk_table_collection_check_integrity(self->tables, 0);
    if (ret != 0) {
        ret = msp_set_tsk_error(ret);
        goto out;
    }
    if (flags & MSP_KEEP_SITES) {
        ret = mutgen_initialise_sites(
            self, discrete_sites && !kept_mutations_before_end_time);
        if (ret != 0) {
            goto out;
        }
    }

    ret = tsk_site_table_clear(&self->tables->sites);
    if (ret != 0) {
        goto out;
    }
    ret = tsk_mutation_table_clear(&self->tables->mutations);
    if (ret != 0) {
        goto out;
    }
    ret = mutgen_place_mutations(self, discrete_sites);
    if (ret != 0) {
        goto out;
    }
    ret = mutgen_apply_mutations(self);
    if (ret != 0) {
        goto out;
    }
    ret = mutgen_populate_tables(self);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}
