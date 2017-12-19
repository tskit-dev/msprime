/*
** Copyright (C) 2016-2017 University of Oxford
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
#include "object_heap.h"
#include "msprime.h"

void
vargen_print_state(vargen_t *self, FILE *out)
{
    fprintf(out, "vargen state\n");
    fprintf(out, "tree_site_index = %d\n", (int) self->tree_site_index);
}

static int
vargen_next_tree(vargen_t *self)
{
    int ret = 0;

    ret = sparse_tree_next(&self->tree);
    if (ret == 0) {
        self->finished = 1;
    } else if (ret < 0) {
        goto out;
    }
    self->tree_site_index = 0;
out:
    return ret;
}

int
vargen_alloc(vargen_t *self, tree_sequence_t *tree_sequence, int flags)
{
    int ret = MSP_ERR_NO_MEMORY;
    size_t max_alleles;

    assert(tree_sequence != NULL);
    memset(self, 0, sizeof(vargen_t));

    /* For now, the logic only supports infinite sites binary mutations. We need to
     * think about how to structure this API to support the general case (lots of
     * mutations happening along the tree) without making it too inefficient and
     * breaking too much code.
     */
    if (tree_sequence_get_alphabet(tree_sequence) != MSP_ALPHABET_BINARY) {
        ret = MSP_ERR_NONBINARY_MUTATIONS_UNSUPPORTED;
        goto out;
    }
    self->num_samples = tree_sequence_get_num_samples(tree_sequence);
    self->num_sites = tree_sequence_get_num_sites(tree_sequence);
    self->tree_sequence = tree_sequence;
    self->flags = flags;
    max_alleles = 1 + tree_sequence_get_max_site_mutations(tree_sequence);
    self->variant.genotypes = malloc(self->num_samples * sizeof(*self->variant.genotypes));
    self->variant.alleles = malloc(max_alleles * sizeof(*self->variant.alleles));
    self->variant.allele_lengths = malloc(max_alleles
            * sizeof(*self->variant.allele_lengths));
    if (self->variant.genotypes == NULL || self->variant.alleles == NULL
            || self->variant.allele_lengths == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = tree_sequence_get_sample_index_map(tree_sequence, &self->sample_index_map);
    if (ret != 0) {
        goto out;
    }
    ret = sparse_tree_alloc(&self->tree, tree_sequence, MSP_SAMPLE_LISTS);
    if (ret != 0) {
        goto out;
    }
    self->finished = 0;
    self->tree_site_index = 0;
    ret = sparse_tree_first(&self->tree);
    if (ret < 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

int
vargen_free(vargen_t *self)
{
    sparse_tree_free(&self->tree);
    msp_safe_free(self->variant.genotypes);
    msp_safe_free(self->variant.alleles);
    msp_safe_free(self->variant.allele_lengths);
    return 0;
}

static int
vargen_apply_tree_site(vargen_t *self, site_t *site, char *genotypes, char state_offset)
{
    int ret = 0;
    node_list_t *w, *tail;
    node_id_t sample_index;
    bool not_done;
    table_size_t j;
    char derived;
    char ancestral = (char) (site->ancestral_state[0] - state_offset);

    memset(genotypes, ancestral, self->num_samples);
    for (j = 0; j < site->mutations_length; j++) {
        derived = (char) (site->mutations[j].derived_state[0] - state_offset);
        ret = sparse_tree_get_sample_list(&self->tree, site->mutations[j].node, &w, &tail);
        if (ret != 0) {
            goto out;
        }
        if (w != NULL) {
            not_done = true;
            while (not_done) {
                assert(w != NULL);
                sample_index = self->sample_index_map[w->node];
                assert(sample_index >= 0);
                if (genotypes[sample_index] == derived) {
                    ret = MSP_ERR_INCONSISTENT_MUTATIONS;
                    goto out;
                }
                genotypes[sample_index] = derived;
                not_done = w != tail;
                w = w->next;
            }
        }
    }
out:
    return ret;
}

int
vargen_next(vargen_t *self, site_t **site, char *genotypes)
{
    int ret = 0;

    bool not_done = true;
    site_t *s;
    char offset = 0;

    if (! (self->flags & MSP_GENOTYPES_AS_CHAR)) {
       offset = '0';
    }
    if (!self->finished) {
        while (not_done && self->tree_site_index == self->tree.sites_length) {
            ret = vargen_next_tree(self);
            if (ret < 0) {
                goto out;
            }
            not_done = ret == 1;
        }
        if (not_done) {
            s = &self->tree.sites[self->tree_site_index];
            ret = vargen_apply_tree_site(self, s, genotypes, offset);
            if (ret != 0) {
                goto out;
            }
            self->tree_site_index++;
            *site = s;
            ret = 1;
        }
    }
out:
    return ret;
}

static int
vargen_update_genotypes(vargen_t *self)
{
    int ret = 0;
    node_list_t *w, *tail;
    node_id_t sample_index;
    bool not_done;
    table_size_t j;
    variant_t *var = &self->variant;
    site_t *site = var->site;
    uint8_t *genotypes = var->genotypes;
    uint8_t derived;

    memset(genotypes, 0, self->num_samples);
    for (j = 0; j < site->mutations_length; j++) {
        assert(j < UINT8_MAX - 1);
        derived = (uint8_t) (j + 1);
        ret = sparse_tree_get_sample_list(&self->tree, site->mutations[j].node, &w, &tail);
        if (ret != 0) {
            goto out;
        }
        if (w != NULL) {
            not_done = true;
            while (not_done) {
                assert(w != NULL);
                sample_index = self->sample_index_map[w->node];
                assert(sample_index >= 0);
                if (genotypes[sample_index] == derived) {
                    ret = MSP_ERR_INCONSISTENT_MUTATIONS;
                    goto out;
                }
                genotypes[sample_index] = derived;
                not_done = w != tail;
                w = w->next;
            }
        }
    }
out:
    return ret;
}

static int
vargen_update_alleles(vargen_t *self)
{
    int ret = 0;
    table_size_t j;
    mutation_t *mutation;
    variant_t *var = &self->variant;

    var->alleles[0] = var->site->ancestral_state;
    var->allele_lengths[0] = var->site->ancestral_state_length;
    var->num_alleles = 1 + var->site->mutations_length;
    for (j = 0; j < var->site->mutations_length; j++) {
        mutation = &var->site->mutations[j];
        var->alleles[j + 1] = mutation->derived_state;
        var->allele_lengths[j + 1] = mutation->derived_state_length;
    }
    return ret;
}

int
vargen_next_dev(vargen_t *self, variant_t **variant)
{
    int ret = 0;

    bool not_done = true;

    if (!self->finished) {
        while (not_done && self->tree_site_index == self->tree.sites_length) {
            ret = vargen_next_tree(self);
            if (ret < 0) {
                goto out;
            }
            not_done = ret == 1;
        }
        if (not_done) {
            self->variant.site = &self->tree.sites[self->tree_site_index];
            ret = vargen_update_alleles(self);
            if (ret != 0) {
                goto out;
            }
            ret = vargen_update_genotypes(self);
            if (ret != 0) {
                goto out;
            }
            self->tree_site_index++;
            *variant = &self->variant;
            ret = 1;
        }
    }
out:
    return ret;
}
