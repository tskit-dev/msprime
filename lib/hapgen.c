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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "err.h"
#include "msprime.h"


int
hapgen_alloc(hapgen_t *self, double mutation_rate, tree_file_t * tree_file,
        long random_seed, size_t max_haplotype_length)
{
    int ret = -1;
    uint32_t j;
    uint32_t N;

    memset(self, 0, sizeof(hapgen_t));
    assert(tree_file != NULL);
    self->mutation_rate = mutation_rate;
    self->tree_file = tree_file;
    self->random_seed = random_seed;
    self->max_haplotype_length = max_haplotype_length;
    self->rng = gsl_rng_alloc(gsl_rng_default);
    if (self->rng == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    gsl_rng_set(self->rng, self->random_seed);
    N  = 2 * tree_file->sample_size;
    self->pi = malloc(N * sizeof(int));
    self->tau = malloc(N * sizeof(float));
    if (self->pi == NULL || self->tau == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->haplotypes = malloc(N * sizeof(char *));
    if (self->haplotypes == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < N; j++) {
        /* TODO make one malloc of this and split it manually */
        self->haplotypes[j] = malloc(self->max_haplotype_length);
        if (self->haplotypes[j]== NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
    }
out:
    return ret;
}

static int
hapgen_generate_haplotypes(hapgen_t *self, uint32_t l)
{
    int ret = -1;
    char *h, *hp;
    int N = 2 * self->tree_file->sample_size;
    int *child = calloc(N, sizeof(int));
    int *sib = calloc(N, sizeof(int));
    int *branch_mutations = malloc(N * sizeof(int));
    int *mutation_sites = malloc(N * sizeof(int));
    int u, v, j, t, s, h_size, not_done;
    double mu = self->mutation_rate * l;

    /* Generate the mutations */
    s = 0;
    for (j = 1; j < N; j++) {
        t = self->tau[self->pi[j]] - self->tau[j];
        branch_mutations[j] = gsl_ran_poisson(self->rng, mu * t);
        mutation_sites[j] = s;
        s += branch_mutations[j];
    }
    if (s > self->max_haplotype_length) {
        ret = MSP_ERR_TOO_MANY_SEG_SITES;
        goto out;
    }
    self->haplotype_length = s;
    /* allocate the haplotypes */
    h_size = s + 1;
    for (j = 0; j < N; j++) {
        memset(self->haplotypes[j], '0', h_size - 1);
        self->haplotypes[j][h_size - 1] = '\0';
    }
    /* Now traverse the tree and apply these mutations */
    for (u = 1; u < N; u++) {
        v = self->pi[u];
        sib[u] = child[v];
        child[v] = u;
    }
    u = child[0];
    while (u != 0) {
        not_done = 1;
        while (not_done) {
            if (child[u] != 0) {
                v = u;
                u = child[u];
                h = self->haplotypes[v];
                hp = self->haplotypes[u];
                strncpy(hp, h, s + 1);
                for (j = 0; j < branch_mutations[u]; j++) {
                    hp[mutation_sites[u] + j] = '1';
                }
            } else {
                not_done = 0;
            }
        }
        not_done = 1;
        while (not_done) {
            if (sib[u] != 0) {
                u = sib[u];
                v = self->pi[u];
                h = self->haplotypes[v];
                hp = self->haplotypes[u];
                strncpy(hp, h, s + 1);
                for (j = 0; j < branch_mutations[u]; j++) {
                    hp[mutation_sites[u] + j] = '1';
                }
                not_done = 0;
            } else {
                u = self->pi[u];
                not_done = u != 0;
            }
        }
    }
    free(child);
    free(sib);
    free(mutation_sites);
    free(branch_mutations);
out:
    return ret;
}


int
hapgen_next(hapgen_t *self, uint32_t *length, char ***haplotypes, size_t *s, int** pi,
        float **tau)
{
    int ret = -1;
    int v = 1;
    uint32_t l;
    coalescence_record_t *cr = &self->next_record;
    uint32_t b = cr->left;


    /* TODO this is completely untested!! Put in an interface in main to exercise
     * the code and see how it goes under valgrind.
     */

    /* We've not read any coalescence records yet, this the first call */
    if (b == 0) {
        v = tree_file_next_record(self->tree_file, cr);
        if (v < 0) {
            ret = v;
            goto out;
        }
        b = cr->left;
    }
    while (cr->left == b && v == 1) {
        self->pi[cr->children[0]] = cr->parent;
        self->pi[cr->children[1]] = cr->parent;
        self->tau[cr->parent] = cr->time;
        v = tree_file_next_record(self->tree_file, cr);
        if (v < 0) {
            ret = v;
            goto out;
        }
    }
    if (v == 0) {
        l = self->tree_file->num_loci - cr->left + 1;
    } else {
        l = cr->left - b;
    }
    ret = hapgen_generate_haplotypes(self, l);
    if (ret < 0) {
        goto out;
    }
    /* if we've read the last coalescence record then we're done; otherwise
     * indicate that there are trees left, as per the protocol for records
     */
    ret = v;
    *haplotypes = self->haplotypes;
    *s = self->haplotype_length;
    *pi = self->pi;
    *tau = self->tau;
out:
    return ret;
}

int
hapgen_free(hapgen_t *self)
{

    return 0;
}
