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
hapgen_alloc(hapgen_t *self, double mutation_rate, const char *tree_file_name,
        long random_seed, size_t max_haplotype_length)
{
    int ret = -1;
    uint32_t j;
    uint32_t N;

    memset(self, 0, sizeof(hapgen_t));
    self->mutation_rate = mutation_rate;
    self->random_seed = random_seed;
    self->max_haplotype_length = max_haplotype_length;
    self->rng = gsl_rng_alloc(gsl_rng_default);
    if (self->rng == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    gsl_rng_set(self->rng, self->random_seed);
    ret = tree_file_open(&self->tree_file, tree_file_name, 'r');
    if (ret != 0) {
        goto out;
    }
    if (!tree_file_issorted(&self->tree_file)) {
        ret = MSP_ERR_TREE_FILE_NOT_SORTED;
        goto out;
    }
    self->sample_size = self->tree_file.sample_size;
    self->num_loci = self->tree_file.num_loci;
    N  = 2 * self->sample_size;
    self->pi = calloc(N, sizeof(int));
    self->tau = calloc(N, sizeof(float));
    self->child = malloc(N * sizeof(int));
    self->sib = malloc(N * sizeof(int));
    self->branch_mutations = malloc(N * sizeof(int));
    self->mutation_sites = malloc(N * sizeof(int));
    if (self->pi == NULL || self->tau == NULL || self->child == NULL
            || self->sib == NULL || self->branch_mutations == NULL
            || self->mutation_sites == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->haplotypes = malloc(N * sizeof(char *));
    self->haplotype_mem = malloc(N * self->max_haplotype_length);
    if (self->haplotypes == NULL || self->haplotype_mem == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < N; j++) {
        self->haplotypes[j] = self->haplotype_mem
            + j * self->max_haplotype_length;
    }
    /* initialise all haplotypes to '0' */
    memset(self->haplotype_mem, '0', N * self->max_haplotype_length);
    ret = 0;
out:
    return ret;
}

int
hapgen_free(hapgen_t *self)
{
    if (self->rng != NULL) {
        gsl_rng_free(self->rng);
    }
    if (self->pi != NULL) {
        free(self->pi);
    }
    if (self->tau != NULL) {
        free(self->tau);
    }
    if (self->child != NULL) {
        free(self->child);
    }
    if (self->sib != NULL) {
        free(self->sib);
    }
    if (self->branch_mutations != NULL) {
        free(self->branch_mutations);
    }
    if (self->mutation_sites != NULL) {
        free(self->mutation_sites);
    }
    if (self->haplotypes != NULL) {
        free(self->haplotypes);
    }
    if (self->haplotype_mem != NULL) {
        free(self->haplotype_mem);
    }
    tree_file_close(&self->tree_file);
    return 0;
}
static int
hapgen_process_tree(hapgen_t *self, uint32_t l)
{
    int ret = -1;
    char *h, *hp;
    int N = 2 * self->sample_size;
    int *child = self->child;
    int *sib = self->sib;
    int *branch_mutations = self->branch_mutations;
    int *mutation_sites = self->mutation_sites;
    int u, v, j, s, not_done;
    double mu = (self->mutation_rate * l) / self->num_loci;
    double t;

    /* This method is can be probably be improved quite a lot for the
     * sparse mutation case.
     */
    memset(child, 0, N * sizeof(int));
    memset(sib, 0, N * sizeof(int));
    memset(branch_mutations, 0, N * sizeof(int));
    memset(mutation_sites, 0, N * sizeof(int));
    /* Generate the mutations */
    s = 0;
    for (j = 1; j < N - 1; j++) {
        t = self->tau[self->pi[j]] - self->tau[j];
        branch_mutations[j] = gsl_ran_poisson(self->rng, mu * t);
        mutation_sites[j] = s;
        s += branch_mutations[j];
    }
    if (s + self->haplotype_length >= self->max_haplotype_length) {
        ret = MSP_ERR_TOO_MANY_SEG_SITES;
        goto out;
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
                h = self->haplotypes[v] + self->haplotype_length;
                hp = self->haplotypes[u] + self->haplotype_length;
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
                h = self->haplotypes[v] + self->haplotype_length;
                hp = self->haplotypes[u] + self->haplotype_length;
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
    self->haplotype_length += s;
    ret = 0;
out:
    return ret;
}


int
hapgen_generate(hapgen_t *self)
{
    int ret = -1;
    uint32_t b, j;
    coalescence_record_t cr;

    b = 1;
    while ((ret = tree_file_next_record(&self->tree_file, &cr)) == 1) {
        if (cr.left != b) {
            ret = hapgen_process_tree(self, cr.left - b);
            if (ret < 0) {
                goto out;
            }
            b = cr.left;
        }
        self->pi[cr.children[0]] = cr.parent;
        self->pi[cr.children[1]] = cr.parent;
        self->tau[cr.parent] = cr.time;
    }
    if (ret != 0) {
        goto out;
    }
    ret = hapgen_process_tree(self, self->num_loci - cr.left + 1);
    if (ret < 0) {
        goto out;
    }
    /* now finish each haplotype with '\0' so its a valid string */
    for (j = 0; j < 2 * self->sample_size; j++) {
        self->haplotypes[j][self->haplotype_length] = '\0';
    }
out:
    return ret;
}

int
hapgen_get_haplotypes(hapgen_t *self, char ***haplotypes, size_t *s)
{
    *haplotypes = self->haplotypes;
    *s = self->haplotype_length;
    return 0;
}
