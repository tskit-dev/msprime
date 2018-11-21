/*
** Copyright (C) 2015 University of Oxford
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

#include <gsl/gsl_math.h>

#include "util.h"
#include "msprime.h"

void
recomb_map_print_state(recomb_map_t *self, FILE *out)
{
    size_t j;

    fprintf(out, "recombination_map (%p):: size = %d\n", (void *) self, (int) self->size);
    fprintf(out, "\tnum_loci = %d\n", recomb_map_get_num_loci(self));
    fprintf(out, "\tsequence_length = %f\n", recomb_map_get_sequence_length(self));
    fprintf(out, "\tper_locus_rate = %f\n",
            recomb_map_get_per_locus_recombination_rate(self));
    fprintf(out, "\tindex\tlocation\trate\n");
    for (j = 0; j < self->size; j++) {
        fprintf(out, "\t%d\t%f\t%f\n", (int) j, self->positions[j], self->rates[j]);
    }
}

int WARN_UNUSED
recomb_map_alloc_uniform(recomb_map_t *self, uint32_t num_loci, double sequence_length,
        double rate)
{
    double positions[] = {0.0, sequence_length};
    double rates[] = {rate, 0.0};

    return recomb_map_alloc(self, num_loci, sequence_length, positions, rates, 2);
}

int WARN_UNUSED
recomb_map_alloc(recomb_map_t *self, uint32_t num_loci, double sequence_length,
        double *positions, double *rates, size_t size)
{
    int ret = MSP_ERR_BAD_RECOMBINATION_MAP;
    double length;
    size_t j;

    memset(self, 0, sizeof(recomb_map_t));
    if (size < 2 || num_loci == 0) {
        goto out;
    }
    /* Check the framing positions */
    if (positions[0] != 0.0 || positions[size - 1] != sequence_length) {
        goto out;
    }
    self->positions = malloc(size * sizeof(double));
    self->rates = malloc(size * sizeof(double));
    if (self->positions == NULL || self->rates == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->total_recombination_rate = 0.0;
    self->size = size;
    self->num_loci = num_loci;
    self->sequence_length = sequence_length;
    for (j = 0; j < size; j++) {
        if (rates[j] < 0 || positions[j] < 0) {
            goto out;
        }
        if (j > 0) {
            /* Coordinates must be sorted */
            if (positions[j] <= positions[j - 1]) {
                goto out;
            }
            length = positions[j] - positions[j - 1];
            self->total_recombination_rate += length * rates[j - 1];
        }
        self->rates[j] = rates[j];
        self->positions[j] = positions[j];
    }
    ret = 0;
out:
    return ret;
}

int
recomb_map_free(recomb_map_t *self)
{
    if (self->positions != NULL) {
        free(self->positions);
    }
    if (self->rates != NULL) {
        free(self->rates);
    }
    return 0;
}

/* Returns the equivalent recombination rate between pairs of
 * adjacent loci.
 */
double
recomb_map_get_total_recombination_rate(recomb_map_t *self)
{
    return self->total_recombination_rate;
}

/* Returns the equivalent recombination rate between pairs of
 * adjacent loci.
 */
double
recomb_map_get_per_locus_recombination_rate(recomb_map_t *self)
{
    double ret = 0.0;
    if (self->num_loci > 1) {
        ret = self->total_recombination_rate / (self->num_loci - 1);
    }
    return ret;
}

/* Returns the total physical length of the sequence.
 */
double
recomb_map_get_sequence_length(recomb_map_t *self)
{
    return self->sequence_length;
}

/* Returns the total number of discrete loci, between which
 * recombination can occur.
 */
uint32_t
recomb_map_get_num_loci(recomb_map_t *self)
{
    return self->num_loci;
}

/* Remaps the specified physical coordinate in the range (0, sequence_length)
 * to the genetic coordinate space in the range (0, num_loci)
 */
double
recomb_map_phys_to_genetic(recomb_map_t *self, double x)
{
    size_t j;
    double s = 0.0;
    double ret = 0.0;
    double rate = 1.0;
    double last_phys_x, phys_x;

    if (self->total_recombination_rate == 0) {
        /* When the total recombination rate is zero anything within the interval
         * maps to 0. L maps to num_loci */
        ret = x >= self->sequence_length? 1: 0;
    } else {
        last_phys_x = 0;
        for (j = 1; j < self->size && x > self->positions[j]; j++) {
            phys_x = self->positions[j];
            rate = self->rates[j - 1];
            s += (phys_x - last_phys_x) * rate;
            last_phys_x = phys_x;
        }
        rate = self->rates[j - 1];
        s += (x - last_phys_x) * rate;
        assert(s >= 0 && s <= self->total_recombination_rate);
        ret = s / self->total_recombination_rate;
    }
    return ret * self->num_loci;
}

/* Remaps the specified coordinate in physical space into a discrete genetic
 * coordinate. */
int
recomb_map_phys_to_discrete_genetic(recomb_map_t *self, double x, uint32_t *locus)
{
    int ret = 0;

    if (x < 0 || x > self->sequence_length) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    *locus = (uint32_t) round(recomb_map_phys_to_genetic(self, x));
out:
    return ret;
}

/* Remaps the specified genetic coordinate in the range (0, num_loci) to
 * the physical coordinate space in the range (0, sequence_length)
 */
double
recomb_map_genetic_to_phys(recomb_map_t *self, double genetic_x)
{
    size_t k;
    double ret = 0.0;
    double x, s, excess;
    double *p = self->positions;
    double *r = self->rates;
    double num_loci = self->num_loci;

    assert(num_loci >= 1);
    assert(genetic_x >= 0);
    assert(genetic_x <= num_loci);
    if (self->total_recombination_rate == 0) {
        /* When we have a 0 total rate, anything other than m maps to 0. */
        ret = genetic_x >= self->num_loci? self->sequence_length: 0;
    } else if (self->size == 2) {
        /* Avoid roundoff when num_loci == self->sequence_length */
        ret = genetic_x;
        if (self->sequence_length != num_loci) {
            ret = (genetic_x / num_loci) * self->sequence_length;
        }
    } else {
        /* genetic_x is in the range [0,num_loci], and so we rescale
         * this into [0,total_recombination_rate] so that we can
         * map back into physical coordinates. */
        x = (genetic_x / num_loci) * self->total_recombination_rate;
        if (x > 0) {
            s = 0;
            k = 0;
            while (s < x && k < self->size - 1) {
                s += (p[k + 1] - p[k]) * r[k];
                k++;
            }
            assert(k > 0);
            excess = 0;
            if (r[k - 1] > 0) {
                excess = (s - x) / r[k - 1];
            }
            ret = p[k] - excess;
        }
    }
    return ret;
}

size_t
recomb_map_get_size(recomb_map_t *self)
{
    return self->size;
}

int
recomb_map_get_positions(recomb_map_t *self, double *positions)
{
    memcpy(positions, self->positions, sizeof(double) * self->size);
    return 0;
}

int recomb_map_get_rates(recomb_map_t *self, double *rates)
{
    memcpy(rates, self->rates, sizeof(double) * self->size);
    return 0;

}
