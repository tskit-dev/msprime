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
#include <gsl/gsl_randist.h>

#include "util.h"
#include "msprime.h"

void
recomb_map_print_state(recomb_map_t *self, FILE *out)
{
    size_t j;

    fprintf(out, "recombination_map (%p):: size = %d\n", (void *) self, (int) self->size);
    fprintf(out, "\tsequence_length = %f\n", recomb_map_get_sequence_length(self));
    fprintf(out, "\tindex\tlocation\trate\n");
    for (j = 0; j < self->size; j++) {
        fprintf(out, "\t%d\t%f\t%f\n", (int) j, self->positions[j], self->rates[j]);
    }
}

static void
recomb_map_init_cumulative_recomb_mass(recomb_map_t *self)
{
    size_t j;
    double s = 0;
    self->cumulative[0] = 0;
    for (j = 1; j < self->size; j++) {
        s += (self->positions[j] - self->positions[j - 1]) * self->rates[j - 1];
        self->cumulative[j] = s;
    }
}

int MSP_WARN_UNUSED
recomb_map_alloc_uniform(recomb_map_t *self, double sequence_length,
        double rate, bool discrete)
{
    double positions[] = {0.0, sequence_length};
    double rates[] = {rate, 0.0};

    return recomb_map_alloc(self, sequence_length, positions, rates, 2, discrete);
}

int MSP_WARN_UNUSED
recomb_map_alloc(recomb_map_t *self, double sequence_length,
        double *positions, double *rates, size_t size, bool discrete)
{
    int ret = MSP_ERR_BAD_RECOMBINATION_MAP;
    double length;
    size_t j;

    memset(self, 0, sizeof(recomb_map_t));
    if (size < 2) {
        goto out;
    }
    /* Check the framing positions */
    if (positions[0] != 0.0 || positions[size - 1] != sequence_length) {
        goto out;
    }
    if (sequence_length < 1 && discrete) {
        goto out;
    }
    self->positions = malloc(size * sizeof(double));
    self->rates = malloc(size * sizeof(double));
    self->cumulative = malloc(size * sizeof(double));

    if (self->positions == NULL || self->rates == NULL || self->cumulative == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->total_recombination_rate = 0.0;
    self->size = size;
    self->sequence_length = sequence_length;
    self->discrete = discrete;
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
    recomb_map_init_cumulative_recomb_mass(self);
out:
    return ret;
}

int
recomb_map_copy(recomb_map_t *to, recomb_map_t *from)
{
    return recomb_map_alloc(to, from->sequence_length,
            from->positions, from->rates, from->size, from->discrete);
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
    if (self->cumulative != NULL) {
        free(self->cumulative);
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

/* Returns the total physical length of the sequence.
 */
double
recomb_map_get_sequence_length(recomb_map_t *self)
{
    return self->sequence_length;
}

/* Returns a boolean indicating whether the recombination map is discrete.
 */
bool
recomb_map_get_discrete(recomb_map_t *self)
{
    return self->discrete;
}

double
recomb_map_mass_between(recomb_map_t *self, double left, double right)
{
    double left_mass = recomb_map_position_to_mass(self, left);
    double right_mass = recomb_map_position_to_mass(self, right);
    return right_mass - left_mass;
}

double
recomb_map_mass_between_left_exclusive(recomb_map_t *self, double left, double right)
{
    double left_bound = self->discrete ? left + 1 : left;
    return recomb_map_mass_between(self, left_bound, right);
}

/* Translates a physical coordinate to the cumulative recombination
 * mass up to (but not including) that position.
 */
double
recomb_map_position_to_mass(recomb_map_t *self, double pos)
{
    double offset;
    size_t index;

    if (pos == self->positions[0]) {
        return 0;
    }
    if (pos >= self->positions[self->size - 1]) {
        return self->cumulative[self->size - 1];
    }
    index = msp_binary_interval_search(pos, self->positions, self->size);
    assert(index > 0);
    index--;
    offset = pos - self->positions[index];

    return self->cumulative[index] + offset * self->rates[index];
}

/* Finds the physical coordinate such that the sequence up to (but not
 * including) that position has the specified recombination mass.
 */
double
recomb_map_mass_to_position(recomb_map_t *self, double mass)
{
    double mass_in_interval, pos;
    size_t index;

    if (mass == 0.0) {
        return self->positions[0];
    }
    index = msp_binary_interval_search(mass, self->cumulative, self->size);
    assert(index > 0);
    index--;
    mass_in_interval = mass - self->cumulative[index];
    pos =  self->positions[index] + mass_in_interval / self->rates[index];

    return self->discrete ? floor(pos) : pos;
}

double
recomb_map_shift_by_mass(recomb_map_t *self, double pos, double mass)
{
    double result_mass = recomb_map_position_to_mass(self, pos) + mass;
    return recomb_map_mass_to_position(self, result_mass);
}

/* Select a position from (start, sequence_length) by sampling
 * a distance in mass to the next recombination site and finding the
 * corresponding position that much after start.
 */
double
recomb_map_sample_poisson(recomb_map_t *self, gsl_rng *rng, double start)
{
    double left_bound, mass_to_next_recomb;

    left_bound = self->discrete ? start + 1 : start;
    mass_to_next_recomb = gsl_ran_exponential(rng, 1.0);

    return recomb_map_shift_by_mass(self, left_bound, mass_to_next_recomb);
}

void
recomb_map_convert_rates(recomb_map_t *self, msp_convert_func convert, void *obj)
{
    size_t i;
    for (i = 0; i < self->size; i++) {
        self->rates[i] = convert(obj, self->rates[i]);
    }
    recomb_map_init_cumulative_recomb_mass(self);
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
