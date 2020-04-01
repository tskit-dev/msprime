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

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "util.h"
#include "msprime.h"

void
recomb_map_print_state(recomb_map_t *self, FILE *out)
{
    fprintf(out, "recombination_map (%p)\n", (void *) self);
    interval_map_print_state(&self->map, out);
}

static int
recomb_map_init_cumulative_recomb_mass(recomb_map_t *self)
{
    int ret = 0;
    size_t j;
    double s = 0;
    const double *position = self->map.position;
    const double *rate = self->map.value;

    self->cumulative[0] = 0;
    for (j = 1; j < self->map.size; j++) {
        if (rate[j - 1] < 0) {
            ret = MSP_ERR_BAD_RECOMBINATION_MAP;
            goto out;
        }
        s += (position[j] - position[j - 1]) * rate[j - 1];
        self->cumulative[j] = s;
    }
out:
    return ret;
}

int MSP_WARN_UNUSED
recomb_map_alloc_uniform(recomb_map_t *self, double sequence_length,
        double rate, bool discrete)
{
    double positions[] = {0.0, sequence_length};
    double rates[] = {rate, 0.0};

    return recomb_map_alloc(self, 2, positions, rates, discrete);
}

int MSP_WARN_UNUSED
recomb_map_alloc(recomb_map_t *self, size_t size, double *position, double *rate,
        bool discrete)
{
    int ret = 0;

    memset(self, 0, sizeof(recomb_map_t));
    self->discrete = discrete;

    ret = interval_map_alloc(&self->map, size, position, rate);
    if (ret != 0) {
        goto out;
    }
    if (interval_map_get_sequence_length(&self->map) < 1 && discrete) {
        ret = MSP_ERR_BAD_RECOMBINATION_MAP;
        goto out;
    }
    self->cumulative = malloc(size * sizeof(double));
    if (self->cumulative == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = recomb_map_init_cumulative_recomb_mass(self);
out:
    return ret;
}

int
recomb_map_copy(recomb_map_t *to, recomb_map_t *from)
{
    return recomb_map_alloc(to, from->map.size,
            from->map.position, from->map.value, from->discrete);
}

int
recomb_map_free(recomb_map_t *self)
{
    interval_map_free(&self->map);
    msp_safe_free(self->cumulative);
    return 0;
}

double
recomb_map_get_total_recombination_rate(recomb_map_t *self)
{
    return self->cumulative[self->map.size - 1];
}

/* Returns the total physical length of the sequence.
 */
double
recomb_map_get_sequence_length(recomb_map_t *self)
{
    return interval_map_get_sequence_length(&self->map);
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
    const double *position = self->map.position;
    const double *rate = self->map.value;
    double offset;
    size_t index;

    if (pos == position[0]) {
        return 0;
    }
    if (pos >= position[self->map.size - 1]) {
        return self->cumulative[self->map.size - 1];
    }
    /* TODO replace this with tsk_search_sorted */
    index = msp_binary_interval_search(pos, position, self->map.size);
    assert(index > 0);
    index--;
    offset = pos - position[index];

    return self->cumulative[index] + offset * rate[index];
}

/* Finds the physical coordinate such that the sequence up to (but not
 * including) that position has the specified recombination mass.
 */
double
recomb_map_mass_to_position(recomb_map_t *self, double mass)
{
    const double *position = self->map.position;
    const double *rate = self->map.value;
    double mass_in_interval, pos;
    size_t index;

    if (mass == 0.0) {
        return position[0];
    }
    index = msp_binary_interval_search(mass, self->cumulative, self->map.size);
    assert(index > 0);
    index--;
    mass_in_interval = mass - self->cumulative[index];
    pos = position[index] + mass_in_interval / rate[index];

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
    do {
        mass_to_next_recomb = gsl_ran_exponential(rng, 1.0);
    } while (mass_to_next_recomb == 0.0);

    return recomb_map_shift_by_mass(self, left_bound, mass_to_next_recomb);
}

void
recomb_map_convert_rates(recomb_map_t *self, msp_convert_func convert, void *obj)
{
    size_t i;
    for (i = 0; i < self->map.size; i++) {
        self->map.value[i] = convert(obj, self->map.value[i]);
    }
    recomb_map_init_cumulative_recomb_mass(self);
}

size_t
recomb_map_get_size(recomb_map_t *self)
{
    return interval_map_get_size(&self->map);
}

int
recomb_map_get_positions(recomb_map_t *self, double *position)
{
    memcpy(position, self->map.position, sizeof(double) * self->map.size);
    return 0;
}

int recomb_map_get_rates(recomb_map_t *self, double *rate)
{
    memcpy(rate, self->map.value, sizeof(double) * self->map.size);
    return 0;
}
