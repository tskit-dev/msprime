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
    size_t j;

    fprintf(out, "recombination_map (%p)\n", (void *) self);
    fprintf(out, "total_mass  = %g\n", self->total_mass);
    fprintf(out, "mass_scale  = %gf\n", self->mass_scale);
    fprintf(
        out, "sequence_length  = %f\n", interval_map_get_sequence_length(&self->map));
    for (j = 0; j < self->map.size; j++) {
        fprintf(out, "\t%d\t%g\t%g\t%" PRId64 "\n", (int) j, self->map.position[j],
            self->map.value[j], self->cumulative_scaled_mass[j]);
    }
    interval_map_print_state(&self->map, out);
}

inline static int64_t
recomb_map_scale_mass(recomb_map_t *self, double mass)
{
    double scaled_value = round(mass * self->mass_scale);
    return (int64_t) scaled_value;
}

inline static double
recomb_map_unscale_mass(recomb_map_t *self, int64_t scaled_mass)
{
    return ((double) scaled_mass) / self->mass_scale;
}

/* Sets scaling constant used to convert floating point recombination
 * mass values into integers.
 */
int MSP_WARN_UNUSED
recomb_map_set_mass_scale(recomb_map_t *self, double mass_scale)
{
    int ret = 0;
    size_t j;
    double mass;
    int64_t sum = 0;
    const double *position = self->map.position;
    const double *rates = self->map.value;

    if (mass_scale <= 0 || !isfinite(mass_scale)) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->mass_scale = mass_scale;

    self->cumulative_scaled_mass[0] = 0;
    for (j = 1; j < self->map.size; j++) {
        mass = (position[j] - position[j - 1]) * rates[j - 1];
        sum += recomb_map_scale_mass(self, mass);
        assert(sum >= 0);
        self->cumulative_scaled_mass[j] = sum;
    }
out:
    return ret;
}

int MSP_WARN_UNUSED
recomb_map_alloc_uniform(
    recomb_map_t *self, double sequence_length, double rate, bool discrete)
{
    double positions[] = { 0.0, sequence_length };
    double rates[] = { rate, 0.0 };

    return recomb_map_alloc(self, 2, positions, rates, discrete);
}

int MSP_WARN_UNUSED
recomb_map_alloc(
    recomb_map_t *self, size_t size, double *position, double *rate, bool discrete)
{
    int ret = 0;
    size_t j;
    double x, mass;

    memset(self, 0, sizeof(recomb_map_t));
    self->discrete = discrete;
    self->mass_scale = 0;

    ret = interval_map_alloc(&self->map, size, position, rate);
    if (ret != 0) {
        goto out;
    }
    if (interval_map_get_sequence_length(&self->map) < 1 && discrete) {
        ret = MSP_ERR_BAD_RECOMBINATION_MAP;
        goto out;
    }
    self->cumulative_scaled_mass = calloc(size, sizeof(*self->cumulative_scaled_mass));
    if (self->cumulative_scaled_mass == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->total_mass = 0;
    for (j = 1; j < self->map.size; j++) {
        x = rate[j - 1];
        if ((!isfinite(x)) || x < 0) {
            ret = MSP_ERR_BAD_RECOMBINATION_MAP;
            goto out;
        }
        mass = (position[j] - position[j - 1]) * x;
        self->total_mass += mass;
    }
out:
    return ret;
}

int
recomb_map_copy(recomb_map_t *to, recomb_map_t *from)
{
    return recomb_map_alloc(
        to, from->map.size, from->map.position, from->map.value, from->discrete);
}

int
recomb_map_free(recomb_map_t *self)
{
    interval_map_free(&self->map);
    msp_safe_free(self->cumulative_scaled_mass);
    return 0;
}

double
recomb_map_get_total_mass(recomb_map_t *self)
{
    return self->total_mass;
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

int64_t
recomb_map_scaled_mass_between(recomb_map_t *self, double left, double right)
{
    int64_t left_scaled_mass = recomb_map_position_to_scaled_mass(self, left);
    int64_t right_scaled_mass = recomb_map_position_to_scaled_mass(self, right);
    return right_scaled_mass - left_scaled_mass;
}

int64_t
recomb_map_scaled_mass_between_left_exclusive(
    recomb_map_t *self, double left, double right)
{
    double left_bound = self->discrete ? left + 1 : left;
    return recomb_map_scaled_mass_between(self, left_bound, right);
}

/* Translates a physical coordinate to the correspoding cumulative scaled
 * mass up to (but not including) that position.
 */
int64_t
recomb_map_position_to_scaled_mass(recomb_map_t *self, double pos)
{
    const double *position = self->map.position;
    const double *rate = self->map.value;
    double mass_diff;
    int64_t scaled_mass_diff;
    size_t index;

    if (pos == position[0]) {
        return 0;
    }
    if (pos >= position[self->map.size - 1]) {
        return self->cumulative_scaled_mass[self->map.size - 1];
    }
    index = msp_binary_search_double(pos, position, self->map.size);
    assert(index > 0);
    index--;
    mass_diff = (pos - position[index]) * rate[index];
    scaled_mass_diff = recomb_map_scale_mass(self, mass_diff);
    return self->cumulative_scaled_mass[index] + scaled_mass_diff;
}

/* Finds the physical coordinate such that the sequence up to (but not
 * including) that position has the specified recombination mass.
 */
double
recomb_map_scaled_mass_to_position(recomb_map_t *self, int64_t mass)
{
    const double *position = self->map.position;
    const double *rate = self->map.value;
    int64_t scaled_mass_in_interval;
    double pos, mass_in_interval;
    size_t index;

    if (mass == 0) {
        return position[0];
    }
    index = msp_binary_search_int64(mass, self->cumulative_scaled_mass, self->map.size);
    assert(index > 0);
    index--;
    scaled_mass_in_interval = mass - self->cumulative_scaled_mass[index];
    mass_in_interval = recomb_map_unscale_mass(self, scaled_mass_in_interval);
    pos = position[index] + mass_in_interval / rate[index];
    return self->discrete ? floor(pos) : pos;
}

double
recomb_map_shift_by_scaled_mass(recomb_map_t *self, double pos, int64_t scaled_mass)
{
    int64_t pos_scaled_mass = recomb_map_position_to_scaled_mass(self, pos);
    uint64_t target_mass = (uint64_t) pos_scaled_mass + (uint64_t) scaled_mass;
    uint64_t max_mass = (uint64_t) self->cumulative_scaled_mass[self->map.size - 1];
    double ret = recomb_map_get_sequence_length(self);

    assert(pos_scaled_mass >= 0);
    assert(scaled_mass >= 0);

    /* Protect against overflow and make sure we return L when very large
     * scaled mass values are provided */
    if (target_mass < max_mass) {
        ret = recomb_map_scaled_mass_to_position(self, (int64_t) target_mass);
    }
    return ret;
}

double
recomb_map_shift_by_mass(recomb_map_t *self, double pos, double mass)
{
    return recomb_map_shift_by_scaled_mass(self, pos, recomb_map_scale_mass(self, mass));
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

int
recomb_map_get_rates(recomb_map_t *self, double *rate)
{
    memcpy(rate, self->map.value, sizeof(double) * self->map.size);
    return 0;
}
