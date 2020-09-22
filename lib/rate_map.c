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

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "util.h"
#include "rate_map.h"

#include <tskit/core.h>

void
rate_map_print_state(rate_map_t *self, FILE *out)
{
    size_t j;
    double rate;

    fprintf(out, "rate_map (%p):: size = %d\n", (void *) self, (int) self->size);
    fprintf(out, "\tsequence_length = %.14g\n", rate_map_get_sequence_length(self));
    fprintf(out, "\tindex\tposition\tvalue\n");
    for (j = 0; j <= self->size; j++) {
        rate = j < self->size ? self->rate[j] : -1;
        fprintf(out, "\t%d\t%.14g\t%.14g\n", (int) j, self->position[j], rate);
    }
}

int MSP_WARN_UNUSED
rate_map_alloc(rate_map_t *self, size_t size, double *position, double *rate)
{
    int ret = 0;
    size_t j;
    double sum;

    memset(self, 0, sizeof(rate_map_t));
    if (size < 1) {
        ret = MSP_ERR_INSUFFICIENT_INTERVALS;
        goto out;
    }
    /* Check the framing position */
    if (position[0] != 0.0) {
        ret = MSP_ERR_INTERVAL_MAP_START_NON_ZERO;
        goto out;
    }
    self->rate = malloc(size * sizeof(*self->rate));
    self->position = malloc((size + 1) * sizeof(*self->position));
    self->cumulative_mass = malloc((size + 1) * sizeof(*self->cumulative_mass));
    if (self->position == NULL || self->rate == NULL || self->cumulative_mass == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->size = size;
    sum = 0;
    for (j = 0; j <= size; j++) {
        if (!isfinite(position[j])) {
            ret = MSP_ERR_NONFINITE_INTERVAL_POSITION;
            goto out;
        }
        self->cumulative_mass[j] = sum;
        self->position[j] = position[j];
        if (j < size) {
            /* This condition and the requirement that position 0 = 0 means that
             * we can't have negative values. */
            if (position[j] >= position[j + 1]) {
                ret = MSP_ERR_INTERVAL_POSITIONS_UNSORTED;
                goto out;
            }
            if ((!isfinite(rate[j])) || rate[j] < 0) {
                ret = MSP_ERR_BAD_RATE_VALUE;
                goto out;
            }
            self->rate[j] = rate[j];
            sum += (position[j + 1] - position[j]) * rate[j];
        }
    }
    fast_search_lookup_alloc(&self->position_lookup, self->position, size + 1);
out:
    return ret;
}

int MSP_WARN_UNUSED
rate_map_alloc_single(rate_map_t *self, double sequence_length, double rate)
{
    double position[2] = { 0, sequence_length };
    return rate_map_alloc(self, 1, position, &rate);
}

int
rate_map_free(rate_map_t *self)
{
    fast_search_lookup_free(&(self->position_lookup));
    msp_safe_free(self->position);
    msp_safe_free(self->rate);
    msp_safe_free(self->cumulative_mass);
    return 0;
}

double
rate_map_get_sequence_length(rate_map_t *self)
{
    return self->position[self->size];
}

size_t
rate_map_get_size(rate_map_t *self)
{
    return self->size;
}

double
rate_map_get_total_mass(rate_map_t *self)
{
    return self->cumulative_mass[self->size];
}

size_t
rate_map_get_index(rate_map_t *self, double x)
{
    size_t index = tsk_search_sorted(self->position, self->size + 1, x);

    if (self->position[index] > x) {
        index--;
    }
    return index;
}

double
rate_map_mass_between(rate_map_t *self, double left, double right)
{
    double left_mass = rate_map_position_to_mass(self, left);
    double right_mass = rate_map_position_to_mass(self, right);
    return right_mass - left_mass;
}

/* Translates a coordinate to the cumulative mass up to that position.
 */
double
rate_map_position_to_mass(rate_map_t *self, double pos)
{
    const double *position = self->position;
    const double *rate = self->rate;
    const double *ptr;
    double offset;
    size_t index;

    /* TODO we should add an LRU cache with ~3 items here to cache recent
     * hits and avoid looking up the map for them. This will make a difference
     * in the large recombination map case.*/
    if (pos == position[0]) {
        return 0;
    }
    if (pos >= position[self->size]) {
        return self->cumulative_mass[self->size];
    }
    ptr = fast_search_lookup_find(&(self->position_lookup), pos);
    index = (size_t)(ptr - position);
    if (index == self->size + 1) {
        index--;
    }
    assert(index > 0);
    index--;
    offset = pos - position[index];

    return self->cumulative_mass[index] + offset * rate[index];
}

/* Finds the physical coordinate such that the sequence up to (but not
 * including) that position has the specified recombination mass.
 */
double
rate_map_mass_to_position(rate_map_t *self, double mass)
{
    const double *position = self->position;
    const double *rate = self->rate;
    double mass_in_interval, pos;
    size_t index;

    if (mass == 0.0) {
        return position[0];
    }
    /* replace with tsk_search_sorted */
    index = msp_binary_interval_search(mass, self->cumulative_mass, self->size + 1);
    assert(index > 0);
    index--;
    mass_in_interval = mass - self->cumulative_mass[index];
    pos = position[index] + mass_in_interval / rate[index];
    return pos;
}

double
rate_map_shift_by_mass(rate_map_t *self, double pos, double mass)
{
    double result_mass = rate_map_position_to_mass(self, pos) + mass;
    return rate_map_mass_to_position(self, result_mass);
}

#ifdef _Static_assert
_Static_assert(FLT_RADIX == 2, "Base 2 floating point types required");
#endif

static int
fast_search_lookup_init_lookups(
    fast_search_lookup_t *self, const double *values, size_t n_values)
{
    int ret = 0;
    const double *lower, *strict_upper;
    uint32_t i;
    double query = -INFINITY;

    lower = values;
    for (i = 0; i < self->num_lookups; i++) {
        while (*lower <= query) {
            lower++;
        }
        self->lookups[i].start = lower;
        strict_upper = lower;
        query = ldexp(i, -(self->power_shift));
        assert(strict_upper < values + n_values);
        while (strict_upper < values + n_values && *strict_upper <= query) {
            strict_upper++;
        }
        const int64_t num = strict_upper - lower;
        if (num > INT32_MAX) {
            ret = MSP_ERR_BAD_PARAM_VALUE; // TODO: consider new err code
            goto out;
        }
        self->lookups[i].num = (size_t) num;
    }
out:
    return ret;
}

static inline long
ceil_to_long(double x)
{
    return (long) ceil(x); // TODO: test use of fast non-portable equivalent
}

int
fast_search_lookup_alloc(
    fast_search_lookup_t *self, const double *values, size_t n_values)
{
    int ret = 0;
    const double max = values[n_values - 1];
    const int power = ilogb((double) n_values) + 1;

    memset(self, 0, sizeof(*self));
    if (*values != 0.0 || !isfinite(max) || max < 0) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    self->power_shift = power - (ilogb(max) + 1);
    const long max_index = ceil_to_long(ldexp(max, self->power_shift));
    self->num_lookups = 1 + (size_t) max_index;
    self->lookups = malloc(self->num_lookups * sizeof(*(self->lookups)));
    if (self->lookups == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = fast_search_lookup_init_lookups(self, values, n_values);
out:
    return ret;
}

int
fast_search_lookup_free(fast_search_lookup_t *self)
{
    msp_safe_free(self->lookups);
    return 0;
}

/*  PRE-CONDITIONS:
 *      1) query >= 0.0
 *      2) self is valid fast_search_lookup_t
 */
const double *
fast_search_lookup_find(fast_search_lookup_t *self, double query)
{
    assert(query >= 0.0);
    search_range_t *back = self->lookups + self->num_lookups - 1;
    if (query > back->start[back->num - 1]) {
        return back->start + back->num;
    }
    size_t i = (size_t) ceil_to_long(ldexp(query, self->power_shift));
    search_range_t *range = self->lookups + i;
    return find_first_upper_bound(range->start, range->num, query);
}
