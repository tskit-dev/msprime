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
    /* search for upper bounds strictly before the max position
       any position greather than or equal to max position returns self->size */
    ret = fast_search_lookup_alloc(&self->position_lookup, self->position, size);
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

    assert(pos >= 0.0);
    if (pos <= 0.0) {
        return 0;
    }
    assert(pos <= position[self->size]);
    ptr = fast_search_lookup_find(&(self->position_lookup), pos);
    /* any `pos` greather than or equal to max position has `index == self->size` */
    index = (size_t)(ptr - position);
    assert(index == msp_binary_interval_search(pos, position, self->size));
    assert(index > 0);
    index--;
    offset = pos - position[index];

    return self->cumulative_mass[index] + offset * rate[index];
}

#ifdef EXTRAVAGANT_ASSERTS
static size_t
slow_emulate_msp_binary_interval_search(
    double query, const double *values, size_t n_values)
{
    fast_search_lookup_t table;
    const double *ptr;

    assert(0 == fast_search_lookup_alloc(&table, values, n_values));
    ptr = fast_search_lookup_find(&table, query);
    fast_search_lookup_free(&table);
    return (size_t)(ptr - values);
}
#endif

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

    assert(mass >= 0.0);
    if (mass <= 0.0) {
        return position[0];
    }
    /* search for upper bounds strictly before the final cum mass
       any mass greather than or equal to final cum mass returns self->size */
    index = msp_binary_interval_search(mass, self->cumulative_mass, self->size);
#ifdef EXTRAVAGANT_ASSERTS
    assert(index
           == slow_emulate_msp_binary_interval_search(
                  mass, self->cumulative_mass, self->size));
#endif
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

#if FLT_RADIX != 2
#error "Base 2 floating point types required"
#endif

static bool
valid_sorted_nonempty_array(const double *array, size_t size)
{
    size_t idx;

    if ((ptrdiff_t) size < 1 || isnan(array[0])) {
        return false;
    }
    for (idx = 1; idx < size; idx++) {
        if (!(array[idx - 1] <= array[idx])) {
            return false;
        }
    }
    return true;
}

static bool
fast_search_lookup_valid(fast_search_lookup_t *self)
{
    const double *start, *stop;
    size_t idx;

    if (!(self->query_multiplier >= 0.0)) { // NaN not valid
        return false;
    }
    if (self->num_lookups < 2) {
        return false;
    }
    for (idx = 1; idx < self->num_lookups; idx++) {
        if (self->lookups[idx - 1] > self->lookups[idx]) {
            return false;
        }
    }
    start = self->lookups[0];
    stop = self->lookups[self->num_lookups - 1];
    if (stop <= start || start[0] != 0.0) {
        return false;
    }
    return valid_sorted_nonempty_array(start, (size_t)(stop - start));
}

static int
fast_search_lookup_init_lookups(
    fast_search_lookup_t *self, const double *elements, size_t n_elements)
{
    int ret = 0;
    const double *ptr, *stop;
    uint32_t i;
    double query = 0.0;

    ptr = elements;
    stop = elements + n_elements;
    for (i = 0; i < self->num_lookups; i++) {
        while (ptr < stop && *ptr < query) {
            ptr++;
        }
        self->lookups[i] = ptr;
        if (ptr < stop) {
            query = (i + 1) / self->query_multiplier;
        }
    }
    if (!fast_search_lookup_valid(self)) {
        ret = MSP_ERR_ASSERTION_FAILED;
        goto out;
    }
out:
    return ret;
}

/* Returns the least power of 2 higher than (or equal to) `x` for positive numbers,
 * otherwise returns zero.
 */
static double
higher_power_of_2(double x)
{
    assert(x >= 0);
    return (x > 0 ? exp2(ceil(logb(x))) : 0.0);
}

/* PRE-CONDITIONS:
 *     1) `elements` must point to array of doubles starting at exactly `0.0`
 *     2) `elements` must be non-empty non-decreasing and with no NaN
 */
int
fast_search_lookup_alloc(
    fast_search_lookup_t *self, const double *elements, size_t n_elements)
{
    int ret = 0;
    double max_element;

    const size_t max_input_size = 1L << (DBL_MANT_DIG - 1); // 4096 terabytes

    memset(self, 0, sizeof(*self));

    if (!valid_sorted_nonempty_array(elements, n_elements)) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (n_elements >= max_input_size) {
        ret = MSP_ERR_ASSERTION_FAILED;
        goto out;
    }
    if (elements[0] != 0.0 || isnan(elements[0])) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    max_element = elements[n_elements - 1];
    if (max_element < 0.0 || !isfinite(max_element)) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }

    self->query_multiplier
        = higher_power_of_2((double) n_elements - 1) / higher_power_of_2(max_element);
    if (!isfinite(self->query_multiplier)) {
        self->query_multiplier = exp2(DBL_MAX_EXP - 1); // largest possible power of 2
    }
    assert(isfinite(self->query_multiplier));

    self->num_lookups = 2 + (size_t)(max_element * self->query_multiplier);

    self->lookups = malloc(self->num_lookups * sizeof(*(self->lookups)));
    if (self->lookups == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    ret = fast_search_lookup_init_lookups(self, elements, n_elements);
out:
    return ret;
}

int
fast_search_lookup_free(fast_search_lookup_t *self)
{
    msp_safe_free(self->lookups);
    return 0;
}

static inline const double *
ptr_first_upper_bound(const double *start, const double *stop, double query)
{
    assert(start <= stop);
    return start + msp_binary_interval_search(query, start, (size_t)(stop - start));
}

/*  PRE-CONDITIONS:
 *      1) query >= 0.0
 *      2) self is valid fast_search_lookup_t
 */
const double *
fast_search_lookup_find(fast_search_lookup_t *self, double query)
{
    const double **lookups = self->lookups;
    size_t idx;
    const double *ret;

    assert(query >= 0.0);
    idx = (size_t)(query * self->query_multiplier);
    if (idx + 1 >= self->num_lookups) {
        ret = lookups[self->num_lookups - 1];
    } else {
        ret = ptr_first_upper_bound(lookups[idx], lookups[idx + 1], query);
    }
    assert(
        ret == ptr_first_upper_bound(lookups[0], lookups[self->num_lookups - 1], query));
    return ret;
}
