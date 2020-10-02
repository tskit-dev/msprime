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
    self->rate = malloc((size + 1) * sizeof(*self->rate));
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
    /* Consider the rate after the max position as zero.
       This trailing zero simplifies calcs like rate_map_position_to_mass by
       gracefully handling positions past the max position */
    self->rate[size] = 0.0;
    ret = fast_search_alloc(&self->position_lookup, self->position, size + 1);
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
    fast_search_free(&self->position_lookup);
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
    double base, offset;
    size_t index;

    if (pos <= 0.0) {
        return 0;
    }
    assert(pos <= position[self->size]);
    index = fast_search_idx_strict_upper(&self->position_lookup, pos);
    assert(index == idx_1st_strict_upper_bound(position, self->size + 1, pos));
    assert(index > 0);
    index--;
    base = self->cumulative_mass[index];
    offset = pos - position[index];
    if (offset <= 0) {
        return base;
    }
    return base + offset * rate[index];
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

    assert(mass >= 0.0);
    if (mass <= 0.0) {
        return position[0];
    }
    /* search for upper bounds strictly before the final cum mass
       any mass greather than or equal to final cum mass returns self->size */
    index = idx_1st_upper_bound(self->cumulative_mass, self->size, mass);
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

/* The `fast_search_t` structure and functions
 *
 * `fast_search_t` is a lookup table pointing to an array of `double`s
 * (the "elements"). It does not own the array but does point to all its memory.
 *
 * Its purpose is to speed up the search for a "query" value in that array of
 * elements. The `fast_search_idx_strict_upper` function is essentially a drop-in
 * replacement for the function `idx_1st_strict_upper_bound` (which in turn is
 * a drop-in replacement for the C++ <algorithm> library function
 * std::upper_bound).
 *
 * The key concept behind `fast_search_t` is that the bits inside a
 * "query" `double` contain a good deal of information about where near-by
 * values in the array of elements might be found in memory. The best
 * performance is when the array of elements have values that increase roughly
 * linearly. But this is not a requirement. The only requirement is that the
 * elements be non-decreasing.
 *
 * Since the element values are not perfectly linear some re-mapping is
 * required.  This is achieved through an array of "lookup" pointers into the
 * array of elements.  Although there need not be a linear mapping between
 * element values and the range of query values, the "lookup" pointer are
 * linear with the range of query values.
 *
 * The `query_multiplier` determines the linear mapping between possible query
 * values and indexes in the array of "lookup" pointers. Because the
 * `query_multiplier` is a power of 2, a predictable range of possible "query"
 * values will truncate to a specific lookup index.
 *
 * The number of lookup pointers could be adjustable, but this simple
 * implementation chooses one simple strategy for the number of lookup
 * pointers. For a perfectly linear sequence of element values, the lookup
 * pointers will be a simple element-address-by-next-element-address
 * progression of pointers. The extremely simple case of element values being
 * the non-negative integers is coded in the `test_fast_search_identity`
 * unit test.  An example of alternative linear sequences with different
 * `query_multiplier` are coded in `test_fast_search_2powers`.
 *
 * Consider example target array of elements (memory address, values):
 *
 * 0xECE0, 0.0
 * 0xECE1, 0.1
 * 0xECE2, 0.8
 * 0xECE3, 2.7
 * 0xECE4, 6.4
 * 0xECE5, ___ (the end, past last element of the array)
 *
 * The power of 2 greater or equal to the maximum element is 8.
 * The number of memory address steps to the maximum element is 4.
 * The power of 2 greater or equal to that number of steps is also 4.
 * Thus the query multiplier is going to be 0.5.
 * The lookup table is going to start with a pointer to the last of the elements.
 * The last pointer in the lookup table will point to the end of the element array.
 * The lookup table should be just the right size so that rescaling (and truncating) a
 * query equaling the maximum element will result in an index to the last pointer
 * in the lookup table. This is `trunc(6.4 * 0.5) = 3`. Thus the lookup table
 * will be `(1 + 1 + 3) = 5` pointers.
 *
 * The lookup table of pointers will handle the following 5 query ranges:
 *
 * [0.0, 2.0)
 * [2.0, 4.0)
 * [4.0, 6.0)
 * [6.0, 8.0)
 * [8.0, INFINITY)
 *
 * with each lookup pointer pointing to the first upper bound to the minimum in those
 * query ranges. Thus we get the following lookup table of pointers:
 *
 * [0] 0xECE0 (0.0 in elements is upper bound to query range minimum 0.0)
 * [1] 0xECE3 (2.7 in elements is upper bound to query range minimum 2.0)
 * [2] 0xECE4 (6.4 in elements is upper bound to query range minimum 4.0)
 * [4] 0xECE4 (6.4 in elements is upper bound to query range minimum 6.0)
 * [5] 0xECE5 (no upper bound in elements to query range minimum 8.0)
 *
 * Consider the search with a query of `3.3`. This query will be multiplied by
 * the query multiplier of 0.5 and then truncated to get lookup index `1`. This
 * corresponds the query range [2.0, 4.0) and all queries in that range when
 * multiplied by the query multiplier and truncated will result in lookup index
 * `1`. This entry in the lookup table maps to 0xECE3 which contains the value
 * 2.7. The next entry in the lookup table corresponds to the query range [4.0, 6.0)
 * and has 0xECE4 pointing to element value 6.4. Thus a much smaller binary
 * search will be performed on the memory range [0xECE3, 0xECE4) which will
 * result in finding 0xECE4 as the first strict upper bound to the query of `3.3`.
 */

#if FLT_RADIX != 2
#error "Base 2 floating point types required"
#endif

#if SIZE_MAX < ULLONG_MAX
#error "size_t must be at least 64 bits"
/* unsigned long long must be at least 64 bits */
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

static int
fast_search_init_lookups(fast_search_t *self, const double *elements, size_t n_elements)
{
    int ret = 0;
    const double *ptr, *stop;
    uint32_t idx;
    double min_query;

    self->lookups[0] = elements;
    ptr = elements;
    stop = elements + n_elements;
    for (idx = 1; idx < self->num_lookups; idx++) {
        /* `min_query` is the smallest possible query that will get rescaled to idx */
        min_query = idx / self->query_multiplier;
        /* move ptr to the first upper bound of min_query in elements */
        while (ptr < stop && *ptr < min_query) {
            ptr++;
        }
        /* The query range of [A, B) maps to `idx` where
         * A = idx/query_multiplier and B = (idx + 1)/query_multiplier).
         * Want `lookup[idx]` to point to the first upper bound of A in the elements.
         */
        self->lookups[idx] = ptr;
    }
    return ret;
}

/* Returns the least power of 2 greater than or equal to `x` for positive numbers,
 * otherwise returns zero.
 */
static double
higher_power_of_2(double x)
{
    assert(x >= 0);
    if (x <= 0) {
        return 0.0;
    }
    double floor = exp2(logb(x));
    return (floor < x ? floor * 2.0 : floor);
}

/* PRE-CONDITIONS:
 *     1) `elements` must point to array of doubles starting at exactly `0.0`
 *     2) `elements` must be non-empty non-decreasing and with no NaN
 */
int
fast_search_alloc(fast_search_t *self, const double *elements, size_t n_elements)
{
    int ret = 0;
    double max_element;

    const size_t max_input_size = 1LL << (DBL_MANT_DIG - 1); // 4096 terabytes

    memset(self, 0, sizeof(*self));

    if (!valid_sorted_nonempty_array(elements, n_elements)) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (elements[0] != 0.0 || n_elements >= max_input_size) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    max_element = elements[n_elements - 1];
    if (max_element < 0.0 || !isfinite(max_element)) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }

    /* query_multiplier will rescale queries to be a desired lookup table index
     * (plus a fraction). To avoid rounding problems query_multiplier is a power of 2.
     */
    self->query_multiplier
        = higher_power_of_2((double) n_elements - 1) / higher_power_of_2(max_element);
    if (!isfinite(self->query_multiplier)) {
        self->query_multiplier = exp2(DBL_MAX_EXP - 1); // largest possible power of 2
    }
    assert(isfinite(self->query_multiplier));

    /* Many different sizes for the lookup table are reasonable. The lookup
     * table size choice here is a simple one that will use roughly the same amount
     * of memory as the target array of element values.
     * The first lookup pointer will point to the start of the elements array.
     * The last lookup pointer will point to the end, just past the last, element of the
     * array. The rest of the lookup pointers point to (max_element * query_multiplier)
     * non-zero element values.
     */
    self->num_lookups = 2 + (size_t)(max_element * self->query_multiplier);

    self->lookups = malloc(self->num_lookups * sizeof(*(self->lookups)));
    if (self->lookups == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }

    ret = fast_search_init_lookups(self, elements, n_elements);
out:
    return ret;
}

int
fast_search_free(fast_search_t *self)
{
    msp_safe_free(self->lookups);
    return 0;
}

static inline const double *
ptr_1st_strict_upper_bound(const double *start, const double *stop, double query)
{
    assert(start <= stop);
    return start + idx_1st_strict_upper_bound(start, (size_t)(stop - start), query);
}

/* PRE-CONDITIONS:
 *   1) query >= 0.0
 * RETURNS:
 *   See idx_1st_strict_upper_bound
 * NOTE:
 *   A non-strict version of this function can be achieved by reimplementing it
 *   to call a non-strict version of `ptr_1st_strict_upper_bound`
 */
size_t
fast_search_idx_strict_upper(fast_search_t *self, double query)
{
    const double **lookups = self->lookups;
    double fidx; // index that can be way larger than max size_t
    size_t idx;
    const double *ptr;

    assert(query >= 0.0);
    /* query_multiplier was set to "rescale" query to the desired lookup table index.
     * Warning: query times query_multiplier might be huge, way bigger than ULLONG_MAX.
     */
    fidx = query * self->query_multiplier;
    if (fidx < self->num_lookups - 1) {
        /* The floating point rescale `fidx` is not huge, it can be truncated now */
        idx = (size_t) fidx;
        /* The query range of [A, B) maps to `idx` where
         * A = idx/query_multiplier and B = (idx + 1)/query_multiplier).
         * The lookup table has been initialized such that
         * `lookup[idx]` and `lookup[idx+1]` points to the first upper bound of A and B
         * respectively in the element values.
         */
        ptr = ptr_1st_strict_upper_bound(lookups[idx], lookups[idx + 1], query);
    } else {
        ptr = lookups[self->num_lookups - 1];
    }
    assert(ptr
           == ptr_1st_strict_upper_bound(
                  lookups[0], lookups[self->num_lookups - 1], query));
    /* function interface is in terms of index to target array of element values */
    return (size_t)(ptr - lookups[0]);
}
