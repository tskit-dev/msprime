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

/* Basic utilities needed in all files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include <tskit/core.h>
#include "util.h"

static const char *
msp_strerror_internal(int err)
{
    const char *ret = "Unknown error";

    switch (err) {
        case 0:
            ret = "Normal exit condition. This is not an error!";
            goto out;
        case MSP_ERR_NO_MEMORY:
            ret = "Out of memory.";
            break;
        case MSP_ERR_GENERIC:
            ret = "Generic error; please file a bug report";
            break;
        case MSP_ERR_BAD_STATE:
            ret = "Bad simulator state. Initialise or reset must be called.";
            break;
        case MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS:
            ret = "Demographic events must be time sorted.";
            break;
        case MSP_ERR_POPULATION_OVERFLOW:
            ret = "Population Overflow occurred.";
            break;
        case MSP_ERR_OUT_OF_BOUNDS:
            ret = "Object reference out of bounds";
            break;
        case MSP_ERR_BAD_PARAM_VALUE:
            ret = "Bad parameter value provided";
            break;
        case MSP_ERR_BAD_POPULATION_CONFIGURATION:
            ret = "Bad population configuration provided.";
            break;
        case MSP_ERR_BAD_POPULATION_SIZE:
            ret = "Bad population size provided. Must be > 0.";
            break;
        case MSP_ERR_POPULATION_OUT_OF_BOUNDS:
            ret = "Population ID out of bounds.";
            break;
        case MSP_ERR_BAD_MIGRATION_MATRIX:
            ret = "Bad migration matrix provided.";
            break;
        case MSP_ERR_BAD_MIGRATION_MATRIX_INDEX:
            ret = "Bad migration matrix index provided.";
            break;
        case MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX:
            ret = "Cannot set diagonal migration matrix elements.";
            break;
        case MSP_ERR_INFINITE_WAITING_TIME:
            ret = "Infinite waiting time until next simulation event.";
            break;
        case MSP_ERR_ASSERTION_FAILED:
            ret = "Internal error; please file a bug report.";
            break;
        case MSP_ERR_SOURCE_DEST_EQUAL:
            ret = "Source and destination populations equal.";
            break;
        case MSP_ERR_BAD_RATE_MAP:
            ret = "Bad rate map provided.";
            break;
        case MSP_ERR_INSUFFICIENT_SAMPLES:
            ret = "At least two samples needed.";
            break;
        case MSP_ERR_BAD_SAMPLES:
            ret = "Bad sample configuration provided.";
            break;
        case MSP_ERR_BAD_MODEL:
            ret = "Model error. Either a bad model, or the requested operation "
                  "is not supported for the current model";
            break;
        case MSP_ERR_INCOMPATIBLE_FROM_TS:
            ret = "The specified tree sequence is not a compatible starting point "
                  "for the current simulation";
            break;
        case MSP_ERR_BAD_START_TIME_FROM_TS:
            ret = "The specified start_time and from_ts are not compatible. All "
                  "node times in the tree sequence must be <= start_time.";
            break;
        case MSP_ERR_BAD_START_TIME:
            ret = "start_time must be >= 0.";
            break;
        case MSP_ERR_BAD_DEMOGRAPHIC_EVENT_TIME:
            ret = "demographic event time must be >= start_time.";
            break;
        case MSP_ERR_TIME_TRAVEL:
            ret = "The simulation model supplied resulted in a parent node having "
                  "a time value <= to its child. This can occur as a result of "
                  "multiple bottlenecks happening at the same time, multiple census "
                  "events at the same time or numerical imprecision with very small"
                  "population sizes.";
            break;
        case MSP_ERR_INTEGRATION_FAILED:
            ret = "GSL numerical integration failed. Please check the stderr for "
                  "details.";
            break;
        case MSP_ERR_BAD_SWEEP_POSITION:
            ret = "Sweep position must be between 0 and sequence length.";
            break;
        case MSP_ERR_BAD_TIME_DELTA:
            ret = "Time delta values must be > 0.";
            break;
        case MSP_ERR_BAD_ALLELE_FREQUENCY:
            ret = "Allele frequency values must be between 0 and 1.";
            break;
        case MSP_ERR_BAD_TRAJECTORY_START_END:
            ret = "Start frequency must be less than end frequency";
            break;
        case MSP_ERR_BAD_SWEEP_GENIC_SELECTION_S:
            ret = "alpha must be > 0";
            break;
        case MSP_ERR_EVENTS_DURING_SWEEP:
            ret = "Demographic and sampling events during a sweep "
                  "are not supported";
            break;
        case MSP_ERR_UNSUPPORTED_OPERATION:
            ret = "Current simulation configuration is not supported.";
            break;
        case MSP_ERR_DTWF_ZERO_POPULATION_SIZE:
            ret = "Population size has decreased to zero individuals.";
            break;
        case MSP_ERR_DTWF_UNSUPPORTED_BOTTLENECK:
            ret = "Bottleneck events are not supported in the DTWF model. "
                  "They can be implemented as population size changes.";
            break;

        case MSP_ERR_BAD_PEDIGREE_NUM_SAMPLES:
            ret = "Not enough sample individuals provided in the pedigree.";
            break;
        case MSP_ERR_EMPTY_PEDIGREE:
            ret = "No individuals in the input pedigree.";
            break;
        case MSP_ERR_OTHER_MODELS_WITH_PED:
            ret = "Cannot combine the through-pedigree simulation with other "
                  "ancestry models";
            break;
        case MSP_ERR_PEDIGREE_IND_NODE_TIME_DISAGREE:
            ret = "The times for the two nodes in a pedigree individual are not equal";
            break;
        case MSP_ERR_PEDIGREE_IND_NODE_POPULATION_DISAGREE:
            ret = "The populations for the two nodes in a pedigree individual "
                  "are not equal";
            break;
        case MSP_ERR_PEDIGREE_TIME_TRAVEL:
            ret = "The time for a parent must be greater than its children";
            break;
        case MSP_ERR_PEDIGREE_IND_NOT_DIPLOID:
            ret = "All individuals in the input pedigree must be associated with "
                  "exactly two nodes";
            break;
        case MSP_ERR_PEDIGREE_IND_NOT_TWO_PARENTS:
            ret = "All individuals in the input pedigree must be associated with "
                  "exactly two parents (can be TSK_NULL, if not known)";
            break;
        case MSP_ERR_PEDIGREE_INTERNAL_SAMPLE:
            ret = "Samples that are internal nodes in the pedigree are not "
                  "currently supported. Please comment on this GitHub issue if you "
                  "would like to see this feature implemented: "
                  "https://github.com/tskit-dev/msprime/issues/1855 ";
            break;

        case MSP_ERR_BAD_PROPORTION:
            ret = "Proportion values must have 0 <= x <= 1";
            break;
        case MSP_ERR_BAD_BETA_MODEL_ALPHA:
            ret = "Bad alpha. Must have 1 < alpha < 2";
            break;
        case MSP_ERR_BAD_TRUNCATION_POINT:
            ret = "Bad truncation_point. Must have 0 < truncation_point.";
            break;
        case MSP_ERR_BAD_RATE_VALUE:
            ret = "Rates must be non-negative and finite";
            break;
        case MSP_ERR_INCOMPATIBLE_MUTATION_MAP_LENGTH:
            ret = "Incompatible mutation map: sequence length differs from that of "
                  "the tree sequence.";
            break;
        case MSP_ERR_INSUFFICIENT_INTERVALS:
            ret = "At least one interval must be specified.";
            break;
        case MSP_ERR_INTERVAL_MAP_START_NON_ZERO:
            ret = "The first interval must start with zero";
            break;
        case MSP_ERR_INTERVAL_POSITIONS_UNSORTED:
            ret = "Interval positions must be listed in increasing order";
            break;
        case MSP_ERR_NONFINITE_INTERVAL_POSITION:
            ret = "Interval positions must be finite.";
            break;
        case MSP_ERR_BAD_C:
            ret = "Bad C. Must have 0 < C ";
            break;
        case MSP_ERR_BAD_PSI:
            ret = "Bad PSI. Must have 0 < PSI <= 1";
            break;
        case MSP_ERR_UNKNOWN_ALLELE:
            ret = "Existing allele(s) incompatible with mutation model alphabet.";
            break;
        case MSP_ERR_INSUFFICIENT_ALLELES:
            ret = "Must have at least two alleles.";
            break;
        case MSP_ERR_BAD_ROOT_PROBABILITIES:
            ret = "Root probabilities must be nonnegative and sum to one.";
            break;
        case MSP_ERR_BAD_TRANSITION_MATRIX:
            ret = "Each row of the transition matrix must be nonnegative and sum to "
                  "one.";
            break;
        case MSP_ERR_BAD_SLIM_PARAMETERS:
            ret = "SLiM mutation IDs and mutation type IDs must be nonnegative.";
            break;
        case MSP_ERR_MUTATION_ID_OVERFLOW:
            ret = "Mutation ID overflow.";
            break;
        case MSP_ERR_BREAKPOINT_MASS_NON_FINITE:
            ret = "An unlikely numerical error occurred computing recombination "
                  "breakpoints (non finite breakpoint mass). Please check your "
                  "parameters, and if they make sense help us fix the problem "
                  "by opening an issue on GitHub.";
            break;
        case MSP_ERR_BREAKPOINT_RESAMPLE_OVERFLOW:
            ret = "An unlikely numerical error occurred computing recombination "
                  "breakpoints (resample overflow). Please check your "
                  "parameters, and if they make sense help us fix the problem "
                  "by opening an issue on GitHub.";
            break;
        case MSP_ERR_TRACTLEN_RESAMPLE_OVERFLOW:
            ret = "An unlikely numerical error occurred computing gene conversion"
                  "tract lengths (resample overflow). Please check your "
                  "parameters, and if they make sense help us fix the problem "
                  "by opening an issue on GitHub.";
            break;
        case MSP_ERR_FENWICK_REBUILD_FAILED:
            ret = "An unlikely numerical error occurred (Fenwick tree rebuild "
                  "did not reduce drift sufficiently). Please check your "
                  "parameters, and if they make sense help us fix the problem "
                  "by opening an issue on GitHub.";
            break;
        case MSP_ERR_BAD_PLOIDY:
            ret = "Ploidy must be at least 1";
            break;
        case MSP_ERR_DTWF_MIGRATION_MATRIX_NOT_STOCHASTIC:
            ret = "The row sums of the migration matrix must not exceed one for "
                  "the discrete time Wright-Fisher model.";
            break;
        case MSP_ERR_DTWF_GC_NOT_SUPPORTED:
            ret = "Gene conversion is not supported in the DTWF model";
            break;
        case MSP_ERR_SWEEPS_GC_NOT_SUPPORTED:
            ret = "Gene conversion is not supported in the selective sweep model";
            break;
        case MSP_ERR_BAD_SEQUENCE_LENGTH:
            ret = "Sequence length must be > 0";
            break;
        case MSP_ERR_ZERO_POPULATIONS:
            ret = "At least one population must be defined";
            break;
        case MSP_ERR_BAD_ANCIENT_SAMPLE_NODE:
            ret = "Only isolated sample nodes are supported as ancient samples";
            break;
        case MSP_ERR_UNKNOWN_TIME_NOT_SUPPORTED:
            ret = "Kept mutations must have known mutation times (which can be added "
                  "using compute_mutation_times).";
            break;
        case MSP_ERR_DTWF_DIPLOID_ONLY:
            ret = "The DTWF model only supports ploidy = 2";
            break;
        case MSP_ERR_TOO_MANY_EVENT_POPULATIONS:
            ret = "Cannot have more than 100 populations in one event. "
                  "If this is something that you need to do, please open an issue "
                  "on GitHub";
            break;
        case MSP_ERR_DUPLICATE_POPULATION:
            ret = "Population IDs must be unique";
            break;
        case MSP_ERR_POPULATION_INACTIVE_MOVE:
            ret = "Attempt to move a lineage into an inactive population";
            break;
        case MSP_ERR_POPULATION_INACTIVE_SAMPLE:
            ret = "Attempt to sample a lineage from an inactive population";
            break;
        case MSP_ERR_POPULATION_PREVIOUSLY_ACTIVE:
            ret = "Attempt to set a previously active population to active";
            break;
        case MSP_ERR_POPULATION_CURRENTLY_ACTIVE:
            ret = "Attempt to set a currently active population to active";
            break;
        case MSP_ERR_DEACTIVATE_INACTIVE_POPULATION:
            ret = "Attempt to set a currently inactive population to inactive";
            break;
        case MSP_ERR_DEACTIVATE_PREVIOUSLY_ACTIVE_POPULATION:
            ret = "Attempt to set a previously active population to inactive";
            break;
        case MSP_ERR_ADMIX_DERIVED_NOT_ACTIVE:
            ret = "The derived population in an admixture must be active";
            break;
        case MSP_ERR_SPLIT_DERIVED_NOT_ACTIVE:
            ret = "The derived population in a population split must be active";
            break;
        case MSP_ERR_POP_SIZE_ZERO_SAMPLE:
            ret = "Attempt to sample lineage in a population with size=0";
            break;
        default:
            ret = "Error occurred generating error string. Please file a bug "
                  "report!";
            break;
    }
out:
    return ret;
}

int
msp_set_tsk_error(int err)
{
    /* Flip this bit. As the error is negative, this sets the bit to 0 */
    return err ^ (1 << MSP_TSK_ERR_BIT);
}

bool
msp_is_tsk_error(int err)
{
    return (err < 0) && !(err & (1 << MSP_TSK_ERR_BIT));
}

const char *
msp_strerror(int err)
{
    if (msp_is_tsk_error(err)) {
        err ^= (1 << MSP_TSK_ERR_BIT);
        return tsk_strerror(err);
    } else {
        return msp_strerror_internal(err);
    }
}

void
__msp_safe_free(void **ptr)
{
    if (ptr != NULL) {
        if (*ptr != NULL) {
            free(*ptr);
            *ptr = NULL;
        }
    }
}

bool
doubles_almost_equal(double a, double b, double eps)
{
    return (fabs(a) < eps && fabs(b) < eps) || gsl_fcmp(a, b, eps) == 0;
}

/* Inputs:
 *   - `lengths` is an array of the real number lengths of `num_intervals`
 *      consecutive half-open intervals on the real number line starting at zero,
 *      i.e. given lengths L_i, the half-open intervals are:
 *
 *           [S_0, S_1), [S_1, S_2), [S_2, S_3), ..., [S_{N-1}, S_N)
 *
 *      where S_0 = 0.0 and S_{i+1} = S_i + L_i and N = `num_intervals`.
 *
 *   - `x` is a number to find in the following extension of those intervals:
 *
 *     [-INFINITY, S_1), [S_1, S_2), [S_2, S_3), ..., [S_{N-1}, S_N), [S_N, INFINITY]
 *
 * Returns:
 *   - the index of the interval in which `x` was found (0 to `num_interval`, inclusive)
 *   - `num_interval` if `x` is NaN
 * Pre-conditions:
 *   - `lengths` must point to `num_interval` non-negative non-NaN doubles
 * Post-condition:
 *   - returned value is less than or equal to `num_intervals` (and non-negative)
 *
 */
static inline size_t
positive_interval_select(double x, size_t num_intervals, const double *lengths)
{
    size_t i = 0;
    double sum = 0.0;
    while (i < num_intervals) {
        assert(lengths[i] >= 0.0);
        sum += lengths[i];
        if (x < sum) {
            break;
        }
        i++;
    }
    return i;
}

/* Convenient helper function to use with random variant `u` and array of probabilities
 * that approximately sum to one. Note the function does not actually use the last
 * probability because its correct value is implied by all the previous probabilities.
 */
size_t
probability_list_select(double u, size_t num_probs, double const *probs)
{
    return (num_probs > 0 ? positive_interval_select(u, num_probs - 1, probs) : 0);
}

/* binary search functions */

/* This function follows standard semantics of:
 *   numpy.searchsorted(..., side='left') and
 *   std::lower_bound() from the standard C++ <algorithm> library
 * PRE-CONDITION:
 *   1) `values` are sorted
 * RETURNS:
 *   First (leftmost) `index` of upper bounds of `query`
 *   **or** `n_values` if all `values` are strict lower bounds of `query`
 * POST-CONDITION:
 *   values[index-1] < query <= values[index]
 *   (ignoring comparisons with invalid [] indexing)
 */
size_t
idx_1st_upper_bound(const double *values, size_t n_values, double query)
{
    size_t l = 0;
    size_t r = n_values;
    size_t m;

    while (l < r) {
        m = (l + r) / 2UL;
        assert(values[l] <= values[m]);
        if (values[m] < query) {
            l = m + 1;
        } else {
            r = m;
        }
    }
    return l;
}

/* This function follows standard semantics of:
 *   std::upper_bound() from the standard C++ <algorithm> library
 *   and numpy.searchsorted(..., side='right') [Caveat]
 * PRE-CONDITION:
 *   1) `values` are sorted
 * RETURNS:
 *   First (leftmost) `index` of strict upper bounds of `query`
 *   **or** `n_values` if all `values` are lower bounds of `query`
 * POST-CONDITION:
 *    values[index-1] <= query < values[index]
 *    (ignoring comparisons with invalid [] indexing)
 * [Caveat]:
 *   This function matches NaN semantics of std::upper_bound, not numpy
 */
size_t
idx_1st_strict_upper_bound(const double *elements, size_t n_elements, double query)
{
    return sub_idx_1st_strict_upper_bound(elements, 0, n_elements, query);
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
 * required.  This is achieved through an lookup array of element indexes.
 * Although there need not be a linear mapping between element values and the
 * range of query values, the lookup values (element indexes) are linear with
 * the range of query values.
 *
 * The `query_multiplier` determines the linear mapping between possible query
 * values and indexes in the lookup array of element indexes. Because the
 * `query_multiplier` is a power of 2, a predictable range of possible "query"
 * values will truncate to a specific lookup index.
 *
 * The size of the lookup array could be adjustable, but this simple
 * implementation chooses one simple strategy for the size of the lookup array.
 * For a perfectly linear sequence of element values, the lookup array will be
 * a simple progression of element indexes.
 * The extremely simple case of element values being the non-negative integers
 * is coded in the `test_fast_search_identity` unit test.
 * An example of alternative linear sequences with different `query_multiplier`
 * are coded in `test_fast_search_2powers`.
 *
 * Consider example target array of elements:
 *
 * (element index, element value)
 * [0], 0.0
 * [1], 0.1
 * [2], 0.8
 * [3], 2.7
 * [4], 6.4
 * [5], ___ (the end, past last element of the array)
 *
 * The power of 2 greater or equal to the maximum element is 8.
 * The number of memory address steps to the maximum element is 4.
 * The power of 2 greater or equal to that number of steps is also 4.
 * Thus the query multiplier is going to be 0.5.
 * The lookup table is going to start with an element index to the last of the elements.
 * The last element index in the lookup table will point to the end of the element array.
 * The lookup table should be just the right size so that rescaling (and truncating) a
 * query equaling the maximum element will result in an index to the last element index
 * in the lookup table. This is `trunc(6.4 * 0.5) = 3`. Thus the lookup table
 * will be `(1 + 1 + 3) = 5` element indexes.
 *
 * The lookup table of element indexes will handle the following 5 query ranges:
 *
 * [0.0, 2.0)
 * [2.0, 4.0)
 * [4.0, 6.0)
 * [6.0, 8.0)
 * [8.0, INFINITY)
 *
 * with each lookup entry pointing to the first upper bound to the minimum in those
 * query ranges. Thus we get the following lookup table of element indexes:
 *
 * (lookup index, element index)
 * [0] 0 (0.0 in elements is upper bound to query range minimum 0.0)
 * [1] 3 (2.7 in elements is upper bound to query range minimum 2.0)
 * [2] 4 (6.4 in elements is upper bound to query range minimum 4.0)
 * [4] 4 (6.4 in elements is upper bound to query range minimum 6.0)
 * [5] 5 (no upper bound in elements to query range minimum 8.0)
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
    unsigned idx;
    double min_query;

    if (n_elements > UINT_MAX) {
        ret = MSP_ERR_BAD_PARAM_VALUE; // LCOV_EXCL_LINE
        goto out;                      // LCOV_EXCL_LINE
    }
    self->elements = elements;
    self->lookups[0] = 0;
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
        self->lookups[idx] = (unsigned) (ptr - elements);
    }
out:
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
    const uint64_t max_input_size = 1ULL << (DBL_MANT_DIG - 1); // 4096 terabytes

    memset(self, 0, sizeof(*self));

    if (!valid_sorted_nonempty_array(elements, n_elements)) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    if (elements[0] != 0.0 || (uint64_t) n_elements >= max_input_size) {
        ret = MSP_ERR_BAD_PARAM_VALUE;
        goto out;
    }
    max_element = elements[n_elements - 1];
    if (!isfinite(max_element)) {
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
     * The first lookup element index will point to the start of the elements array.
     * The last lookup element index will point to the end, just past the last, element
     * of the array. The rest of the lookup element indexes point to (max_element *
     * query_multiplier) non-zero element values.
     */
    self->num_lookups = 2 + (size_t)(max_element * self->query_multiplier);

    self->query_cutoff = ((double) self->num_lookups - 1) / self->query_multiplier;

    self->lookups = malloc(self->num_lookups * sizeof(*(self->lookups)));
    if (self->lookups == NULL) {
        ret = MSP_ERR_NO_MEMORY; // LCOV_EXCL_LINE
        goto out;                // LCOV_EXCL_LINE
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

/*******************************
 *  `extern inline` declarations
 *  Due to compiler/linker limitations of C99, `inline` function declarations
 *  must be restated with the `extern` keyword in one ".c" file
 *  (or never declared in header files and declared with `static`)
 */

extern inline size_t sub_idx_1st_strict_upper_bound(
    const double *base, size_t start, size_t stop, double query);
extern inline size_t fast_search_idx_strict_upper(fast_search_t *self, double query);

/*   The gsl_ran_flat() function is supposed to output lo<=x<hi, but
 *   sometimes (with probability <1e-9) we find x=hi.
 *   This function includes a while loop to ensure x<hi.
 *   https://github.com/tskit-dev/msprime/issues/1997
 */
double
msp_gsl_ran_flat(gsl_rng *rng, double lo, double hi)
{
    double x;

    /* The logic here requires that lo <= hi, but */
    /* gsl_ran_flat produces output in this case. Keep */
    /* this assert to make sure that msprime's usage */
    /* has the required properties.  */
    tsk_bug_assert(lo <= hi);
    if (lo == hi) {
        x = lo;
    } else {
        do {
            x = gsl_ran_flat(rng, lo, hi);
        } while (x >= hi);
    }
    return x;
}
