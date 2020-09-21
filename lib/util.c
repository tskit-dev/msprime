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
#include <assert.h>

#include <gsl/gsl_math.h>

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
                  "a time value <= to its child. This can occur either as a result "
                  "of multiple bottlenecks happening at the same time or because of "
                  "numerical imprecision with very small population sizes.";
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
        case MSP_ERR_BAD_SWEEP_GENIC_SELECTION_ALPHA:
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
            ret = "The number of haploid lineages denoted by sample_size must "
                  "be divisible by ploidy (default 2)";
            break;
        case MSP_ERR_BAD_PEDIGREE_ID:
            ret = "Individual IDs in pedigrees must be strictly > 0.";
            break;
        case MSP_ERR_BAD_PROPORTION:
            ret = "Proportion values must have 0 <= x <= 1";
            break;
        case MSP_ERR_BAD_BETA_MODEL_ALPHA:
            ret = "Bad alpha. Must have 1 < alpha < 2";
            break;
        case MSP_ERR_BAD_TRUNCATION_POINT:
            ret = "Bad truncation_point. Must have 0 < truncation_point <= 1";
            break;
        case MSP_ERR_BAD_RATE_VALUE:
            ret = "Rates must be non-negative and finite";
            break;
        case MSP_ERR_INCOMPATIBLE_MUTATION_MAP:
            ret = "Mutation map is not compatible with specified tables.";
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
        case MSP_ERR_MUTATION_GENERATION_OUT_OF_ORDER:
            ret = "Tree sequence contains mutations that would descend from "
                  "existing mutations: finite site mutations must be generated on "
                  "older time periods first.";
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
            ret = "An unlikely numerical error occured computing recombination "
                  "breakpoints (non finite breakpoint mass). Please check your "
                  "parameters, and if they make sense help us fix the problem "
                  "by opening an issue on GitHub.";
            break;
        case MSP_ERR_BREAKPOINT_RESAMPLE_OVERFLOW:
            ret = "An unlikely numerical error occured computing recombination "
                  "breakpoints (resample overflow). Please check your "
                  "parameters, and if they make sense help us fix the problem "
                  "by opening an issue on GitHub.";
            break;
        case MSP_ERR_TRACKLEN_RESAMPLE_OVERFLOW:
            ret = "An unlikely numerical error occured computing gene conversion"
                  "track lengths (resample overflow). Please check your "
                  "parameters, and if they make sense help us fix the problem "
                  "by opening an issue on GitHub.";
            break;
        case MSP_ERR_FENWICK_REBUILD_FAILED:
            ret = "An unlikely numerical error occured (Fenwick tree rebuild "
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

/* Find the `index` of the interval within `values` the `query` fits, such that
 * values[index-1] < query <= values[index]
 * Will find the leftmost such index
 * Assumes `values` are sorted
 */
size_t
msp_binary_interval_search(double query, const double *values, size_t n_values)
{
    if (n_values == 0) {
        return 0;
    }
    size_t l = 0;
    size_t r = n_values - 1;
    size_t m;

    while (l < r) {
        m = (l + r) / 2UL;

        if (values[m] < query) {
            l = m + 1;
        } else {
            r = m;
        }
    }
    return l;
}

bool
doubles_almost_equal(double a, double b, double eps)
{
    return (fabs(a) < eps && fabs(b) < eps) || gsl_fcmp(a, b, eps) == 0;
}

/* Inputs:
 *   - `lengths` is an array of the real number lengths of `num_intervals`
 *      consecutive half-open intervals on the real number line start at zero,
 *      i.e. given lengths L_i, the intervals are:
 *          [0.0, L_0), [L_0, L_0 + L_1), [L_0 + L_1, L_0 + L_1 + L_2), ...
 *   - `x` is a real number to find in one of the intervals or after all of them
 * Returns:
 *   - the index (starting at zero) of the interval in which `x` was found
 *   - `num_intervals` if `x` is at or past the end of the last interval
 *     (the sum of all the lengths)
 *   - zero if `x` is negative
 *   - `num_interval` if `x` is NaN
 * Pre-conditions:
 *   - `lengths` must point to `num_interval` non-negative non-NaN doubles
 * Post-condition:
 *   - returned value is less than or equal to `num_intervals` (and non-negative)
 *
 */
static inline size_t
positive_interval_select(double x, size_t num_intervals, double const *restrict lengths)
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
