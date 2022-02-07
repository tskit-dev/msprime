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
#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include <gsl/gsl_randist.h>

#include "tskit.h"

#ifdef __GNUC__
/*
 * raise a compiler warning if a potentially error raising function's return
 * value is not used.
 */
#define MSP_WARN_UNUSED __attribute__((warn_unused_result))
/* Annotate a function parameter as unused */
#define MSP_UNUSED(x) MSP_UNUSED_##x __attribute__((__unused__))
#else
#define MSP_WARN_UNUSED
#define MSP_UNUSED(x) MSP_UNUSED_##x
/* Don't bother with restrict for MSVC */
#define restrict
#endif

/* clang-format off */
/* Error codes */
#define MSP_ERR_GENERIC                                             -1
#define MSP_ERR_NO_MEMORY                                           -2
#define MSP_ERR_BAD_STATE                                           -3
#define MSP_ERR_BAD_PARAM_VALUE                                     -4
#define MSP_ERR_OUT_OF_BOUNDS                                       -5
#define MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS                         -6
#define MSP_ERR_POPULATION_OVERFLOW                                 -7
#define MSP_ERR_POPULATION_OUT_OF_BOUNDS                            -8
#define MSP_ERR_BAD_POPULATION_CONFIGURATION                        -9
#define MSP_ERR_BAD_MIGRATION_MATRIX                                -10
#define MSP_ERR_BAD_MIGRATION_MATRIX_INDEX                          -11
#define MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX                     -12
#define MSP_ERR_INFINITE_WAITING_TIME                               -13
#define MSP_ERR_ASSERTION_FAILED                                    -14
#define MSP_ERR_SOURCE_DEST_EQUAL                                   -15
#define MSP_ERR_BAD_RATE_MAP                                        -16
#define MSP_ERR_BAD_POPULATION_SIZE                                 -17
#define MSP_ERR_BAD_SAMPLES                                         -18
#define MSP_ERR_BAD_MODEL                                           -19
#define MSP_ERR_INSUFFICIENT_SAMPLES                                -20
#define MSP_ERR_DUPLICATE_SITE_POSITION                             -21
#define MSP_ERR_UNDEFINED_MULTIPLE_MERGER_COALESCENT                -22
#define MSP_ERR_INCOMPATIBLE_FROM_TS                                -23
#define MSP_ERR_BAD_START_TIME_FROM_TS                              -24
#define MSP_ERR_BAD_START_TIME                                      -25
#define MSP_ERR_BAD_DEMOGRAPHIC_EVENT_TIME                          -26
#define MSP_ERR_TIME_TRAVEL                                         -27
#define MSP_ERR_INTEGRATION_FAILED                                  -28
#define MSP_ERR_BAD_SWEEP_POSITION                                  -29
#define MSP_ERR_BAD_TIME_DELTA                                      -30
#define MSP_ERR_BAD_ALLELE_FREQUENCY                                -31
#define MSP_ERR_BAD_TRAJECTORY_START_END                            -32
#define MSP_ERR_BAD_SWEEP_GENIC_SELECTION_S                         -33
#define MSP_ERR_EVENTS_DURING_SWEEP                                 -34
#define MSP_ERR_UNSUPPORTED_OPERATION                               -35
#define MSP_ERR_DTWF_ZERO_POPULATION_SIZE                           -36
#define MSP_ERR_DTWF_UNSUPPORTED_BOTTLENECK                         -37
#define MSP_ERR_BAD_PROPORTION                                      -38
/* REUSE 39 and 40 */
#define MSP_ERR_BAD_BETA_MODEL_ALPHA                                -41
#define MSP_ERR_BAD_TRUNCATION_POINT                                -42
#define MSP_ERR_BAD_RATE_VALUE                                      -43
#define MSP_ERR_INCOMPATIBLE_MUTATION_MAP_LENGTH                    -44
#define MSP_ERR_INSUFFICIENT_INTERVALS                              -45
#define MSP_ERR_INTERVAL_MAP_START_NON_ZERO                         -46
#define MSP_ERR_NONFINITE_INTERVAL_POSITION                         -47
#define MSP_ERR_INTERVAL_POSITIONS_UNSORTED                         -48
#define MSP_ERR_BAD_C                                               -49
#define MSP_ERR_BAD_PSI                                             -50
#define MSP_ERR_UNKNOWN_ALLELE                                      -51
#define MSP_ERR_INSUFFICIENT_ALLELES                                -52
#define MSP_ERR_BAD_ROOT_PROBABILITIES                              -53
#define MSP_ERR_BAD_TRANSITION_MATRIX                               -54
#define MSP_ERR_BAD_SLIM_PARAMETERS                                 -56
#define MSP_ERR_MUTATION_ID_OVERFLOW                                -57
#define MSP_ERR_BREAKPOINT_MASS_NON_FINITE                          -58
#define MSP_ERR_BREAKPOINT_RESAMPLE_OVERFLOW                        -59
#define MSP_ERR_TRACTLEN_RESAMPLE_OVERFLOW                          -60
#define MSP_ERR_FENWICK_REBUILD_FAILED                              -61
#define MSP_ERR_BAD_PLOIDY                                          -62
#define MSP_ERR_DTWF_MIGRATION_MATRIX_NOT_STOCHASTIC                -63
#define MSP_ERR_DTWF_GC_NOT_SUPPORTED                               -64
#define MSP_ERR_SWEEPS_GC_NOT_SUPPORTED                             -65
#define MSP_ERR_BAD_SEQUENCE_LENGTH                                 -66
#define MSP_ERR_ZERO_POPULATIONS                                    -67
#define MSP_ERR_BAD_ANCIENT_SAMPLE_NODE                             -68
#define MSP_ERR_UNKNOWN_TIME_NOT_SUPPORTED                          -69
#define MSP_ERR_DTWF_DIPLOID_ONLY                                   -70
#define MSP_ERR_BAD_ANCESTRAL_MUTATION                              -71
#define MSP_ERR_TOO_MANY_EVENT_POPULATIONS                          -72
#define MSP_ERR_DUPLICATE_POPULATION                                -73
#define MSP_ERR_POPULATION_INACTIVE_MOVE                            -74
#define MSP_ERR_POPULATION_INACTIVE_SAMPLE                          -75
#define MSP_ERR_POPULATION_PREVIOUSLY_ACTIVE                        -76
#define MSP_ERR_SPLIT_DERIVED_NOT_ACTIVE                            -77
#define MSP_ERR_ADMIX_DERIVED_NOT_ACTIVE                            -78
#define MSP_ERR_POP_SIZE_ZERO_SAMPLE                                -79
#define MSP_ERR_POPULATION_CURRENTLY_ACTIVE                         -80
#define MSP_ERR_DEACTIVATE_INACTIVE_POPULATION                      -81
#define MSP_ERR_DEACTIVATE_PREVIOUSLY_ACTIVE_POPULATION             -82
#define MSP_ERR_BAD_PEDIGREE_NUM_SAMPLES                            -83
#define MSP_ERR_OTHER_MODELS_WITH_PED                               -84
#define MSP_ERR_EMPTY_PEDIGREE                                      -85
#define MSP_ERR_PEDIGREE_IND_NODE_TIME_DISAGREE                     -86
#define MSP_ERR_PEDIGREE_IND_NODE_POPULATION_DISAGREE               -87
#define MSP_ERR_PEDIGREE_TIME_TRAVEL                                -88
#define MSP_ERR_PEDIGREE_IND_NOT_DIPLOID                            -89
#define MSP_ERR_PEDIGREE_IND_NOT_TWO_PARENTS                        -90
#define MSP_ERR_PEDIGREE_INTERNAL_SAMPLE                            -91

/* clang-format on */
/* This bit is 0 for any errors originating from tskit */
#define MSP_TSK_ERR_BIT 13

int msp_set_tsk_error(int err);
bool msp_is_tsk_error(int err);
const char *msp_strerror(int err);
void __msp_safe_free(void **ptr);

#define msp_safe_free(pointer) __msp_safe_free((void **) &(pointer))

bool doubles_almost_equal(double a, double b, double eps);

size_t probability_list_select(double u, size_t num_probs, double const *probs);

/* binary search functions */

size_t idx_1st_upper_bound(const double *values, size_t n_values, double query);
size_t idx_1st_strict_upper_bound(
    const double *elements, size_t n_elements, double query);

typedef struct {
    const double *elements;
    double query_multiplier;
    double query_cutoff;
    size_t num_lookups;
    unsigned *lookups;
} fast_search_t;

int fast_search_alloc(fast_search_t *self, const double *values, size_t n_values);
int fast_search_free(fast_search_t *self);
inline size_t fast_search_idx_strict_upper(fast_search_t *self, double query);

double msp_gsl_ran_flat(gsl_rng *rng, double lo, double hi);

/***********************************
 * INLINE FUNCTION IMPLEMENTATIONS *
 ***********************************/

inline size_t
sub_idx_1st_strict_upper_bound(
    const double *base, size_t start, size_t stop, double query)
{
    while (start < stop) {
        size_t mid = (start + stop) / 2;
        assert(base[start] <= base[mid]);
        if (!(base[mid] > query)) { // match NaN logic of std::upper_bound
            start = mid + 1;
        } else {
            stop = mid;
        }
    }
    return stop;
}

/* PRE-CONDITIONS:
 *   1) query >= 0.0
 * RETURNS:
 *   See idx_1st_strict_upper_bound
 * NOTE:
 *   A non-strict version of this function can be achieved by reimplementing it
 *   to call a non-strict version of `idx_1st_strict_upper_bound`
 */
inline size_t
fast_search_idx_strict_upper(fast_search_t *self, double query)
{
    size_t ret;
    unsigned *lookups = self->lookups;
    assert(query >= 0.0);
    if (query < self->query_cutoff) {
        int64_t idx = (int64_t)(query * self->query_multiplier);
        /* The query range of [A, B) maps to `idx` where
         * A = idx/query_multiplier and B = (idx + 1)/query_multiplier).
         * The lookup table has been initialized such that
         * `lookup[idx]` and `lookup[idx+1]` points to the first upper bound of A and B
         * respectively in the element values.
         */
        ret = sub_idx_1st_strict_upper_bound(
            self->elements, lookups[idx], lookups[idx + 1], query);
    } else {
        ret = lookups[self->num_lookups - 1];
    }
    assert(ret
           == idx_1st_strict_upper_bound(
                  self->elements, lookups[self->num_lookups - 1], query));
    /* function interface is in terms of index to target array of element values */
    return ret;
}

#endif /*__UTIL_H__*/
