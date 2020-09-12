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
#define MSP_ERR_BAD_SWEEP_GENIC_SELECTION_ALPHA                     -33
#define MSP_ERR_EVENTS_DURING_SWEEP                                 -34
#define MSP_ERR_UNSUPPORTED_OPERATION                               -35
#define MSP_ERR_DTWF_ZERO_POPULATION_SIZE                           -36
#define MSP_ERR_DTWF_UNSUPPORTED_BOTTLENECK                         -37
#define MSP_ERR_BAD_PROPORTION                                      -38
#define MSP_ERR_BAD_PEDIGREE_NUM_SAMPLES                            -39
#define MSP_ERR_BAD_PEDIGREE_ID                                     -40
#define MSP_ERR_BAD_BETA_MODEL_ALPHA                                -41
#define MSP_ERR_BAD_TRUNCATION_POINT                                -42
#define MSP_ERR_BAD_RATE_VALUE                                      -43
#define MSP_ERR_INCOMPATIBLE_MUTATION_MAP                           -44
#define MSP_ERR_INSUFFICIENT_INTERVALS                              -45
#define MSP_ERR_INTERVAL_MAP_START_NON_ZERO                         -46
#define MSP_ERR_NONFINITE_INTERVAL_POSITION                         -47
#define MSP_ERR_INTERVAL_POSITIONS_UNSORTED                         -48
#define MSP_ERR_BAD_C                                               -49
#define MSP_ERR_BAD_PSI                                             -50
#define MSP_ERR_UNKNOWN_ALLELE                                      -51
#define MSP_ERR_MUTATION_GENERATION_OUT_OF_ORDER                    -52
#define MSP_ERR_INSUFFICIENT_ALLELES                                -53
#define MSP_ERR_BAD_ROOT_PROBABILITIES                              -54
#define MSP_ERR_BAD_TRANSITION_MATRIX                               -55
#define MSP_ERR_BAD_SLIM_PARAMETERS                                 -57
#define MSP_ERR_MUTATION_ID_OVERFLOW                                -58
#define MSP_ERR_BREAKPOINT_MASS_NON_FINITE                          -59
#define MSP_ERR_BREAKPOINT_RESAMPLE_OVERFLOW                        -60
#define MSP_ERR_TRACKLEN_RESAMPLE_OVERFLOW                          -61
#define MSP_ERR_FENWICK_REBUILD_FAILED                              -62
#define MSP_ERR_BAD_PLOIDY                                          -63
#define MSP_ERR_DTWF_MIGRATION_MATRIX_NOT_STOCHASTIC                -64
#define MSP_ERR_DTWF_GC_NOT_SUPPORTED                               -65
#define MSP_ERR_SWEEPS_GC_NOT_SUPPORTED                             -66
#define MSP_ERR_BAD_SEQUENCE_LENGTH                                 -67
#define MSP_ERR_ZERO_POPULATIONS                                    -68
#define MSP_ERR_BAD_ANCIENT_SAMPLE_NODE                             -69

/* clang-format on */
/* This bit is 0 for any errors originating from tskit */
#define MSP_TSK_ERR_BIT 13

int msp_set_tsk_error(int err);
bool msp_is_tsk_error(int err);
const char *msp_strerror(int err);
void __msp_safe_free(void **ptr);

#define msp_safe_free(pointer) __msp_safe_free((void **) &(pointer))

size_t msp_binary_interval_search(double query, const double *values, size_t n_values);
bool doubles_almost_equal(double a, double b, double eps);

size_t probability_list_select(double u, size_t num_probs, double const *probs);

#endif /*__UTIL_H__*/
