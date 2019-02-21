/*
** Copyright (C) 2015-2018 University of Oxford
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
	#define MSP_WARN_UNUSED __attribute__ ((warn_unused_result))
    /* Annotate a function parameter as unused */
	#define MSP_UNUSED(x) MSP_UNUSED_ ## x __attribute__((__unused__))
#else
	#define MSP_WARN_UNUSED
	#define MSP_UNUSED(x) MSP_UNUSED_ ## x
#endif

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
#define MSP_ERR_BAD_RECOMBINATION_MAP                               -16
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
#define MSP_ERR_RECOMB_MAP_TOO_COARSE                               -27
#define MSP_ERR_TIME_TRAVEL                                         -28
#define MSP_ERR_INTEGRATION_FAILED                                  -29
#define MSP_ERR_BAD_SWEEP_POSITION                                  -30
#define MSP_ERR_BAD_TIME_DELTA                                      -31
#define MSP_ERR_BAD_ALLELE_FREQUENCY                                -32
#define MSP_ERR_BAD_TRAJECTORY_START_END                            -33
#define MSP_ERR_EVENTS_DURING_SWEEP                                 -34
#define MSP_ERR_UNSUPPORTED_OPERATION                               -35
#define MSP_ERR_DTWF_ZERO_POPULATION_SIZE                           -36
#define MSP_ERR_DTWF_UNSUPPORTED_BOTTLENECK                         -37
#define MSP_ERR_BAD_PROPORTION                                      -38
#define MSP_ERR_BAD_PEDIGREE_NUM_SAMPLES                            -39
#define MSP_ERR_BAD_PEDIGREE_ID                                     -40
#define MSP_ERR_BAD_ALPHA                                           -41
#define MSP_ERR_BAD_TRUNCATION_POINT                                -42
#define MSP_ERR_BAD_MUTATION_MAP_SIZE                               -43
#define MSP_ERR_BAD_MUTATION_MAP_POSITION                           -44
#define MSP_ERR_BAD_MUTATION_MAP_RATE                               -45
#define MSP_ERR_INCOMPATIBLE_MUTATION_MAP                           -46

/* This bit is 0 for any errors originating from tskit */
#define MSP_TSK_ERR_BIT 13

int msp_set_tsk_error(int err);
bool msp_is_tsk_error(int err);
const char * msp_strerror(int err);
void __msp_safe_free(void **ptr);

#define msp_safe_free(pointer) __msp_safe_free((void **) &(pointer))

size_t msp_binary_interval_search(double query, double *values, size_t n_values);

#endif /*__UTIL_H__*/
