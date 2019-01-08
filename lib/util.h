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

/*
 * raise a compiler warning if a potentially error raising function's return
 * value is not used.
 */
#ifdef __GNUC__
    #define WARN_UNUSED __attribute__ ((warn_unused_result))
#else
    #define WARN_UNUSED
    /* Don't bother with restrict for MSVC */
    #define restrict
#endif

// BCH 30 May 2018: probably the WARN_UNUSED macro above should be retired
// in favor of MSP_WARN_UNUSED (adapted from KAS_WARN_UNUSED), for clean
// naming, but I didn't make that change since it would be extensive
#ifdef __GNUC__
	#define MSP_WARN_UNUSED __attribute__ ((warn_unused_result))
	#define MSP_UNUSED(x) MSP_UNUSED_ ## x __attribute__((__unused__))
#else
	#define MSP_WARN_UNUSED
	#define MSP_UNUSED(x) MSP_UNUSED_ ## x
#endif

/* This sets up MSP_DBL_DECIMAL_DIG, which can then be used as a
 * precision specifier when writing out doubles, if you want sufficient
 * decimal digits to be written to guarantee a lossless round-trip
 * after being read back in.  Usage:
 *
 *     printf("%.*g", MSP_DBL_DECIMAL_DIG, foo);
 *
 * See https://stackoverflow.com/a/19897395/2752221
 */
#ifdef DBL_DECIMAL_DIG
#define MSP_DBL_DECIMAL_DIG (DBL_DECIMAL_DIG)
#else
#define MSP_DBL_DECIMAL_DIG (DBL_DIG + 3)
#endif

/* Node flags */
#define MSP_NODE_IS_SAMPLE 1u

/* The root node indicator */
#define MSP_NULL_NODE (-1)
/* Indicates the that the population ID has not been set. */
#define MSP_NULL_POPULATION (-1)
/* There is no parent for a given mutation */
#define MSP_NULL_MUTATION (-1)
/* Indicates that no individual has been set */
#define MSP_NULL_INDIVIDUAL (-1)

/* Flags for simplify() */
#define MSP_FILTER_SITES                 (1 << 0)
#define MSP_REDUCE_TO_SITE_TOPOLOGY      (1 << 1)
#define MSP_FILTER_POPULATIONS           (1 << 2)
#define MSP_FILTER_INDIVIDUALS           (1 << 3)

/* Flags for check_integrity */
#define MSP_CHECK_OFFSETS                (1 << 0)
#define MSP_CHECK_EDGE_ORDERING          (1 << 1)
#define MSP_CHECK_SITE_ORDERING          (1 << 2)
#define MSP_CHECK_SITE_DUPLICATES        (1 << 3)
#define MSP_CHECK_MUTATION_ORDERING      (1 << 4)
#define MSP_CHECK_INDEXES                (1 << 5)
#define MSP_CHECK_ALL                    \
    (MSP_CHECK_OFFSETS | MSP_CHECK_EDGE_ORDERING | MSP_CHECK_SITE_ORDERING | \
     MSP_CHECK_SITE_DUPLICATES | MSP_CHECK_MUTATION_ORDERING | MSP_CHECK_INDEXES)

/* Flags for dump tables */
#define MSP_ALLOC_TABLES 1

/* Flags for load tables */
#define MSP_BUILD_INDEXES 1

#define MSP_LOAD_EXTENDED_CHECKS  1

#define MSP_FILE_FORMAT_NAME          "tskit.trees"
#define MSP_FILE_FORMAT_NAME_LENGTH   11
#define MSP_FILE_FORMAT_VERSION_MAJOR 12
#define MSP_FILE_FORMAT_VERSION_MINOR 0

/* Error codes */
#define MSP_ERR_GENERIC                                             -1
#define MSP_ERR_NO_MEMORY                                           -2
#define MSP_ERR_BAD_STATE                                           -3
#define MSP_ERR_BAD_PARAM_VALUE                                     -4
#define MSP_ERR_OUT_OF_BOUNDS                                       -5
#define MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS                         -6
#define MSP_ERR_POPULATION_OVERFLOW                                 -7
#define MSP_ERR_LINKS_OVERFLOW                                      -8
#define MSP_ERR_POPULATION_OUT_OF_BOUNDS                            -9
#define MSP_ERR_BAD_POPULATION_CONFIGURATION                        -10
#define MSP_ERR_BAD_MIGRATION_MATRIX                                -11
#define MSP_ERR_BAD_MIGRATION_MATRIX_INDEX                          -12
#define MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX                     -13
#define MSP_ERR_INFINITE_WAITING_TIME                               -14
#define MSP_ERR_ASSERTION_FAILED                                    -15
#define MSP_ERR_SOURCE_DEST_EQUAL                                   -16
#define MSP_ERR_BAD_RECOMBINATION_MAP                               -17
#define MSP_ERR_BAD_POPULATION_SIZE                                 -18
#define MSP_ERR_BAD_SAMPLES                                         -19
#define MSP_ERR_BAD_MODEL                                           -20
#define MSP_ERR_INSUFFICIENT_SAMPLES                                -21
#define MSP_ERR_DUPLICATE_SITE_POSITION                             -22
#define MSP_ERR_UNDEFINED_MULTIPLE_MERGER_COALESCENT                -23
#define MSP_ERR_INCOMPATIBLE_FROM_TS                                -24
#define MSP_ERR_BAD_START_TIME_FROM_TS                              -25
#define MSP_ERR_BAD_START_TIME                                      -26
#define MSP_ERR_BAD_DEMOGRAPHIC_EVENT_TIME                          -27
#define MSP_ERR_RECOMB_MAP_TOO_COARSE                               -28
#define MSP_ERR_TIME_TRAVEL                                         -29
#define MSP_ERR_INTEGRATION_FAILED                                  -30
#define MSP_ERR_BAD_SWEEP_LOCUS                                     -31
#define MSP_ERR_BAD_TRAJECTORY_TIME                                 -32
#define MSP_ERR_BAD_TRAJECTORY_ALLELE_FREQUENCY                     -33
#define MSP_ERR_EMPTY_TRAJECTORY                                    -34

/* This bit is 0 for any errors originating from tskit */
#define MSP_TSK_ERR_BIT 13

int msp_set_tsk_error(int err);
bool msp_is_tsk_error(int err);
const char * msp_strerror(int err);
void __msp_safe_free(void **ptr);

#define msp_safe_free(pointer) __msp_safe_free((void **) &(pointer))

#endif /*__UTIL_H__*/
