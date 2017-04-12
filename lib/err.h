/*
** Copyright (C) 2015-2016 University of Oxford
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
#ifndef __ERR_H__
#define __ERR_H__

/*
 * raise a compiler warning if a potentially error raising function's return
 * value is not used.
 */
#ifdef __GNUC__
    #define WARN_UNUSED __attribute__ ((warn_unused_result))
#else
    #define WARN_UNUSED
#endif

#define MSP_ERR_GENERIC                                             -1
#define MSP_ERR_NO_MEMORY                                           -2
#define MSP_ERR_IO                                                  -3
#define MSP_ERR_FILE_FORMAT                                         -4
#define MSP_ERR_REFCOUNT_NONZERO                                    -5
#define MSP_ERR_BAD_STATE                                           -6
#define MSP_ERR_BAD_PARAM_VALUE                                     -7
#define MSP_ERR_OUT_OF_BOUNDS                                       -8
#define MSP_ERR_NEWICK_OVERFLOW                                     -9
#define MSP_ERR_UNSORTED_DEMOGRAPHIC_EVENTS                         -10
#define MSP_ERR_POPULATION_OVERFLOW                                 -11
#define MSP_ERR_LINKS_OVERFLOW                                      -12
#define MSP_ERR_HDF5                                                -13
#define MSP_ERR_BAD_POPULATION_ID                                   -14
#define MSP_ERR_DUPLICATE_SAMPLE                                    -15
#define MSP_ERR_BAD_ORDERING                                        -16
#define MSP_ERR_BAD_MUTATION                                        -17
#define MSP_ERR_UNSUPPORTED_OPERATION                               -18
#define MSP_ERR_BAD_POPULATION_CONFIGURATION                        -19
#define MSP_ERR_BAD_MIGRATION_MATRIX                                -20
#define MSP_ERR_BAD_MIGRATION_MATRIX_INDEX                          -21
#define MSP_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX                     -22
#define MSP_ERR_INFINITE_WAITING_TIME                               -23
#define MSP_ERR_ASSERTION_FAILED                                    -24
#define MSP_ERR_SOURCE_DEST_EQUAL                                   -25
#define MSP_ERR_BAD_RECOMBINATION_MAP                               -26
#define MSP_ERR_BAD_COALESCENCE_RECORDS                             -27
#define MSP_ERR_BAD_SAMPLES                                         -28
#define MSP_ERR_NONBINARY_NEWICK                                    -29
#define MSP_ERR_FILE_VERSION_TOO_OLD                                -30
#define MSP_ERR_FILE_VERSION_TOO_NEW                                -31
#define MSP_ERR_CANNOT_SIMPLIFY                                     -32
#define MSP_ERR_BAD_MODEL                                           -33
#define MSP_ERR_ZERO_CHILDREN                                       -34
#define MSP_ERR_UNSORTED_CHILDREN                                   -35
#define MSP_ERR_NULL_NODE_IN_RECORD                                 -36
#define MSP_ERR_BAD_NODE_TIME_ORDERING                              -37
#define MSP_ERR_RECORDS_NOT_TIME_SORTED                             -38
#define MSP_ERR_INCONSISTENT_NODE_TIMES                             -39
#define MSP_ERR_INCONSISTENT_POPULATION_IDS                         -40
#define MSP_ERR_BAD_RECORD_INTERVAL                                 -41
#define MSP_ERR_ZERO_RECORDS                                        -42
#define MSP_ERR_NOT_INITIALISED                                     -43
/* #define MSP_ERR_MUTATIONS_NOT_POSITION_SORTED                       -44 */
#define MSP_ERR_UNSORTED_MUTATION_NODES                             -45
#define MSP_ERR_DUPLICATE_MUTATION_NODES                            -46
#define MSP_ERR_NONBINARY_MUTATIONS_UNSUPPORTED                     -47
#define MSP_ERR_INCONSISTENT_MUTATIONS                              -48
#define MSP_ERR_BAD_COALESCENCE_RECORDS_SAMPLE_SIZE                 -49
#define MSP_ERR_BAD_COALESCENCE_RECORDS_SEQUENCE_LENGTH             -50
#define MSP_ERR_BAD_COALESCENCE_RECORD_NONMATCHING_RIGHT            -51
#define MSP_ERR_COORDINATE_NOT_FOUND                                -52
#define MSP_ERR_BAD_NODES_ARRAY                                     -53
#define MSP_ERR_BAD_CHILDREN_ARRAY                                  -54
#define MSP_ERR_SITE_OUT_OF_BOUNDS                                  -55
#define MSP_ERR_NODE_OUT_OF_BOUNDS                                  -56
#define MSP_ERR_LENGTH_MISMATCH                                     -57
#define MSP_ERR_BAD_ALPHABET                                        -58
#define MSP_ERR_NON_SINGLE_CHAR_MUTATION                            -59
#define MSP_ERR_UNSORTED_SITES                                      -60
#define MSP_ERR_BAD_SITE_POSITION                                   -61
#define MSP_ERR_UNSORTED_MUTATIONS                                  -62
#define MSP_ERR_UNDEFINED_MULTIPLE_MERGER_COALESCENT                -63

/* TODO need to go through these and rename all the COALESCENCE_RECORD ones to
 * have better names. There are also unused error numbers here which should be
 * reused.
 */

#endif /*__ERR_H__*/
