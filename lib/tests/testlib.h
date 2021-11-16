/*
** Copyright (C) 2016-2020 University of Oxford
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

#ifndef __TESTLIB_H__
#define __TESTLIB_H__

#include "msprime.h"
#include "likelihood.h"

#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <CUnit/Basic.h>

#ifdef __SSE2__
/* On CPUs lacking SSE2 instructions we can't do exact floating
 * point comparisons. See https://gcc.gnu.org/wiki/FloatingPointMath
 */
#define MSP_TEST_EXACT_FLOAT_COMPARISONS
#endif

#define ALPHABET_BINARY 0
#define ALPHABET_NUCLEOTIDE 1

/* Global variables used for test in state in the test suite */
extern char *_tmp_file_name;
extern FILE *_devnull;

typedef struct {
    population_id_t population;
    double time;
} sample_t;

int build_sim(msp_t *msp, tsk_table_collection_t *tables, gsl_rng *rng,
    double sequence_length, size_t num_populations, sample_t *samples,
    size_t num_samples);
int build_pedigree_sim(msp_t *msp, tsk_table_collection_t *tables, gsl_rng *rng,
    double sequence_length, size_t ploidy, size_t num_individuals, tsk_id_t *parents,
    double *time, tsk_flags_t *is_sample, tsk_id_t *population);
gsl_rng *safe_rng_alloc(void);

int test_main(CU_TestInfo *tests, int argc, char **argv);

#endif
