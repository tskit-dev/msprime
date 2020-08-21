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

#define ALPHABET_BINARY 0
#define ALPHABET_NUCLEOTIDE 1

/* Global variables used for test in state in the test suite */
char *_tmp_file_name;
FILE *_devnull;

int test_main(CU_TestInfo *tests, int argc, char **argv);

#endif
