/*
** Copyright (C) 2015-2019 University of Oxford
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

#include "tsk_core.h"
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
        case MSP_ERR_LINKS_OVERFLOW:
            ret = "Links Overflow occurred.";
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
        case MSP_ERR_BAD_RECOMBINATION_MAP:
            ret = "Bad recombination map provided.";
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
        case MSP_ERR_RECOMB_MAP_TOO_COARSE:
            ret = "The specified recombination map is cannot translate the coordinates"
                "for the specified tree sequence. It is either too coarse (num_loci "
                "is too small) or contains zero recombination rates. Please either "
                "increase the number of loci or recombination rate";
            break;
        case MSP_ERR_TIME_TRAVEL:
            ret = "The simulation model supplied resulted in a parent node having "
                "a time value <= to its child. This can occur either as a result "
                "of multiple bottlenecks happening at the same time or because of "
                "numerical imprecision with very small population sizes.";
            break;
        case MSP_ERR_INTEGRATION_FAILED:
            ret = "GSL numerical integration failed. Please check the stderr for details.";
            break;
        case MSP_ERR_BAD_SWEEP_LOCUS:
            ret = "Sweep locus must be between 0 and num_loci.";
            break;
        case MSP_ERR_BAD_TRAJECTORY_TIME:
            ret = "Time values must be > 0 and in increasing order.";
            break;
        case MSP_ERR_BAD_TRAJECTORY_ALLELE_FREQUENCY:
            ret = "Allele frequency values must be between 0 and 1.";
            break;
        case MSP_ERR_EMPTY_TRAJECTORY:
            ret = "Trajectory must contain at least one time point.";
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
    return !(err & (1 << MSP_TSK_ERR_BIT));
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
__msp_safe_free(void **ptr) {
    if (ptr != NULL) {
        if (*ptr != NULL) {
            free(*ptr);
            *ptr = NULL;
        }
    }
}
