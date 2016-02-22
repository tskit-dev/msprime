/*
** Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/utsname.h>

#include <hdf5.h>
#include <gsl/gsl_version.h>

#include "err.h"
#include "msprime.h"

int
msp_encode_environment(char **result)
{
    int ret = -1;
    /* TODO add more environment: endianess, and word size at a minimum */
    const char *pattern = "{"
        "\"msprime_version\":\"%s\", "
        "\"hdf5_version\":\"%d.%d.%d\", "
        "\"gsl_version\":\"%d.%d\", "
        "\"kernel_name\":\"%s\", "
        "\"kernel_release\":\"%s\", "
        "\"kernel_version\":\"%s\", "
        "\"hardware_identifier\":\"%s\""
        "}";
    herr_t status;
    unsigned int major, minor, release;
    int written;
    size_t size;
    char *str;
    struct utsname system_info;

    if (uname(&system_info) < 0) {
        ret = MSP_ERR_IO;
        goto out;
    }
    status = H5get_libversion(&major, &minor, &release);
    if (status != 0) {
        goto out;
    }
    size = 1 + (size_t) snprintf(NULL, 0, pattern,
            MSP_LIBRARY_VERSION_STR,
            major, minor, release,
            GSL_MAJOR_VERSION, GSL_MINOR_VERSION,
            system_info.sysname, system_info.release, system_info.version,
            system_info.machine);
    str = malloc(size);
    if (str == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    written = snprintf(str, size, pattern,
            MSP_LIBRARY_VERSION_STR,
            major, minor, release,
            GSL_MAJOR_VERSION, GSL_MINOR_VERSION,
            system_info.sysname, system_info.release, system_info.version,
            system_info.machine);
    if (written < 0) {
        ret = MSP_ERR_IO;
        goto out;
    }
    assert(written == (int) size - 1);
    *result = str;
    ret = 0;
out:
    return ret;
}
