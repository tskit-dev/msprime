/*
** Copyright (C) 2018 University of Oxford
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
#include <stdlib.h>
#include <stdint.h>

#include "util.h"
#include "uuid.h"


#define NUM_BYTES 16

static int WARN_UNUSED
get_random_bytes(uint8_t *buf)
{
    int ret = MSP_ERR_GENERATE_UUID;
    FILE *f = fopen("/dev/urandom", "r");

    if (f == NULL) {
        goto out;
    }
    if (fread(buf, NUM_BYTES, 1, f) != 1) {
        goto out;
    }
    if (fclose(f) != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

/* Generate a new UUID4 using a system-generated source of randomness.
 * Note that this function writes a NULL terminator to the end of this
 * string, so that the total length of the buffer must be 37 bytes.
 */
int
tsk_generate_uuid(char *dest, int MSP_UNUSED(flags))
{
    int ret = 0;
    uint8_t buf[NUM_BYTES];
    const char *pattern =
        "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x";

    ret = get_random_bytes(buf);
    if (ret != 0) {
        goto out;
    }
    if (snprintf(dest, TSK_UUID_SIZE + 1, pattern,
            buf[0], buf[1], buf[2], buf[3],
            buf[4], buf[5], buf[6], buf[7],
            buf[8], buf[9], buf[10], buf[11],
            buf[12], buf[13], buf[14], buf[15]) < 0) {
        ret = MSP_ERR_GENERATE_UUID;
        goto out;
    }
out:
    return ret;
}
