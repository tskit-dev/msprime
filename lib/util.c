/*
** Copyright (C) 2008 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software 
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/
/*
 * Collection of miscellaneous functions shared throughout source.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "util.h"

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <stdarg.h>


void *
xmalloc(size_t size)
{
    register void *value = malloc(size);
    if (value == NULL) {
        perror("virtual memory exhausted");
        abort();
    }
    return value;
}

void *
xcalloc(size_t count, size_t eltsize)
{
    register void *value = calloc(count, eltsize);
    if (value == NULL) {
        perror("virtual memory exhausted");
        abort();
    }
    return value;
}

void *
xrealloc(void *ptr, size_t size)
{
    register void *value = realloc(ptr, size);
    if (value == NULL) {
        perror("Virtual memory exhausted");
        abort();
    }
    return value;
}

/*
 * Parses the specified string into a double and assigns the value into 
 * the specified pointer. Returns EINVAL if the string cannot be 
 * converted to double or if min <= x <= max does not hold; returns 0 if 
 * the value is converted successfully.
 */
int 
parse_double(const char *str, double *value, const double min, 
        const double max)
{
    int ret = 0;
    double x;
    char *tail; 
    x = strtod(str, &tail);
    if (tail[0] != '\0') {
        ret = EINVAL; 
    } else if (min > x || max < x) {
        ret = EINVAL;
    } else {
        *value = x;
    }
    return ret;
}

/*
 * Parses the specified string into a long and assigns the value into 
 * the specified pointer. Returns EINVAL if the string cannot be 
 * converted to double or if min <= x <= max does not hold; returns 0 if 
 * the value is converted successfully.
 */
int 
parse_long(const char *str, long *value, const long min, 
        const long max)
{
    int ret = 0;
    long x;
    char *tail; 
    x = strtol(str, &tail, 10);
    if (tail[0] != '\0') {
        ret = EINVAL; 
    } else if (min > x || max < x) {
        ret = EINVAL;
    } else {
        *value = x;
    }
    return ret;
}

void 
fatal_error(const char *msg, ...)
{
    va_list argp;
    fprintf(stderr, "sms:");
    va_start(argp, msg);
    vfprintf(stderr, msg, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}

