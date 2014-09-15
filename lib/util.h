/*
** Copyright (C) 2009 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
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

#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdlib.h>

#ifdef DEBUG_ERRORS
#define ERROR_CHECK(r, label) \
    printf("ERROR_CHECK: %d: %d: %s\n", r, __LINE__, __FILE__); \
    if (r < 0) { goto label;}
#else
#define ERROR_CHECK(r, label) \
    if (r < 0) { goto label;}
#endif



void * xmalloc(size_t size);
void * xcalloc(size_t count, size_t eltsize);
void * xrealloc(void *ptr, size_t size);
        
int parse_double(const char *str, double *value, const double min, 
        const double max);
int parse_long(const char *str, long *value, const long min, 
        const long max);

void fatal_error(const char *msg, ...);

#endif /*__UTIL_H__*/
