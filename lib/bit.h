/*
** Copyright (C) 2014 Jerome Kelleher <jerome.kelleher@ed.ac.uk>
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

#ifndef __BIT_H__
#define __BIT_H__

typedef struct bit_t_t {
    unsigned int max_index;
    unsigned int log_max_index;
    long long *tree;
    long long *values;
} bit_t;


void bit_alloc(bit_t *);
void bit_free(bit_t *);
long long bit_get_total(bit_t *);
void bit_increment(bit_t *, unsigned int, long long);
void bit_set_value(bit_t *, unsigned int, long long);
long long bit_get_cumulative_sum(bit_t *, unsigned int);
long long bit_get_value(bit_t *, unsigned int);
unsigned int bit_find(bit_t *, long long);

#endif /*__BIT_H__*/
