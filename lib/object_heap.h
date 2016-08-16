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

#include "msprime.h"

size_t object_heap_get_num_allocated(object_heap_t *self);
void object_heap_print_state(object_heap_t *self, FILE *out);
int object_heap_expand(object_heap_t *self);
void * object_heap_get_object(object_heap_t *self, size_t index);
int object_heap_empty(object_heap_t *self);
void * object_heap_alloc_object(object_heap_t *self);
void object_heap_free_object(object_heap_t *self, void *obj);
int object_heap_init(object_heap_t *self, size_t object_size, size_t block_size,
        void (*init_object)(void **, size_t));
void object_heap_free(object_heap_t *self);
