/*
** Copyright (C) 2014 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_int.h>

#include "err.h"
#include "avl.h"
#include "fenwick.h"
#include "msprime.h"

static int
cmp_individual(const void *a, const void *b) {
    const segment_t *ia = (const segment_t *) a;
    const segment_t *ib = (const segment_t *) b;
    return (ia->index > ib->index) - (ia->index < ib->index);
}

static int
cmp_node_mapping(const void *a, const void *b) {
    const node_mapping_t *ia = (const node_mapping_t *) a;
    const node_mapping_t *ib = (const node_mapping_t *) b;
    return (ia->left > ib->left) - (ia->left < ib->left);
}

static inline unsigned int
num_links(segment_t *u)
{
    return u->right - u->left;
}

static void
segment_init(void **obj, size_t index)
{
    segment_t *seg = (segment_t *) obj;
    seg->index = index + 1;
    printf("set seg %p to %d\n", seg, (int) index + 1);
}

/* memory heap manager */

void
object_heap_print_state(object_heap_t *self)
{
    memory_block_t *mb;
    int num_blocks = 0;

    mb = self->mem_blocks;
    while (mb != NULL) {
        mb = mb->next;
        num_blocks++;
    }
    assert(self->num_blocks == num_blocks);
    printf("object heap %p::\n", self);
    printf("\tsize = %d\n", (int) self->size);
    printf("\ttop = %d\n", (int) self->top);
    printf("\tblock_size = %d\n", (int) self->block_size);
    printf("\tnum_blocks = %d\n", num_blocks);
}

static void
object_heap_add_block(object_heap_t *self, memory_block_t *block)
{
    size_t j, index;

    for (j = 0; j < self->block_size; j++) {
        self->heap[j] = block->mem + j * self->object_size;
        if (self->init_object != NULL) {
            index = j + (self->num_blocks - 1) * self->block_size;
            self->init_object(self->heap[j], index);
        }
    }
    self->top = self->block_size;
}

static int WARN_UNUSED
object_heap_expand(object_heap_t *self)
{
    int ret = -1;
    memory_block_t *mb = malloc(sizeof(memory_block_t));
    memory_block_t *prev;

    if (mb == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    prev = self->mem_blocks;
    while (prev->next != NULL) {
        prev = prev->next;
    }
    prev->next = mb;
    mb->next = NULL;
    mb->mem = malloc(self->block_size * self->object_size);
    if (mb->mem == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->num_blocks++;
    /* Now we increase the size of the heap. Since it is currently empty,
     * we avoid the copying cost of realloc and free before making a new
     * heap.
     */
    free(self->heap);
    self->heap = NULL;
    self->size += self->block_size;
    self->heap = malloc(self->size * sizeof(void **));
    if (self->heap == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    object_heap_add_block(self, mb);
    ret = 0;
out:
    return ret;
}

/*
 * Returns the jth object in the memory buffers.
 */
static inline void * WARN_UNUSED
object_heap_get_object(object_heap_t *self, size_t index)
{
    void *ret = NULL;
    size_t block, obj, j;
    memory_block_t *mb = self->mem_blocks;

    block = index / self->block_size;
    obj = index % self->block_size;
    j = 0;
    while (j < block && mb != NULL) {
        j++;
        mb = mb->next;
    }
    if (j == block) {
        ret = mb->mem + obj * self->object_size;
    }

    return ret;
}

static inline void * WARN_UNUSED
object_heap_alloc_object(object_heap_t *self)
{
    void *ret = NULL;

    if (self->top == 0) {
        if (object_heap_expand(self) != 0) {
            goto out;
        }
    }
    self->top--;
    ret = self->heap[self->top];
out:
    return ret;
}

static inline void
object_heap_free_object(object_heap_t *self, void *obj)
{
    assert(self->top < self->size);
    self->heap[self->top] = obj;
    self->top++;
}

static int WARN_UNUSED
object_heap_init(object_heap_t *self, size_t object_size, size_t block_size,
        void (*init_object)(void **, size_t))
{
    int ret = -1;

    memset(self, 0, sizeof(object_heap_t));
    self->block_size = block_size;
    self->size = block_size;
    self->object_size = object_size;
    self->init_object = init_object;
    self->num_blocks = 1;
    self->heap = malloc(self->size * sizeof(void **));
    self->mem_blocks = malloc(sizeof(memory_block_t));
    if (self->heap == NULL || self->mem_blocks == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->mem_blocks->mem = malloc(self->size * self->object_size);
    if (self->mem_blocks->mem == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->top = 0;
    self->mem_blocks->next = NULL;
    object_heap_add_block(self, self->mem_blocks);
    ret = 0;
out:
    return ret;
}

static void
object_heap_free(object_heap_t *self)
{
    memory_block_t *mb = self->mem_blocks;
    memory_block_t *tmp;

    while (mb != NULL) {
        tmp = mb;
        if (mb->mem != NULL) {
            free(mb->mem);
        }
        mb = mb->next;
        free(tmp);
    }
    if (self->heap != NULL) {
        free(self->heap);
    }
}

/* population models */

static double
constant_population_model_get_size(population_model_t *self, double t)
{
    return self->param;
}

static double
constant_population_model_get_waiting_time(population_model_t *self,
        double lambda_coancestry, double t, gsl_rng* rng)
{
    return self->param * gsl_ran_exponential(rng, 1.0 / lambda_coancestry);
}

static double
exponential_population_model_get_size(population_model_t *self, double t)
{
    double alpha = self->param;
    return exp(-alpha * (t - self->start_time));
}

static double
exponential_population_model_get_waiting_time(population_model_t *self,
        double lambda_coancestry, double t, gsl_rng* rng)
{
    double ret = DBL_MAX;
    double alpha = self->param;
    double dt = t - self->start_time;
    double z = alpha * self->initial_size * exp(-alpha * dt);
    z = 1 + z * gsl_ran_exponential(rng, 1.0 / lambda_coancestry);
    /* if z is <= 0 no coancestry can occur */
    if (z > 0) {
        ret = log(z) / alpha;
    }
    return ret;
}

static int WARN_UNUSED
msp_add_population_model(msp_t *self, population_model_t *model)
{
    int ret = -1;
    population_model_t *m = self->population_models;
    if (m == NULL) {
        self->population_models = model;
    } else {
        // TODO check for sortedness
        while (m->next != NULL) {
            m = m->next;
        }
        m->next = model;
    }
    model->next = NULL;
    ret = 0;

    return ret;
}

int
msp_add_constant_population_model(msp_t *self, double time, double size)
{
    int ret = -1;
    population_model_t *model = malloc(sizeof(population_model_t));

    // TODO check for model specific restrictions.
    if (model == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    model->start_time = time;
    model->param = size;
    model->type = POP_MODEL_CONSTANT;
    model->get_size = constant_population_model_get_size;
    model->get_waiting_time = constant_population_model_get_waiting_time;
    ret = msp_add_population_model(self, model);
    if (ret != 0) {
        free(model);
    }
out:
    return ret;
}

int
msp_add_exponential_population_model(msp_t *self, double time, double alpha)
{
    int ret = -1;
    population_model_t *model = malloc(sizeof(population_model_t));

    // TODO check for model specific restrictions.
    if (model == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    model->start_time = time;
    model->param = alpha;
    model->type = POP_MODEL_EXPONENTIAL;
    model->get_size = exponential_population_model_get_size;
    model->get_waiting_time = exponential_population_model_get_waiting_time;
    ret = msp_add_population_model(self, model);
    if (ret != 0) {
        free(model);
    }
out:
    return ret;
}

static memory_block_t *
msp_alloc_node_mapping_block(msp_t *self)
{
    memory_block_t *ret = NULL;
    memory_block_t *mb;

    mb = malloc(sizeof(memory_block_t));
    if (mb == NULL) {
        goto out;
    }
    mb->next = NULL;
    mb->mem = malloc(self->node_mapping_block_size * sizeof(node_mapping_t));
    if (mb->mem == NULL) {
        free(mb);
        goto out;
    }
    ret = mb;
out:
    return ret;
}

int
msp_alloc(msp_t *self)
{
    int ret = -1;
    int err, n;
    //segment_t *seg;

    /* turn off GSL error handler so we don't abort on memory error */
    gsl_set_error_handler_off();
    self->rng = gsl_rng_alloc(gsl_rng_default);
    if (self->rng == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    gsl_rng_set(self->rng, self->random_seed);
    /* Allocate the memory heaps */
    ret = object_heap_init(&self->avl_node_heap, sizeof(avl_node_t),
           self->avl_node_block_size, NULL);
    if (ret != 0) {
        goto out;
    }
    /* set up the node mapping allocator */
    self->current_node_mapping_block = msp_alloc_node_mapping_block(self);
    if (self->current_node_mapping_block == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->node_mapping_blocks = self->current_node_mapping_block;
    self->next_node_mapping = 0;

    ret = object_heap_init(&self->segment_heap, sizeof(segment_t),
           self->max_segments, segment_init);
    if (ret != 0) {
        goto out;
    }

    /* old memory management */
    n = self->max_segments;
    /*
    self->segment_mem = malloc(n * sizeof(segment_t));
    self->segment_heap = malloc(n * sizeof(segment_t *));
    if (self->segment_mem == NULL || self->segment_heap == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    for (j = 0; j < n; j++) {
        seg = &self->segment_mem[j];
        self->segment_heap[j] = seg;
        seg->index = j + 1;
    }
    self->segment_heap_top = (int) n - 1;
    */
    self->links = calloc(1, sizeof(fenwick_t));
    if (self->links == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    self->links->max_index = n;
    err = fenwick_alloc(self->links);
    if (err != 0) {
        ret = err;
        goto out;
    }
    /* set up the AVL trees */
    self->ancestral_population = malloc(sizeof(avl_tree_t));
    self->breakpoints = malloc(sizeof(avl_tree_t));
    if (self->ancestral_population == NULL || self->breakpoints == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_tree(self->ancestral_population, cmp_individual, NULL);
    avl_init_tree(self->breakpoints, cmp_node_mapping, NULL);
    ret = 0;
out:
    return ret;
}

int
msp_free(msp_t *self)
{
    int ret = -1;
    population_model_t *u, *v;
    memory_block_t *mb, *prev;

    u = self->population_models;
    while (u != NULL) {
        v = u->next;
        free(u);
        u = v;
    }
    if (self->rng != NULL) {
        gsl_rng_free(self->rng);
    }
    if (self->links != NULL) {
        fenwick_free(self->links);
        free(self->links);
    }
    if (self->ancestral_population != NULL) {
        free(self->ancestral_population);
    }
    if (self->breakpoints != NULL) {
        free(self->breakpoints);
    }
    /* free the object heaps */
    object_heap_free(&self->avl_node_heap);
    object_heap_free(&self->segment_heap);
    mb = self->node_mapping_blocks;
    while (mb != NULL) {
        prev = mb;
        mb = mb->next;
        free(prev->mem);
        free(prev);
    }

    /* old */

    /*
    if (self->segment_mem != NULL) {
        free(self->segment_mem);
    }
    if (self->segment_heap != NULL) {
        free(self->segment_heap);
    }
    */
    if (self->coalescence_record_file != NULL) {
        // TODO check for errors here?
        fclose(self->coalescence_record_file);
    }
    ret = 0;
    return ret;
}


static inline avl_node_t * WARN_UNUSED
msp_alloc_avl_node(msp_t *self)
{
    return (avl_node_t *) object_heap_alloc_object(&self->avl_node_heap);
}

static inline void
msp_free_avl_node(msp_t *self, avl_node_t *node)
{
    object_heap_free_object(&self->avl_node_heap, node);
}

static inline node_mapping_t *
msp_alloc_node_mapping(msp_t *self)
{
    node_mapping_t *ret = NULL;
    memory_block_t *mb = self->current_node_mapping_block;

    if (self->next_node_mapping == self->node_mapping_block_size) {
        mb = msp_alloc_node_mapping_block(self);
        if (mb == NULL) {
            goto out;
        }
        self->current_node_mapping_block->next = mb;
        self->current_node_mapping_block = mb;
        self->next_node_mapping = 0;
    }
    ret = mb->mem + self->next_node_mapping * sizeof(node_mapping_t);
    self->next_node_mapping++;
out:
    return ret;
}


static segment_t * WARN_UNUSED
msp_alloc_segment(msp_t *self, int left, int right, int value, segment_t *prev,
        segment_t *next)
{
    /*
    segment_t *seg = self->segment_heap[self->segment_heap_top];
    self->segment_heap_top--;
    if (self->segment_heap_top < 0) {
        seg = NULL;
        goto out;
    }
    */
    segment_t *seg = (segment_t *) object_heap_alloc_object(
            &self->segment_heap);
    if (seg == NULL) {
        goto out;
    }
    seg->prev = prev;
    seg->next = next;
    seg->left = left;
    seg->right = right;
    seg->value = value;
out:
    return seg;
}

/*
 * Returns the segment with the specified index.
 */
static segment_t *
msp_get_segment(msp_t *self, unsigned int index)
{
    //segment_t *u = &self->segment_mem[index - 1];
    segment_t *u = object_heap_get_object(&self->segment_heap, index - 1);
    assert(u != NULL);
    printf("looking for %d, got %d %p\n", index, u->index, u);
    assert(u->index == index);
    return u;
}

static void
msp_free_segment(msp_t *self, segment_t *node)
{
    /*
    self->segment_heap_top++;
    assert(self->segment_heap_top < (int) self->max_segments);
    self->segment_heap[self->segment_heap_top] = node;
    */
    object_heap_free_object(&self->segment_heap, node);
    fenwick_set_value(self->links, node->index, 0);
}

static int WARN_UNUSED
msp_open_coalescence_record_file(msp_t *self)
{
    int ret = 0;
    FILE *f = fopen(self->coalescence_record_filename, "w");
    if (f == NULL) {
        ret = MSP_ERR_IO;
        goto out;
    }
    self->coalescence_record_file = f;
out:
    return ret;
}


static inline int WARN_UNUSED
msp_insert_individual(msp_t *self, segment_t *u)
{
    int ret = 0;
    avl_node_t *node;

    assert(u != NULL);
    node = msp_alloc_avl_node(self);
    if (node == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    avl_init_node(node, u);
    node = avl_insert_node(self->ancestral_population, node);
    assert(node != NULL);
out:
    return ret;
}

static void
msp_print_segment_chain(msp_t *self, segment_t *head)
{
    int j;
    segment_t *s = head;
    j = 1;
    while (s != NULL) {
        while (j < s->left) {
            printf("%4c", '-');
            j++;
        }
        while (j <= s->right) {
            printf("%4d", s->value);
            j++;
        }
        s = s->next;
    }
    while (j <= self->num_loci) {
        printf("%4c", '-');
        j++;
    }
    printf("\n");
}


static void
msp_verify(msp_t *self)
{
    long long s, ss, total_links;
    unsigned int total_segments = 0;
    //unsigned int total_avl_nodes = 0;
    avl_node_t *node;
    segment_t *u;

    total_links = 0;
    node = self->ancestral_population->head;
    while (node != NULL) {
        u = (segment_t *) node->item;
        assert(u->prev == NULL);
        while (u != NULL) {
            total_segments++;
            assert(u->left <= u->right);
            /*
            for (j = u->left; j <= u->right; j++) {
                msp_get_tree(self, j, pi, tau);
                assert(pi[u->value] == 0);
            }
            */
            if (u->prev != NULL) {
                s = u->right - u->prev->right;
            } else {
                s = u->right - u->left;
            }
            ss = fenwick_get_value(self->links, u->index);
            total_links += ss;
            assert(s == ss);
            u = u->next;
        }
        node = node->next;
    }
    assert(total_links == fenwick_get_total(self->links));
    /*
    assert(total_segments == self->max_segments - self->segment_heap_top - 1);
    total_avl_nodes = avl_count(self->ancestral_population)
            + avl_count(self->breakpoints);
    assert(total_avl_nodes == self->max_avl_nodes -
            self->avl_node_heap_top - 1);
    */
}

int
msp_print_state(msp_t *self)
{
    int ret = 0;
    avl_node_t *node;
    node_mapping_t *nm;
    segment_t *u;
    long long v;
    int j, k, bp, num_blocks;
    int *pi;
    float *tau;
    coalescence_record_t cr;
    memory_block_t *mb;
    FILE *f = fopen(self->coalescence_record_filename, "r");
    tree_viewer_t *tv = msp_get_tree_viewer(self);

    if (f == NULL) {
        ret = MSP_ERR_IO;
        goto out;
    }
    if (tv == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    ret = tree_viewer_init(tv);
    if (ret != 0) {
        goto out;
    }
    printf("n = %d\n", self->sample_size);
    printf("m = %d\n", self->num_loci);
    printf("random seed = %ld\n", self->random_seed);
    printf("num_links = %lld\n", fenwick_get_total(self->links));
    printf("population = %d\n", avl_count(self->ancestral_population));
    printf("time = %f\n", self->time);
    node = self->ancestral_population->head;
    while (node != NULL) {
        u = (segment_t *) node->item;
        printf("\t");
        msp_print_segment_chain(self, u);
        node = node->next;
    }
    printf("trees = %d\n", tv->num_trees);
    bp = 0;
    pi = NULL;
    tau = NULL;
    for (j = 0; j < tv->num_trees; j++) {
        ret = tree_viewer_get_tree(tv, j, &bp, &pi, &tau);
        printf("\t%4d:", bp);
        for (k = 0; k < 2 * self->sample_size; k++) {
            printf("%3d ", pi[k]);
        }
        printf("|");
        for (k = 0; k < 2 * self->sample_size; k++) {
            printf("%.3f ", tau[k]);
        }
        printf("\n");

    }
    printf("Fenwick tree\n");
    for (j = 1; j <= self->max_segments; j++) {
        u = msp_get_segment(self, j);
        v = fenwick_get_value(self->links, j);
        if (v != 0) {
            printf("\t%lld\tl=%d r=%d v=%d prev=%p next=%p\n", v, u->left,
                    u->right, u->value, u->prev, u->next);
        }
    }
    printf("Breakpoints = %d\n", avl_count(self->breakpoints));
    for (node = self->breakpoints->head; node != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        printf("\t%d -> %d\n", nm->left, nm->value);
    }
    printf("Coalescence records = %d\n", self->num_coalescence_records);
    /* seek to the start of file and read all records */
    if (fseek(f, 0, SEEK_SET) != 0) {
        ret = MSP_ERR_IO;
        goto out;
    }
    while ((ret = fread(&cr, sizeof(cr), 1, f)) == 1) {
        printf("\t(%d, %d) -- (%d, %d)->%d @ %f\n", cr.left, cr.right,
                cr.children[0], cr.children[1], cr.parent, cr.time);
    }
    if (!feof(f)) {
        ret = MSP_ERR_IO;
        goto out;
    }
    printf("Memory heaps\n");
    printf("avl_node_heap:");
    object_heap_print_state(&self->avl_node_heap);
    printf("segment_heap:");
    object_heap_print_state(&self->segment_heap);

    mb = self->node_mapping_blocks;
    num_blocks = 0;
    while (mb != NULL) {
        mb = mb->next;
        num_blocks++;
    }
    printf("node_mapping stack:\n");
    printf("\tblock size = %d\n", (int) self->node_mapping_block_size);
    printf("\tnext = %d\n", (int) self->next_node_mapping);
    printf("\tnum_blocks = %d\n", num_blocks);
    msp_verify(self);
out:
    if (f != NULL) {
        fclose(f);
    }
    if (tv != NULL) {
        tree_viewer_free(tv);
    }
    return ret;
}

/*
 * Inserts a new breakpoint at the specified locus left, mapping to the
 * specified tree node v.
 */
static int WARN_UNUSED
msp_insert_breakpoint(msp_t *self, unsigned int left, int v)
{
    int ret = 0;
    avl_node_t *node = msp_alloc_avl_node(self);
    node_mapping_t *m = msp_alloc_node_mapping(self);

    if (node == NULL || m == NULL) {
        ret = MSP_ERR_NO_MEMORY;
        goto out;
    }
    m->left = left;
    m->value = v;
    avl_init_node(node, m);
    node = avl_insert_node(self->breakpoints, node);
    assert(node != NULL);
out:
    return ret;
}

/*
 * Inserts a new breakpoint at the specified locus, and copies it's
 * node mapping from the containing breakpoint.
 */
static int WARN_UNUSED
msp_copy_breakpoint(msp_t *self, unsigned int k)
{
    int ret;
    node_mapping_t search, *nm;
    avl_node_t *node;

    search.left = k;
    avl_search_closest(self->breakpoints, &search, &node);
    assert(node != NULL);
    nm = (node_mapping_t *) node->item;
    if (nm->left > k) {
        node = node->prev;
        assert(node != NULL);
        nm = (node_mapping_t *) node->item;
    }
    ret = msp_insert_breakpoint(self, k, nm->value);
    return ret;
}

static int WARN_UNUSED
msp_record_coalescence(msp_t *self, unsigned int left, unsigned int right,
        int child1, int child2, int parent)
{
    int ret = 0;
    coalescence_record_t r;

    r.left = left;
    r.right = right;
    r.children[0] = child1;
    r.children[1] = child2;
    r.parent = parent;
    r.time = self->time;
    self->num_coalescence_records++;
    if (fwrite(&r, sizeof(r), 1, self->coalescence_record_file) != 1) {
        ret = MSP_ERR_IO;
        goto out;
    }
out:
    return ret;

}

static int WARN_UNUSED
msp_recombination_event(msp_t *self)
{
    int ret = 0;
    long long l, t, gap, k;
    unsigned int j;
    node_mapping_t search;
    segment_t *x, *y, *z;
    long long num_links = fenwick_get_total(self->links);

    self->num_re_events++;
    /* We can't use the GSL integer generator here as the range is too large */
    l = 1 + (long long) (gsl_rng_uniform(self->rng) * num_links);
    assert(l > 0 && l <= num_links);
    j = fenwick_find(self->links, l);
    t = fenwick_get_cumulative_sum(self->links, j);
    gap = t - l;
    y = msp_get_segment(self, j);
    x = y->prev;
    k = y->right - gap - 1;
    if (y->left <= k) {
        z = msp_alloc_segment(self, k + 1, y->right, y->value, NULL, y->next);
        if (z == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        fenwick_increment(self->links, z->index, z->right - k - 1);
        if (y->next != NULL) {
            y->next->prev = z;
        }
        y->next = NULL;
        y->right = k;
        fenwick_increment(self->links, y->index, k - z->right);
        search.left = k + 1;
        if (avl_search(self->breakpoints, &search) == NULL) {
            ret = msp_copy_breakpoint(self, k + 1);
            if (ret != 0) {
                goto out;
            }
        }
    } else {
        x->next = NULL;
        y->prev = NULL;
        t = x->right;
        t -= y->left;
        fenwick_increment(self->links, y->index, t);
        z = y;
        y = x;
        self->num_trapped_re_events++;
    }
    ret = msp_insert_individual(self, z);
out:
    return ret;
}

static int WARN_UNUSED
msp_coancestry_event(msp_t *self)
{
    int ret = 0;
    int eta;
    unsigned int j, n, l, r, r_max;
    avl_node_t *node;
    node_mapping_t *nm, search;
    segment_t *x, *y, *z, *alpha, *beta;

    self->num_ca_events++;
    /* Choose x and y */
    n = avl_count(self->ancestral_population);
    j = gsl_rng_uniform_int(self->rng, n);
    node = avl_at(self->ancestral_population, j);
    assert(node != NULL);
    x = (segment_t *) node->item;
    avl_unlink_node(self->ancestral_population, node);
    msp_free_avl_node(self, node);
    j = gsl_rng_uniform_int(self->rng, n - 1);
    node = avl_at(self->ancestral_population, j);
    assert(node != NULL);
    y = (segment_t *) node->item;
    avl_unlink_node(self->ancestral_population, node);
    msp_free_avl_node(self, node);

    /* update num_links and get ready for loop */
    z = NULL;
    while (x != NULL || y != NULL) {
        alpha = NULL;
        if (x == NULL || y == NULL) {
            if (x != NULL) {
                alpha = x;
                x = NULL;
            }
            if (y != NULL) {
                alpha = y;
                y = NULL;
            }
        } else {
            if (y->left < x->left) {
                beta = x;
                x = y;
                y = beta;
            }
            if (x->right < y->left) {
                alpha = x;
                x = x->next;
                alpha->next = NULL;
            } else if (x->left != y->left) {
                alpha = msp_alloc_segment(self, x->left, y->left - 1, x->value,
                        NULL, NULL);
                if (alpha == NULL) {
                    ret = MSP_ERR_NO_MEMORY;
                    goto out;
                }
                x->left = y->left;
            } else {
                l = x->left;
                search.left = l;
                node = avl_search(self->breakpoints, &search);
                assert(node != NULL);
                nm = (node_mapping_t *) node->item;
                eta = nm->value;
                r_max = GSL_MIN(x->right, y->right);
                nm->value++;
                node = node->next;
                nm = (node_mapping_t *) node->item;
                r = nm->left;
                while (nm->value == eta && r < r_max) {
                    nm->value++;
                    node = node->next;
                    nm = (node_mapping_t *) node->item;
                    r = nm->left;
                }
                r--;
                ret = msp_record_coalescence(self, l, r, x->value, y->value,
                        eta);
                if (ret != 0) {
                    goto out;
                }
                if (eta < 2 * self->sample_size - 1) {
                    alpha = msp_alloc_segment(self, l, r, eta, NULL, NULL);
                    if (alpha == NULL) {
                        ret = MSP_ERR_NO_MEMORY;
                        goto out;
                    }
                }
                if (x->right == r) {
                    beta = x;
                    x = x->next;
                    msp_free_segment(self, beta);
                } else {
                    x->left = r + 1;
                }
                if (y->right == r) {
                    beta = y;
                    y = y->next;
                    msp_free_segment(self, beta);
                } else {
                    y->left = r + 1;
                }
            }
        }
        if (alpha != NULL) {
            l = alpha->left;
            if (z == NULL) {
                ret = msp_insert_individual(self, alpha);
                if (ret != 0) {
                    goto out;
                }
            } else {
                z->next = alpha;
                l = z->right;
            }
            alpha->prev = z;
            z = alpha;
            fenwick_set_value(self->links, alpha->index, alpha->right - l);
        }
    }
out:
    return ret;
}


/*
 * Sets up the initial population and trees.
 */
int WARN_UNUSED
msp_initialise(msp_t *self)
{
    int j;
    int ret = 0;
    segment_t *u;

    /* zero the counters */
    self->num_re_events = 0;
    self->num_ca_events = 0;
    self->num_coalescence_records = 0;
    self->num_trapped_re_events = 0;;
    for (j = 1; j <= self->sample_size; j++) {
        u = msp_alloc_segment(self, 1, self->num_loci, j, NULL, NULL);
        if (u == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        ret = msp_insert_individual(self, u);
        if (ret != 0) {
            goto out;
        }
        fenwick_increment(self->links, u->index, self->num_loci - 1);
    }
    ret = msp_insert_breakpoint(self, 1, self->sample_size + 1);
    if (ret != 0) {
        goto out;
    }
    ret = msp_insert_breakpoint(self, self->num_loci + 1, 0);
    if (ret != 0) {
        goto out;
    }
    /* insert the sentinel population model */
    ret = msp_add_constant_population_model(self, DBL_MAX, 0.0);
    if (ret != 0) {
        goto out;
    }
    self->current_population_model = self->population_models;
    self->time = 0.0;
    ret = msp_open_coalescence_record_file(self);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int WARN_UNUSED
msp_run(msp_t *self, double max_time, unsigned long max_events)
{
    int ret = -1;
    double lambda_c, lambda_r, t_c, t_r, t_wait, pop_size;
    long long num_links;
    unsigned int n = avl_count(self->ancestral_population);
    int (*event_method)(msp_t *);
    population_model_t *pop_model = self->current_population_model;
    unsigned long events = 0;

    while (n > 1 && self->time < max_time && events < max_events) {
        events++;
        num_links = fenwick_get_total(self->links);
        lambda_r = num_links * self->recombination_rate;
        t_r = DBL_MAX;
        if (lambda_r != 0.0) {
            t_r = gsl_ran_exponential(self->rng, 1.0 / lambda_r);
        }
        /* Need to perform n * (n - 1) as a double due to overflow */
        lambda_c = n * ((double) n - 1.0);
        t_c = pop_model->get_waiting_time(pop_model, lambda_c, self->time,
                self->rng);
        if (t_r < t_c) {
            t_wait = t_r;
            event_method = msp_recombination_event;
        } else {
            t_wait = t_c;
            event_method = msp_coancestry_event;
        }
        // TODO check for infinite waiting time
        if (self->time + t_wait >= pop_model->next->start_time) {
            pop_size = pop_model->get_size(pop_model, self->time);
            pop_model = pop_model->next;
            pop_model->initial_size = pop_size;
            self->time = pop_model->start_time;
            self->current_population_model = pop_model;
            assert(pop_model->next != NULL);
        } else {
            self->time += t_wait;
            ret = event_method(self);
            if (ret != 0) {
                goto out;
            }
            n = avl_count(self->ancestral_population);
        }
    }
    if (fflush(self->coalescence_record_file) != 0) {
        ret = MSP_ERR_IO;
        goto out;
    }
    // TODO we probably want a different protocol to indicate if max_time
    // has been exceeded or max_events.
    ret = n;
out:
    return ret;
}


tree_viewer_t *
msp_get_tree_viewer(msp_t *self)
{
    int j = 0;
    avl_node_t *node;
    node_mapping_t *nm;
    tree_viewer_t *tv = calloc(1, sizeof(tree_viewer_t));

    if (tv == NULL) {
        goto out;
    }
    tv->sample_size = self->sample_size;
    tv->coalescence_record_filename = malloc(1 +
            strlen(self->coalescence_record_filename));
    if (tv->coalescence_record_filename == NULL) {
        tree_viewer_free(tv);
        goto out;
    }
    strcpy(tv->coalescence_record_filename, self->coalescence_record_filename);
    /* Skip the last tree as it's not real */
    tv->num_trees = avl_count(self->breakpoints) - 1;
    tv->breakpoints = malloc(tv->num_trees * sizeof(int));
    if (tv->breakpoints == NULL) {
        tree_viewer_free(tv);
        tv = NULL;
        goto out;
    }
    for (node = self->breakpoints->head; node->next != NULL; node = node->next) {
        nm = (node_mapping_t *) node->item;
        tv->breakpoints[j] = nm->left;
        j++;
    }
    tv->pi = calloc(2 * tv->sample_size * tv->num_trees, sizeof(int));
    if (tv->pi == NULL) {
        tree_viewer_free(tv);
        tv = NULL;
        goto out;
    }
    tv->tau = calloc(2 * tv->sample_size * tv->num_trees, sizeof(float));
    if (tv->tau == NULL) {
        tree_viewer_free(tv);
        tv = NULL;
    }
out:
    return tv;
}

int
tree_viewer_init(tree_viewer_t *self)
{
    int ret = -1;
    FILE *f = NULL;
    int j, l;
    int *pi;
    float *tau;
    coalescence_record_t cr;

    f = fopen(self->coalescence_record_filename, "r");
    if (f == NULL) {
        ret = MSP_ERR_IO;
        goto out;
    }
    /* seek to the start of file and read all records */
    if (fseek(f, 0, SEEK_SET) != 0) {
        ret = MSP_ERR_IO;
        goto out;
    }
    while ((ret = fread(&cr, sizeof(cr), 1, f)) == 1) {
        /* TODO this is very inefficient - should be using binary search */
        for (j = 0; j < self->num_trees; j++) {
            l = self->breakpoints[j];
            if (cr.left <= l && l <= cr.right) {
                pi = self->pi + j * 2 * self->sample_size;
                tau = self->tau + j * 2 * self->sample_size;
                pi[cr.children[0]] = cr.parent;
                pi[cr.children[1]] = cr.parent;
                tau[cr.parent] = cr.time;
            }
        }
    }
    if (!feof(f)) {
        ret = MSP_ERR_IO;
        goto out;
    }
    ret = 0;
out:
    if (f != NULL) {
        fclose(f);
    }
    return ret;
}

int
tree_viewer_get_tree(tree_viewer_t *self, int j, int *breakpoint, int **pi,
        float **tau)
{
    int ret = 0;
    *pi = self->pi + j * 2 * self->sample_size;
    *tau = self->tau + j * 2 * self->sample_size;
    *breakpoint = self->breakpoints[j];
    return ret;
}

int
tree_viewer_free(tree_viewer_t *self)
{
    if (self->breakpoints != NULL) {
        free(self->breakpoints);
    }
    if (self->pi != NULL) {
        free(self->pi);
    }
    if (self->tau != NULL) {
        free(self->tau);
    }
    if (self->coalescence_record_filename != NULL) {
        free(self->coalescence_record_filename);
    }
    // TODO standardise this protocol; we should always free self in these
    // destructors or never.
    free(self);
    return 0;
}

char *
msp_strerror(int err)
{
    char *ret = "Unknown error";
    if (err == MSP_ERR_NO_MEMORY) {
        ret = "Out of memory";
    } else if (err == MSP_ERR_GENERIC) {
        ret = "Generic error; please file a bug report";
    } else if (err == MSP_ERR_IO) {
        if (errno != 0) {
            ret = strerror(errno);
        } else {
            ret = "Unspecified IO error";
        }
    }
    return ret;
}
