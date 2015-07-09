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

#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <structmember.h>
#include <float.h>

#include <hdf5.h>
#include <gsl/gsl_version.h>

#include "msprime.h"

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

#define MODULE_DOC \
"Low level interface for msprime"

static PyObject *MsprimeInputError;
static PyObject *MsprimeLibraryError;

typedef struct {
    PyObject_HEAD
    msp_t *sim;
} Simulator;

typedef struct {
    PyObject_HEAD
    tree_sequence_t *tree_sequence;
} TreeSequence;

typedef struct {
    PyObject_HEAD
    sparse_tree_t *sparse_tree;
} SparseTree;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    tree_diff_iterator_t *tree_diff_iterator;
} TreeDiffIterator;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    newick_converter_t *newick_converter;
} NewickConverter;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    hapgen_t *haplotype_generator;
} HaplotypeGenerator;;

static void
handle_library_error(int err)
{
    PyErr_SetString(MsprimeLibraryError, msp_strerror(err));
}

static void
handle_input_error(int err)
{
    PyErr_SetString(MsprimeInputError, msp_strerror(err));
}

static PyObject *
convert_coalescence_record(coalescence_record_t *cr)
{
    return Py_BuildValue("III(II)d",
            (unsigned int) cr->left, (unsigned int) cr->right,
            (unsigned int) cr->node, (unsigned int) cr->children[0],
            (unsigned int) cr->children[1], cr->time);
}

/*
 * Retrieves a number value with the specified key from the specified
 * dictionary.
 */
static PyObject *
get_dict_number(PyObject *dict, const char *key_str)
{
    PyObject *ret = NULL;
    PyObject *value;
    PyObject *key = Py_BuildValue("s", key_str);
    if (!PyDict_Contains(dict, key)) {
        PyErr_Format(MsprimeInputError, "'%s' not specified", key_str);
        goto out;
    }
    value = PyDict_GetItem(dict, key);
    if (!PyNumber_Check(value)) {
        PyErr_Format(MsprimeInputError, "'%s' is not number", key_str);
        goto out;
    }
    ret = value;
out:
    Py_DECREF(key);
    return ret;
}

static PyObject *
convert_breakpoints(uint32_t *breakpoints, size_t num_breakpoints)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_int = NULL;
    size_t j;

    l = PyList_New(num_breakpoints);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < num_breakpoints; j++) {
        py_int = Py_BuildValue("I", (unsigned int) breakpoints[j]);
        if (py_int == NULL) {
            Py_DECREF(l);
            goto out;
        }
        PyList_SET_ITEM(l, j, py_int);
    }
    ret = l;
out:
    return ret;
}


/*===================================================================
 * Simulator
 *===================================================================
 */

static int
Simulator_check_sim(Simulator *self)
{
    int ret = 0;
    if (self->sim == NULL) {
        PyErr_SetString(PyExc_SystemError, "simulator not initialised");
        ret = -1;
    }
    return ret;
}

static int
Simulator_parse_population_models(Simulator *self, PyObject *py_pop_models)
{
    int ret = -1;
    Py_ssize_t j;
    double start_time, size, alpha;
    long type;
    int err;
    PyObject *item, *value;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    for (j = 0; j < PyList_Size(py_pop_models); j++) {
        item = PyList_GetItem(py_pop_models, j);
        if (!PyDict_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "not a dictionary");
            goto out;
        }
        value = get_dict_number(item, "start_time");
        if (value == NULL) {
            goto out;
        }
        start_time = PyFloat_AsDouble(value);
        value = get_dict_number(item, "type");
        if (value == NULL) {
            goto out;
        }
        type = PyLong_AsLong(value);
        if (type == POP_MODEL_CONSTANT) {
            value = get_dict_number(item, "size");
            if (value == NULL) {
                goto out;
            }
            size = PyFloat_AsDouble(value);
            err = msp_add_constant_population_model(self->sim, start_time, size);
        } else if (type == POP_MODEL_EXPONENTIAL) {
            value = get_dict_number(item, "alpha");
            if (value == NULL) {
                goto out;
            }
            alpha = PyFloat_AsDouble(value);
            err = msp_add_exponential_population_model(self->sim, start_time, alpha);
        } else {
            PyErr_SetString(MsprimeInputError,
                    "Invalid population model type");
            goto out;
        }
        if (err != 0) {
            handle_input_error(err);
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

static void
Simulator_dealloc(Simulator* self)
{
    if (self->sim != NULL) {
        msp_free(self->sim);
        PyMem_Free(self->sim);
        self->sim = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
Simulator_init(Simulator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int sim_ret;
    static char *kwlist[] = {"sample_size", "random_seed",
        "num_loci", "scaled_recombination_rate", "population_models",
        "max_memory", "avl_node_block_size", "segment_block_size",
        "node_mapping_block_size", "coalescence_record_block_size", NULL};
    PyObject *population_models = NULL;
    /* parameter defaults */
    unsigned int sample_size = 2;
    unsigned int num_loci = 1;
    unsigned long random_seed = 1;
    double scaled_recombination_rate = 0.0;
    Py_ssize_t max_memory = 10 * 1024 * 1024;
    Py_ssize_t avl_node_block_size = 10;
    Py_ssize_t segment_block_size = 10;
    Py_ssize_t node_mapping_block_size = 10;
    Py_ssize_t coalescence_record_block_size = 10;

    self->sim = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "Il|IdO!nnnnn", kwlist,
            &sample_size, &random_seed, &num_loci,
            &scaled_recombination_rate,
            &PyList_Type, &population_models,
            &max_memory, &avl_node_block_size, &segment_block_size,
            &node_mapping_block_size, &coalescence_record_block_size)) {
        goto out;
    }
    self->sim = PyMem_Malloc(sizeof(msp_t));
    if (self->sim == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    sim_ret = msp_alloc(self->sim, sample_size);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    sim_ret = msp_set_num_loci(self->sim, (uint32_t) num_loci);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    sim_ret = msp_set_random_seed(self->sim, (unsigned long) random_seed);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    sim_ret = msp_set_scaled_recombination_rate(self->sim, scaled_recombination_rate);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    sim_ret = msp_set_max_memory(self->sim, (size_t) max_memory);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    sim_ret = msp_set_avl_node_block_size(self->sim,
            (size_t) avl_node_block_size);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    sim_ret = msp_set_segment_block_size(self->sim,
            (size_t) segment_block_size);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    sim_ret = msp_set_node_mapping_block_size(self->sim,
            (size_t) node_mapping_block_size);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    sim_ret = msp_set_coalescence_record_block_size(self->sim,
            (size_t) coalescence_record_block_size);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    if (population_models != NULL) {
        if (Simulator_parse_population_models(self, population_models) != 0) {
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

static PyMemberDef Simulator_members[] = {
    {NULL}  /* Sentinel */
};



static PyObject *
Simulator_get_num_loci(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sim->num_loci);
out:
    return ret;
}

static PyObject *
Simulator_get_sample_size(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sim->sample_size);
out:
    return ret;
}

static PyObject *
Simulator_get_scaled_recombination_rate(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sim->scaled_recombination_rate);
out:
    return ret;
}


static PyObject *
Simulator_get_random_seed(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("k", self->sim->random_seed);
out:
    return ret;
}

static PyObject *
Simulator_get_max_memory(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sim->max_memory);
out:
    return ret;
}

static PyObject *
Simulator_get_segment_block_size(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sim->segment_block_size);
out:
    return ret;
}

static PyObject *
Simulator_get_avl_node_block_size(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sim->avl_node_block_size);
out:
    return ret;
}

static PyObject *
Simulator_get_node_mapping_block_size(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sim->node_mapping_block_size);
out:
    return ret;
}

static PyObject *
Simulator_get_coalescence_record_block_size(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n",
            (Py_ssize_t) self->sim->coalescence_record_block_size);
out:
    return ret;
}

static PyObject *
Simulator_get_time(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sim->time);
out:
    return ret;
}

static PyObject *
Simulator_get_num_ancestors(Simulator *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_num_ancestors(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_coancestry_events(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sim->num_ca_events);
out:
    return ret;
}

static PyObject *
Simulator_get_num_recombination_events(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sim->num_re_events);
out:
    return ret;
}

static PyObject *
Simulator_get_num_avl_node_blocks(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_num_avl_node_blocks(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_node_mapping_blocks(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_num_node_mapping_blocks(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_segment_blocks(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_num_segment_blocks(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_coalescence_record_blocks(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n",
            (Py_ssize_t) msp_get_num_coalescence_record_blocks(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_used_memory(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_used_memory(self->sim));
out:
    return ret;
}


static PyObject *
Simulator_get_num_breakpoints(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_num_breakpoints(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_coalescence_records(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n",
            (Py_ssize_t) msp_get_num_coalescence_records(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_individual_to_python(Simulator *self, segment_t *ind)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *t = NULL;
    size_t num_segments, j;
    segment_t *u;

    num_segments = 0;
    u = ind;
    while (u != NULL) {
        num_segments++;
        u = u->next;
    }
    l = PyList_New(num_segments);
    if (l == NULL) {
        goto out;
    }
    u = ind;
    j = 0;
    while (u != NULL) {
        t = Py_BuildValue("(I,I,i)", u->left, u->right, u->value);
        if (t == NULL) {
            Py_DECREF(l);
            goto out;
        }
        PyList_SET_ITEM(l, j, t);
        j++;
        u = u->next;
    }
    ret = l;
out:
    return ret;
}

static PyObject *
Simulator_get_ancestors(Simulator *self)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_ind = NULL;
    segment_t **ancestors = NULL;
    size_t num_ancestors, j;
    int err;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    num_ancestors = msp_get_num_ancestors(self->sim);
    ancestors = PyMem_Malloc(num_ancestors * sizeof(segment_t *));
    if (ancestors == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = msp_get_ancestors(self->sim, ancestors);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    l = PyList_New(num_ancestors);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < num_ancestors; j++) {
        py_ind = Simulator_individual_to_python(self, ancestors[j]);
        if (py_ind == NULL) {
            Py_DECREF(l);
            goto out;
        }
        PyList_SET_ITEM(l, j, py_ind);
    }
    ret = l;
out:
    if (ancestors != NULL) {
        PyMem_Free(ancestors);
    }
    return ret;
}

static PyObject *
Simulator_get_breakpoints(Simulator *self, PyObject *args)
{
    PyObject *ret = NULL;
    uint32_t *breakpoints = NULL;
    size_t num_breakpoints;
    int err;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    num_breakpoints = msp_get_num_breakpoints(self->sim);
    breakpoints = PyMem_Malloc(num_breakpoints * sizeof(uint32_t));
    if (breakpoints == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = msp_get_breakpoints(self->sim, breakpoints);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = convert_breakpoints(breakpoints, num_breakpoints);
out:
    if (breakpoints != NULL) {
        PyMem_Free(breakpoints);
    }
    return ret;
}

static PyObject *
Simulator_get_coalescence_records(Simulator *self, PyObject *args)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_cr = NULL;
    coalescence_record_t *coalescence_records = NULL;
    coalescence_record_t *cr;
    size_t num_coalescence_records, j;
    int err;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    num_coalescence_records = msp_get_num_coalescence_records(self->sim);
    coalescence_records = PyMem_Malloc(
            num_coalescence_records * sizeof(coalescence_record_t));
    if (coalescence_records == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = msp_get_coalescence_records(self->sim, coalescence_records);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    l = PyList_New(num_coalescence_records);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < num_coalescence_records; j++) {
        cr = &coalescence_records[j];
        py_cr = convert_coalescence_record(cr);
        if (py_cr == NULL) {
            Py_DECREF(l);
            goto out;
        }
        PyList_SET_ITEM(l, j, py_cr);
    }
    ret = l;
out:
    if (coalescence_records != NULL) {
        PyMem_Free(coalescence_records);
    }
    return ret;
}

static PyObject *
Simulator_get_population_models(Simulator *self)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *d = NULL;
    population_model_t *models = NULL;
    population_model_t *m;
    const char *param_name;
    size_t j = 0;
    size_t num_models;
    int sim_ret = 0;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    num_models = msp_get_num_population_models(self->sim);
    models = PyMem_Malloc(num_models * sizeof(population_model_t));
    if (models == NULL) {
        ret = PyErr_NoMemory();
        goto out;
    }
    sim_ret = msp_get_population_models(self->sim, models);
    if (sim_ret != 0) {
        handle_library_error(sim_ret);
        goto out;
    }
    l = PyList_New(num_models);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < num_models; j++) {
        m = &models[j];
        if (m->type == POP_MODEL_CONSTANT) {
            param_name = "size";
        } else if (m->type == POP_MODEL_EXPONENTIAL) {
            param_name = "alpha";
        } else {
            PyErr_SetString(PyExc_SystemError, "Unexpected pop model");
            goto out;
        }
        d = Py_BuildValue("{s:I,s:d,s:d}", "type", m->type, "start_time",
                m->start_time, param_name, m->param);
        if (d == NULL) {
            goto out;
        }
        PyList_SET_ITEM(l, j, d);
    }
    ret = l;
    l = NULL;
out:
    Py_XDECREF(l);
    if (models != NULL) {
        PyMem_Free(models);
    }
    return ret;
}

static PyObject *
Simulator_run(Simulator *self, PyObject *args)
{
    PyObject *ret = NULL;
    int status, not_done, coalesced;
    uint64_t chunk = 1024;
    double max_time = DBL_MAX;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "|d", &max_time)) {
        goto out;
    }
    not_done = 1;
    while (not_done) {
        Py_BEGIN_ALLOW_THREADS
        status = msp_run(self->sim, max_time, chunk);
        Py_END_ALLOW_THREADS
        if (status < 0) {
            handle_library_error(status);
            goto out;
        }
        not_done = status == 1;
        if (PyErr_CheckSignals() < 0) {
            goto out;
        }
    }
    coalesced = status == 0;
    /* return True if complete coalescence has occured */
    ret = coalesced ? Py_True : Py_False;
    Py_INCREF(ret);
out:
    return ret;
}

static PyObject *
Simulator_run_event(Simulator *self, PyObject *args)
{
    PyObject *ret = NULL;
    int status, coalesced;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    status = msp_run(self->sim, DBL_MAX, 1);
    if (status < 0) {
        handle_library_error(status);
        goto out;
    }
    coalesced = status == 0;
    /* return True if complete coalescence has occured */
    ret = coalesced ? Py_True : Py_False;
    Py_INCREF(ret);
out:
    return ret;
}


static PyMethodDef Simulator_methods[] = {
    {"get_num_loci", (PyCFunction) Simulator_get_num_loci, METH_NOARGS,
            "Returns the number of loci" },
    {"get_sample_size", (PyCFunction) Simulator_get_sample_size, METH_NOARGS,
            "Returns the sample size" },
    {"get_scaled_recombination_rate",
            (PyCFunction) Simulator_get_scaled_recombination_rate, METH_NOARGS,
            "Returns the scaled recombination rate." },
    {"get_random_seed", (PyCFunction) Simulator_get_random_seed, METH_NOARGS,
            "Returns the random seed" },
    {"get_max_memory", (PyCFunction) Simulator_get_max_memory, METH_NOARGS,
            "Returns the maximum memory used by the simulator" },
    {"get_segment_block_size",
            (PyCFunction) Simulator_get_segment_block_size, METH_NOARGS,
            "Returns segment block size." },
    {"get_avl_node_block_size",
            (PyCFunction) Simulator_get_avl_node_block_size, METH_NOARGS,
            "Returns avl_node block size" },
    {"get_node_mapping_block_size",
            (PyCFunction) Simulator_get_node_mapping_block_size, METH_NOARGS,
            "Returns node_mapping block size" },
    {"get_coalescence_record_block_size",
            (PyCFunction) Simulator_get_coalescence_record_block_size,
            METH_NOARGS, "Returns the coalescent record block size" },
    {"get_time", (PyCFunction) Simulator_get_time, METH_NOARGS,
            "Returns the current simulation time" },
    {"get_num_ancestors", (PyCFunction) Simulator_get_num_ancestors, METH_NOARGS,
            "Returns the number of ancestors" },
    {"get_num_coancestry_events",
            (PyCFunction) Simulator_get_num_coancestry_events, METH_NOARGS,
            "Returns the number of coancestry_events" },
    {"get_num_recombination_events",
            (PyCFunction) Simulator_get_num_recombination_events, METH_NOARGS,
            "Returns the number of recombination_events" },
    {"get_num_avl_node_blocks",
            (PyCFunction) Simulator_get_num_avl_node_blocks, METH_NOARGS,
            "Returns the number of avl_node memory blocks"},
    {"get_num_node_mapping_blocks",
            (PyCFunction) Simulator_get_num_node_mapping_blocks, METH_NOARGS,
            "Returns the number of node_mapping memory blocks"},
    {"get_num_segment_blocks",
            (PyCFunction) Simulator_get_num_segment_blocks, METH_NOARGS,
            "Returns the number of segment memory blocks"},
    {"get_num_coalescence_record_blocks",
            (PyCFunction) Simulator_get_num_coalescence_record_blocks, METH_NOARGS,
            "Returns the number of coalescence record memory blocks"},
    {"get_num_breakpoints", (PyCFunction) Simulator_get_num_breakpoints,
            METH_NOARGS, "Returns the number of recombination breakpoints" },
    {"get_num_coalescence_records",
            (PyCFunction) Simulator_get_num_coalescence_records,
            METH_NOARGS, "Returns the number of coalescence records" },
    {"get_used_memory", (PyCFunction) Simulator_get_used_memory,
            METH_NOARGS, "Returns the approximate amount of memory used." },
    {"get_ancestors", (PyCFunction) Simulator_get_ancestors, METH_NOARGS,
            "Returns the ancestors" },
    {"get_breakpoints", (PyCFunction) Simulator_get_breakpoints,
            METH_NOARGS, "Returns the list of breakpoints." },
    {"get_coalescence_records", (PyCFunction) Simulator_get_coalescence_records,
            METH_NOARGS, "Returns the coalescence records." },
    {"get_population_models", (PyCFunction) Simulator_get_population_models,
            METH_VARARGS, "Returns the population models"},
    {"run", (PyCFunction) Simulator_run, METH_VARARGS,
            "Simulates until at most the specified time. Returns True\
            if sample has coalesced and False otherwise." },
    {"run_event", (PyCFunction) Simulator_run_event, METH_NOARGS,
            "Simulates exactly one event. Returns True\
            if sample has coalesced and False otherwise." },
    {NULL}  /* Sentinel */
};


static PyTypeObject SimulatorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.Simulator",             /* tp_name */
    sizeof(Simulator),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Simulator_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "Simulator objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    Simulator_methods,             /* tp_methods */
    Simulator_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Simulator_init,      /* tp_init */
};

/*===================================================================
 * TreeSequence
 *===================================================================
 */

static int
TreeSequence_check_tree_sequence(TreeSequence *self)
{
    int ret = 0;
    if (self->tree_sequence == NULL) {
        PyErr_SetString(PyExc_ValueError, "tree_sequence not initialised");
        ret = -1;
    }
    return ret;
}

static void
TreeSequence_dealloc(TreeSequence* self)
{
    if (self->tree_sequence != NULL) {
        tree_sequence_free(self->tree_sequence);
        PyMem_Free(self->tree_sequence);
        self->tree_sequence = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
TreeSequence_init(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    self->tree_sequence = NULL;
    return 0;
}

static PyObject *
TreeSequence_create(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Simulator *sim = NULL;

    if (self->tree_sequence != NULL) {
        PyErr_SetString(PyExc_ValueError, "tree_sequence already created");
        goto out;
    }
    if (!PyArg_ParseTuple(args, "O!", &SimulatorType, &sim)) {
        goto out;
    }
    if (Simulator_check_sim(sim) != 0) {
        goto out;
    }
    self->tree_sequence = PyMem_Malloc(sizeof(tree_sequence_t));
    if (self->tree_sequence == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->tree_sequence, 0, sizeof(tree_sequence_t));
    err = tree_sequence_create(self->tree_sequence, sim->sim);
    if (err != 0) {
        PyMem_Free(self->tree_sequence);
        self->tree_sequence = NULL;
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
TreeSequence_dump(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    int err;
    char *path;
    PyObject *ret = NULL;
    int zlib_compression = 0;
    int skip_h5close = 0;
    int flags = 0;
    static char *kwlist[] = {"path", "zlib_compression", "skip_h5close",
        NULL};

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|ii", kwlist,
                &path, &zlib_compression, &skip_h5close)) {
        goto out;
    }
    if (zlib_compression) {
        flags = MSP_ZLIB_COMPRESSION;
    }
    if (skip_h5close) {
        flags |= MSP_SKIP_H5CLOSE;
    }
    /* Silence the low-level error reporting HDF5 */
    if (H5Eset_auto(H5E_DEFAULT, NULL, NULL) < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Error silencing HDF5 errors");
        goto out;
    }
    err = tree_sequence_dump(self->tree_sequence, path, flags);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
TreeSequence_load(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    int err;
    char *path;
    int flags = 0;
    int skip_h5close = 0;
    PyObject *ret = NULL;
    static char *kwlist[] = {"path", "skip_h5close", NULL};

    if (self->tree_sequence != NULL) {
        PyErr_SetString(PyExc_ValueError, "TreeSequence already initialised");
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|i", kwlist,
                &path, &skip_h5close)) {
        goto out;
    }
    self->tree_sequence = PyMem_Malloc(sizeof(tree_sequence_t));
    if (self->tree_sequence == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->tree_sequence, 0, sizeof(tree_sequence_t));
    if (skip_h5close) {
        flags |= MSP_SKIP_H5CLOSE;
    }
    /* Silence the low-level error reporting HDF5 */
    if (H5Eset_auto(H5E_DEFAULT, NULL, NULL) < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Error silencing HDF5 errors");
        goto out;
    }
    err = tree_sequence_load(self->tree_sequence, path, flags);
    if (err != 0) {
        PyMem_Free(self->tree_sequence);
        self->tree_sequence = NULL;
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
TreeSequence_generate_mutations(TreeSequence *self,
        PyObject *args, PyObject *kwds)
{
    int err;
    PyObject *ret = NULL;
    static char *kwlist[] = {"mutation_rate", "random_seed", NULL};
    double mutation_rate;
    unsigned long random_seed;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dk", kwlist,
            &mutation_rate, &random_seed)) {
        goto out;
    }
    err = tree_sequence_generate_mutations(self->tree_sequence,
            mutation_rate, random_seed);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
TreeSequence_get_record(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    size_t num_records;
    int order = MSP_ORDER_TIME;
    Py_ssize_t record_index;
    coalescence_record_t cr;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n|i", &record_index, &order)) {
        goto out;
    }
    num_records = tree_sequence_get_num_coalescence_records(self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tree_sequence_get_record(self->tree_sequence,
            (size_t) record_index, &cr, order);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = convert_coalescence_record(&cr);
out:
    return ret;
}

static PyObject *
TreeSequence_get_mutations(TreeSequence *self, PyObject *args)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_mutation = NULL;
    uint32_t *nodes = NULL;
    double *positions = NULL;
    size_t j, num_mutations;
    int err;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_mutations = tree_sequence_get_num_mutations(self->tree_sequence);
    nodes = PyMem_Malloc(num_mutations * sizeof(uint32_t));
    if (nodes == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    positions = PyMem_Malloc(num_mutations * sizeof(double));
    if (positions == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tree_sequence_get_mutations(self->tree_sequence, nodes, positions);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    l = PyList_New(num_mutations);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < num_mutations; j++) {
        py_mutation = Py_BuildValue("Id", (unsigned int) nodes[j],
                positions[j]);
        if (py_mutation == NULL) {
            Py_DECREF(l);
            goto out;
        }
        PyList_SET_ITEM(l, j, py_mutation);
    }
    ret = l;
out:
    if (nodes != NULL) {
        PyMem_Free(nodes);
    }
    if (positions != NULL) {
        PyMem_Free(positions);
    }
    return ret;
}

static PyObject *
TreeSequence_get_num_records(TreeSequence *self, PyObject *args)
{
    PyObject *ret = NULL;
    size_t num_records;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_records = tree_sequence_get_num_coalescence_records(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_records);
out:
    return ret;
}

static PyObject *
TreeSequence_get_num_loci(TreeSequence  *self)
{
    PyObject *ret = NULL;
    size_t num_loci;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_loci = tree_sequence_get_num_loci(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_loci);
out:
    return ret;
}

static PyObject *
TreeSequence_get_sample_size(TreeSequence  *self)
{
    PyObject *ret = NULL;
    size_t sample_size;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    sample_size = tree_sequence_get_sample_size(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) sample_size);
out:
    return ret;
}

static PyObject *
TreeSequence_get_num_nodes(TreeSequence  *self)
{
    PyObject *ret = NULL;
    size_t num_nodes;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_nodes = tree_sequence_get_num_nodes(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_nodes);
out:
    return ret;
}


static PyObject *
TreeSequence_get_num_mutations(TreeSequence  *self)
{
    PyObject *ret = NULL;
    size_t num_mutations;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_mutations = tree_sequence_get_num_mutations(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_mutations);
out:
    return ret;
}

static PyObject *
TreeSequence_get_simulation_parameters(TreeSequence  *self)
{
    PyObject *ret = NULL;
    char *str = NULL;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    str = tree_sequence_get_simulation_parameters(self->tree_sequence);
    ret = Py_BuildValue("s", str);
out:
    return ret;
}

static PyObject *
TreeSequence_get_mutation_parameters(TreeSequence  *self)
{
    PyObject *ret = NULL;
    char *str = NULL;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    str = tree_sequence_get_mutation_parameters(self->tree_sequence);
    ret = Py_BuildValue("s", str);
out:
    return ret;
}


static PyMemberDef TreeSequence_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef TreeSequence_methods[] = {
    {"create", (PyCFunction) TreeSequence_create, METH_VARARGS,
        "Creates a new TreeSequence from the specified simulator."},
    {"dump", (PyCFunction) TreeSequence_dump,
        METH_VARARGS|METH_KEYWORDS,
        "Writes the tree sequence out to the specified path."},
    {"load", (PyCFunction) TreeSequence_load,
        METH_VARARGS|METH_KEYWORDS,
        "Loads a tree sequence from the specified path."},
    {"generate_mutations", (PyCFunction) TreeSequence_generate_mutations,
        METH_VARARGS|METH_KEYWORDS,
        "Generates mutations under the infinite sites model"},
    {"get_record", (PyCFunction) TreeSequence_get_record, METH_VARARGS,
        "Returns the record at the specified index."},
    {"get_mutations", (PyCFunction) TreeSequence_get_mutations,
        METH_NOARGS, "Returns the list of mutations"},
    {"get_num_records", (PyCFunction) TreeSequence_get_num_records,
            METH_NOARGS, "Returns the number of coalescence records." },
    {"get_num_loci", (PyCFunction) TreeSequence_get_num_loci, METH_NOARGS,
            "Returns the number of loci" },
    {"get_num_mutations", (PyCFunction) TreeSequence_get_num_mutations, METH_NOARGS,
            "Returns the number of loci" },
    {"get_num_nodes", (PyCFunction) TreeSequence_get_num_nodes, METH_NOARGS,
            "Returns the number of unique nodes in the tree sequence." },
    {"get_sample_size", (PyCFunction) TreeSequence_get_sample_size, METH_NOARGS,
            "Returns the sample size" },
    {"get_simulation_parameters",
            (PyCFunction) TreeSequence_get_simulation_parameters, METH_NOARGS,
            "Returns the simulation parameters encoded as JSON." },
    {"get_mutation_parameters",
            (PyCFunction) TreeSequence_get_mutation_parameters, METH_NOARGS,
            "Returns the mutation parameters encoded as JSON." },
    {NULL}  /* Sentinel */
};

static PyTypeObject TreeSequenceType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.TreeSequence",             /* tp_name */
    sizeof(TreeSequence),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)TreeSequence_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "TreeSequence objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    TreeSequence_methods,             /* tp_methods */
    TreeSequence_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)TreeSequence_init,      /* tp_init */
};

/*===================================================================
 * SparseTree
 *===================================================================
 */

static int
SparseTree_check_sparse_tree(SparseTree *self)
{
    int ret = 0;
    if (self->sparse_tree == NULL) {
        PyErr_SetString(PyExc_ValueError, "sparse_tree not initialised");
        ret = -1;
    }
    return ret;
}

static int
SparseTree_check_bounds(SparseTree *self, unsigned int node)
{
    int ret = 0;
    if (node > self->sparse_tree->num_nodes) {
        PyErr_SetString(PyExc_ValueError, "Node index out of bounds");
        ret = -1;
    }
    return ret;
}

static void
SparseTree_dealloc(SparseTree* self)
{
    if (self->sparse_tree != NULL) {
        sparse_tree_free(self->sparse_tree);
        PyMem_Free(self->sparse_tree);
        self->sparse_tree = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
SparseTree_init(SparseTree *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"num_nodes", NULL};
    Py_ssize_t num_nodes;

    self->sparse_tree = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "n", kwlist, &num_nodes)) {
        goto out;
    }
    self->sparse_tree = PyMem_Malloc(sizeof(sparse_tree_t));
    if (self->sparse_tree == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = sparse_tree_alloc(self->sparse_tree, (size_t) num_nodes);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    err = sparse_tree_clear(self->sparse_tree);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}




static PyObject *
SparseTree_get_sample_size(SparseTree  *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sparse_tree->sample_size);
out:
    return ret;
}

static PyObject *
SparseTree_get_num_nodes(SparseTree  *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sparse_tree->num_nodes);
out:
    return ret;
}

static PyObject *
SparseTree_get_root(SparseTree  *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sparse_tree->root);
out:
    return ret;
}

static PyObject *
SparseTree_get_left(SparseTree  *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sparse_tree->left);
out:
    return ret;
}


static PyObject *
SparseTree_get_right(SparseTree  *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sparse_tree->right);
out:
    return ret;
}

static PyObject *
SparseTree_get_parent(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    unsigned int parent;
    unsigned int node;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &node)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, node)) {
        goto out;
    }
    parent = self->sparse_tree->parent[node];
    ret = Py_BuildValue("I", parent);
out:
    return ret;
}

static PyObject *
SparseTree_get_time(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    double time;
    unsigned int node;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &node)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, node)) {
        goto out;
    }
    time = self->sparse_tree->time[node];
    ret = Py_BuildValue("d", time);
out:
    return ret;
}

static PyObject *
SparseTree_get_children(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    unsigned int children[2];
    unsigned int node;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &node)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, node)) {
        goto out;
    }
    children[0] = self->sparse_tree->children[2 * node];
    children[1] = self->sparse_tree->children[2 * node + 1];
    ret = Py_BuildValue("II", children[0], children[1]);
out:
    return ret;
}



static PyMemberDef SparseTree_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef SparseTree_methods[] = {
    {"get_num_nodes", (PyCFunction) SparseTree_get_num_nodes, METH_NOARGS,
            "Returns the number of nodes in the sparse tree." },
    {"get_sample_size", (PyCFunction) SparseTree_get_sample_size, METH_NOARGS,
            "Returns the sample size" },
    {"get_root", (PyCFunction) SparseTree_get_root, METH_NOARGS,
            "Returns the root of the tree." },
    {"get_left", (PyCFunction) SparseTree_get_left, METH_NOARGS,
            "Returns the left-most coordinate (inclusive)." },
    {"get_right", (PyCFunction) SparseTree_get_right, METH_NOARGS,
            "Returns the right-most coordinate (exclusive)." },
    {"get_parent", (PyCFunction) SparseTree_get_parent, METH_VARARGS,
            "Returns the parent of node u" },
    {"get_time", (PyCFunction) SparseTree_get_time, METH_VARARGS,
            "Returns the time of node u" },
    {"get_children", (PyCFunction) SparseTree_get_children, METH_VARARGS,
            "Returns the children of node u" },
    {NULL}  /* Sentinel */
};

static PyTypeObject SparseTreeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.SparseTree",             /* tp_name */
    sizeof(SparseTree),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)SparseTree_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "SparseTree objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    SparseTree_methods,             /* tp_methods */
    SparseTree_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)SparseTree_init,      /* tp_init */
};



/*===================================================================
 * TreeDiffIterator
 *===================================================================
 */

static int
TreeDiffIterator_check_state(TreeDiffIterator *self)
{
    int ret = 0;
    if (self->tree_diff_iterator == NULL) {
        PyErr_SetString(PyExc_SystemError, "iterator not initialised");
        ret = -1;
    }
    return ret;
}

static void
TreeDiffIterator_dealloc(TreeDiffIterator* self)
{
    if (self->tree_diff_iterator != NULL) {
        tree_diff_iterator_free(self->tree_diff_iterator);
        PyMem_Free(self->tree_diff_iterator);
        self->tree_diff_iterator = NULL;
    }
    Py_XDECREF(self->tree_sequence);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
TreeDiffIterator_init(TreeDiffIterator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"tree_sequence", NULL};
    TreeSequence *tree_sequence;

    self->tree_diff_iterator = NULL;
    self->tree_sequence = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
            &TreeSequenceType, &tree_sequence)) {
        goto out;
    }
    self->tree_sequence = tree_sequence;
    Py_INCREF(self->tree_sequence);
    if (TreeSequence_check_tree_sequence(self->tree_sequence) != 0) {
        goto out;
    }
    self->tree_diff_iterator = PyMem_Malloc(sizeof(tree_diff_iterator_t));
    if (self->tree_diff_iterator == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->tree_diff_iterator, 0, sizeof(tree_diff_iterator_t));
    err = tree_diff_iterator_alloc(self->tree_diff_iterator,
            self->tree_sequence->tree_sequence);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
TreeDiffIterator_next(TreeDiffIterator  *self)
{
    PyObject *ret = NULL;
    PyObject *out_list = NULL;
    PyObject *in_list = NULL;
    PyObject *value = NULL;
    int err;
    uint32_t length;
    size_t list_size, j;
    tree_node_t *nodes_out, *nodes_in, *node;

    if (TreeDiffIterator_check_state(self) != 0) {
        goto out;
    }
    err = tree_diff_iterator_next(self->tree_diff_iterator, &length,
            &nodes_out, &nodes_in);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    if (err == 1) {
        /* out nodes */
        node = nodes_out;
        list_size = 0;
        while (node != NULL) {
            list_size++;
            node = node->next;
        }
        out_list = PyList_New(list_size);
        if (out_list == NULL) {
            goto out;
        }
        node = nodes_out;
        j = 0;
        while (node != NULL) {
            value = Py_BuildValue("I(II)d", (unsigned int) node->id,
                    (unsigned int) node->children[0],
                    (unsigned int) node->children[1], node->time);
            if (value == NULL) {
                goto out;
            }
            PyList_SET_ITEM(out_list, j, value);
            node = node->next;
            j++;
        }
        /* in nodes */
        node = nodes_in;
        list_size = 0;
        while (node != NULL) {
            list_size++;
            node = node->next;
        }
        in_list = PyList_New(list_size);
        if (in_list == NULL) {
            goto out;
        }
        node = nodes_in;
        j = 0;
        while (node != NULL) {
            value = Py_BuildValue("I(II)d", (unsigned int) node->id,
                    (unsigned int) node->children[0],
                    (unsigned int) node->children[1], node->time);
            if (value == NULL) {
                goto out;
            }
            PyList_SET_ITEM(in_list, j, value);
            node = node->next;
            j++;
        }
        ret = Py_BuildValue("IOO", (unsigned int) length, out_list, in_list);
    }
out:
    Py_XDECREF(out_list);
    Py_XDECREF(in_list);
    return ret;
}

static PyMemberDef TreeDiffIterator_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef TreeDiffIterator_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject TreeDiffIteratorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.TreeDiffIterator",             /* tp_name */
    sizeof(TreeDiffIterator),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)TreeDiffIterator_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "TreeDiffIterator objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    PyObject_SelfIter,                    /* tp_iter */
    (iternextfunc) TreeDiffIterator_next, /* tp_iternext */
    TreeDiffIterator_methods,             /* tp_methods */
    TreeDiffIterator_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)TreeDiffIterator_init,      /* tp_init */
};

/*===================================================================
 * NewickConverter
 *===================================================================
 */

static int
NewickConverter_check_state(NewickConverter *self)
{
    int ret = 0;
    if (self->newick_converter == NULL) {
        PyErr_SetString(PyExc_SystemError, "converter not initialised");
        ret = -1;
    }
    return ret;
}

static void
NewickConverter_dealloc(NewickConverter* self)
{
    if (self->newick_converter != NULL) {
        newick_converter_free(self->newick_converter);
        PyMem_Free(self->newick_converter);
        self->newick_converter = NULL;
    }
    Py_XDECREF(self->tree_sequence);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
NewickConverter_init(NewickConverter *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"tree_sequence", "precision", NULL};
    int precision = 3;
    TreeSequence *tree_sequence;

    self->newick_converter = NULL;
    self->tree_sequence = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|i", kwlist,
            &TreeSequenceType, &tree_sequence, &precision)) {
        goto out;
    }
    self->tree_sequence = tree_sequence;
    Py_INCREF(self->tree_sequence);
    if (TreeSequence_check_tree_sequence(self->tree_sequence) != 0) {
        goto out;
    }
    if (precision < 0 || precision > 16) {
        PyErr_SetString(PyExc_ValueError,
                "precision value out of range (0, 16)");
        goto out;
    }
    self->newick_converter = PyMem_Malloc(sizeof(newick_converter_t));
    if (self->newick_converter == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->newick_converter, 0, sizeof(newick_converter_t));
    err = newick_converter_alloc(self->newick_converter,
            self->tree_sequence->tree_sequence, (size_t) precision);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
NewickConverter_next(NewickConverter  *self)
{
    PyObject *ret = NULL;
    uint32_t length;
    char *tree;
    int err;

    if (NewickConverter_check_state(self) != 0) {
        goto out;
    }
    err = newick_converter_next(self->newick_converter, &length, &tree);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    if (err == 1) {
        ret = Py_BuildValue("Is", (unsigned int) length, tree);
    }
out:
    return ret;
}

static PyMemberDef NewickConverter_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef NewickConverter_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject NewickConverterType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.NewickConverter",             /* tp_name */
    sizeof(NewickConverter),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)NewickConverter_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "NewickConverter objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    PyObject_SelfIter,                    /* tp_iter */
    (iternextfunc) NewickConverter_next, /* tp_iternext */
    NewickConverter_methods,             /* tp_methods */
    NewickConverter_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)NewickConverter_init,      /* tp_init */
};

/*===================================================================
 * HaplotypeGenerator
 *===================================================================
 */

static int
HaplotypeGenerator_check_state(HaplotypeGenerator *self)
{
    int ret = 0;
    if (self->haplotype_generator == NULL) {
        PyErr_SetString(PyExc_SystemError, "converter not initialised");
        ret = -1;
    }
    return ret;
}

static void
HaplotypeGenerator_dealloc(HaplotypeGenerator* self)
{
    if (self->haplotype_generator != NULL) {
        hapgen_free(self->haplotype_generator);
        PyMem_Free(self->haplotype_generator);
        self->haplotype_generator = NULL;
    }
    Py_XDECREF(self->tree_sequence);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
HaplotypeGenerator_init(HaplotypeGenerator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"tree_sequence", NULL};
    TreeSequence *tree_sequence;

    self->haplotype_generator = NULL;
    self->tree_sequence = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
            &TreeSequenceType, &tree_sequence)) {
        goto out;
    }
    self->tree_sequence = tree_sequence;
    Py_INCREF(self->tree_sequence);
    if (TreeSequence_check_tree_sequence(self->tree_sequence) != 0) {
        goto out;
    }
    self->haplotype_generator = PyMem_Malloc(sizeof(hapgen_t));
    if (self->haplotype_generator == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->haplotype_generator, 0, sizeof(hapgen_t));
    err = hapgen_alloc(self->haplotype_generator,
            self->tree_sequence->tree_sequence);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
HaplotypeGenerator_get_haplotype(HaplotypeGenerator *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    char *haplotype;
    unsigned int sample_id;

    if (HaplotypeGenerator_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &sample_id)) {
        goto out;
    }
    err = hapgen_get_haplotype(self->haplotype_generator,
            (uint32_t) sample_id, &haplotype);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("s", haplotype);
out:
    return ret;
}

static PyMemberDef HaplotypeGenerator_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef HaplotypeGenerator_methods[] = {
    {"get_haplotype", (PyCFunction) HaplotypeGenerator_get_haplotype,
        METH_VARARGS, "Returns the haplotype for the specified sample"},
    {NULL}  /* Sentinel */
};

static PyTypeObject HaplotypeGeneratorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.HaplotypeGenerator",             /* tp_name */
    sizeof(HaplotypeGenerator),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)HaplotypeGenerator_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "HaplotypeGenerator objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    HaplotypeGenerator_methods,             /* tp_methods */
    HaplotypeGenerator_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)HaplotypeGenerator_init,      /* tp_init */
};

/*===================================================================
 * Module level functions
 *===================================================================
 */

static PyObject *
msprime_get_gsl_version(PyObject *self)
{
    return Py_BuildValue("ii", GSL_MAJOR_VERSION, GSL_MINOR_VERSION);
}

static PyObject *
msprime_get_hdf5_version(PyObject *self)
{
    herr_t status;
    PyObject *ret = NULL;
    unsigned int major, minor, release;

    status = H5get_libversion(&major, &minor, &release);
    if (status != 0) {
        PyErr_SetString(PyExc_SystemError, "Error getting HDF5 version");
        goto out;
    }
    ret = Py_BuildValue("III", major, minor, release);
out:
    return ret;
}

static PyMethodDef msprime_methods[] = {
    {"get_gsl_version", (PyCFunction) msprime_get_gsl_version, METH_NOARGS,
            "Returns the version of GSL we are linking against." },
    {"get_hdf5_version", (PyCFunction) msprime_get_hdf5_version, METH_NOARGS,
            "Returns the version of HDF5 we are linking against." },
    {NULL}        /* Sentinel */
};

/* Initialisation code supports Python 2.x and 3.x. The framework uses the
 * recommended structure from http://docs.python.org/howto/cporting.html.
 * I've ignored the point about storing state in globals, as the examples
 * from the Python documentation still use this idiom.
 */

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef msprimemodule = {
    PyModuleDef_HEAD_INIT,
    "_msprime",   /* name of module */
    MODULE_DOC, /* module documentation, may be NULL */
    -1,
    msprime_methods,
    NULL, NULL, NULL, NULL
};

#define INITERROR return NULL

PyObject *
PyInit__msprime(void)

#else
#define INITERROR return

void
init_msprime(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&msprimemodule);
#else
    PyObject *module = Py_InitModule3("_msprime", msprime_methods, MODULE_DOC);
#endif
    if (module == NULL) {
        INITERROR;
    }
    /* Simulator type */
    SimulatorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&SimulatorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&SimulatorType);
    PyModule_AddObject(module, "Simulator", (PyObject *) &SimulatorType);
    /* TreeSequence type */
    TreeSequenceType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&TreeSequenceType) < 0) {
        INITERROR;
    }
    Py_INCREF(&TreeSequenceType);
    PyModule_AddObject(module, "TreeSequence", (PyObject *) &TreeSequenceType);
    /* SparseTree type */
    SparseTreeType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&SparseTreeType) < 0) {
        INITERROR;
    }
    Py_INCREF(&SparseTreeType);
    PyModule_AddObject(module, "SparseTree", (PyObject *) &SparseTreeType);
    /* TreeDiffIterator type */
    TreeDiffIteratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&TreeDiffIteratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&TreeDiffIteratorType);
    PyModule_AddObject(module, "TreeDiffIterator",
            (PyObject *) &TreeDiffIteratorType);
    /* NewickConverter type */
    NewickConverterType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&NewickConverterType) < 0) {
        INITERROR;
    }
    Py_INCREF(&NewickConverterType);
    PyModule_AddObject(module, "NewickConverter",
            (PyObject *) &NewickConverterType);
    /* HaplotypeGenerator type */
    HaplotypeGeneratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&HaplotypeGeneratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&HaplotypeGeneratorType);
    PyModule_AddObject(module, "HaplotypeGenerator",
            (PyObject *) &HaplotypeGeneratorType);
     /* Errors and constants */
    MsprimeInputError = PyErr_NewException("_msprime.InputError", NULL, NULL);
    Py_INCREF(MsprimeInputError);
    PyModule_AddObject(module, "InputError", MsprimeInputError);
    MsprimeLibraryError = PyErr_NewException("_msprime.LibraryError", NULL,
            NULL);
    Py_INCREF(MsprimeLibraryError);
    PyModule_AddObject(module, "LibraryError", MsprimeLibraryError);

    PyModule_AddIntConstant(module, "POP_MODEL_CONSTANT", POP_MODEL_CONSTANT);
    PyModule_AddIntConstant(module, "POP_MODEL_EXPONENTIAL",
            POP_MODEL_EXPONENTIAL);
    PyModule_AddIntConstant(module, "MSP_ORDER_TIME", MSP_ORDER_TIME);
    PyModule_AddIntConstant(module, "MSP_ORDER_LEFT", MSP_ORDER_LEFT);
    PyModule_AddIntConstant(module, "MSP_ORDER_RIGHT", MSP_ORDER_RIGHT);
    /* turn off GSL error handler so we don't abort on memory error */
    gsl_set_error_handler_off();

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}


