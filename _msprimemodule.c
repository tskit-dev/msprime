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
    char *tree_file_name;
    msp_t *sim;
} Simulator;

typedef struct {
    PyObject_HEAD
    tree_file_t *tree_file;
} TreeFile;

typedef struct {
    PyObject_HEAD
    hapgen_t *hapgen;
} HaplotypeGenerator;


static void
handle_library_error(int err)
{
    PyErr_SetString(MsprimeLibraryError, msp_strerror(err));
}

static void
handle_input_error(const char *err)
{
    PyErr_SetString(MsprimeInputError, err);
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
            PyErr_SetString(MsprimeInputError, "not a dictionary");
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
            handle_library_error(err);
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
    if (self->tree_file_name != NULL) {
        PyMem_Free(self->tree_file_name);
        self->tree_file_name = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
Simulator_check_input(Simulator *self)
{
    int ret = -1;
    msp_t *sim = self->sim;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (sim->num_loci == 0) {
        handle_input_error("must have num_loci > 0");
        goto out;
    }
    if (sim->scaled_recombination_rate < 0) {
        handle_input_error("must have 0 <= recombination_rate");
        goto out;
    }
    if (sim->sample_size < 2) {
        handle_input_error("sample_size must be > 1");
        goto out;
    }
    if (sim->max_memory == 0) {
        handle_input_error("max_memory must be > 0");
        goto out;
    }
    if (sim->avl_node_block_size == 0) {
        handle_input_error("avl_node_block_size must be > 0");
        goto out;
    }
    if (sim->segment_block_size == 0) {
        handle_input_error("segment_block_size must be > 0");
        goto out;
    }
    if (sim->node_mapping_block_size == 0) {
        handle_input_error("node_mapping_block_size must be > 0");
        goto out;
    }
    if (strlen(sim->tree_file_name) == 0) {
        handle_input_error("Cannot use empty string as filename");
        goto out;
    }
    /* TODO more checks! */
    ret = 0;
out:
    return ret;
}

static int
Simulator_init(Simulator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int sim_ret;
    static char *kwlist[] = {"sample_size", "random_seed",
        "tree_file_name", "num_loci", "scaled_recombination_rate",
        "population_models", "max_memory", "avl_node_block_size",
        "segment_block_size", "node_mapping_block_size", NULL};
    PyObject *population_models = NULL;
    /* parameter defaults */
    unsigned int sample_size = 2;
    unsigned int num_loci = 1;
    long random_seed = 1;
    double scaled_recombination_rate = 0.0;
    Py_ssize_t max_memory = 10 * 1024 * 1024;
    Py_ssize_t avl_node_block_size = 10;
    Py_ssize_t segment_block_size = 10;
    Py_ssize_t node_mapping_block_size = 10;
    const char *tree_file_name;
    Py_ssize_t tree_file_name_len;
    msp_t *sim = PyMem_Malloc(sizeof(msp_t));

    self->sim = sim;
    self->tree_file_name = NULL;
    if (self->sim == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "Ils#|IdO!nnnn", kwlist,
            &sample_size, &random_seed, &tree_file_name,
            &tree_file_name_len,
            &num_loci, &scaled_recombination_rate,
            &PyList_Type, &population_models, &max_memory,
            &avl_node_block_size, &segment_block_size,
            &node_mapping_block_size)) {
        goto out;
    }
    memset(self->sim, 0, sizeof(msp_t));
    sim->sample_size = (uint32_t) sample_size;
    sim->num_loci = (uint32_t) num_loci;
    sim->random_seed = random_seed;
    sim->scaled_recombination_rate = scaled_recombination_rate;
    sim->max_memory = (size_t) max_memory;
    sim->avl_node_block_size = (size_t) avl_node_block_size;
    sim->segment_block_size = (size_t) segment_block_size;
    sim->node_mapping_block_size = (size_t) node_mapping_block_size;
    self->tree_file_name = PyMem_Malloc(tree_file_name_len + 1);
    if (self->tree_file_name == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    strcpy(self->tree_file_name, tree_file_name);
    sim->tree_file_name = self->tree_file_name;
    /* TODO this is very nasty and must be moved into the msprime
     * code when the refactoring is done.
     */
    sim_ret = msp_add_constant_population_model(sim, -1.0, 1.0);
    if (sim_ret != 0) {
        handle_library_error(sim_ret);
        goto out;
    }
    if (population_models != NULL) {
        if (Simulator_parse_population_models(self, population_models) != 0) {
            goto out;
        }
    }
    if (Simulator_check_input(self) != 0) {
        goto out;
    }
    sim_ret = msp_alloc(self->sim);
    if (sim_ret != 0) {
        handle_library_error(sim_ret);
        goto out;
    }
    sim_ret = msp_initialise(self->sim);
    if (sim_ret != 0) {
        handle_library_error(sim_ret);
        goto out;
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
    ret = Py_BuildValue("I", self->sim->num_loci);
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
    ret = Py_BuildValue("I", self->sim->sample_size);
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
Simulator_get_tree_file_name(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("s", self->sim->tree_file_name);
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
    ret = Py_BuildValue("I", msp_get_num_ancestors(self->sim));
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
    ret = Py_BuildValue("I", self->sim->num_ca_events);
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
    ret = Py_BuildValue("I", self->sim->num_re_events);
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
    ret = Py_BuildValue("I", msp_get_num_avl_node_blocks(self->sim));
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
    ret = Py_BuildValue("I", msp_get_num_node_mapping_blocks(self->sim));
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
    ret = Py_BuildValue("I", msp_get_num_segment_blocks(self->sim));
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
    ret = Py_BuildValue("I", msp_get_num_breakpoints(self->sim));
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
    ancestors = PyMem_Malloc(num_ancestors * sizeof(segment_t **));
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
Simulator_get_population_models(Simulator *self)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *d = NULL;
    population_model_t *m;
    const char *param_name;
    size_t j = 0;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    /* TODO this is poor API, need to abstract this somehow. */
    m = self->sim->population_models;
    while (m != NULL) {
        j++;
        m = m ->next;
    }
    l = PyList_New(j);
    if (l == NULL) {
        goto out;
    }
    m = self->sim->population_models;
    j = 0;
    while (m != NULL) {
        if (m->type == POP_MODEL_CONSTANT) {
            param_name = "size";
        } else if (m->type == POP_MODEL_EXPONENTIAL) {
            param_name = "alpha";
        } else {
            PyErr_SetString(PyExc_SystemError, "Unexpected pop model");
            goto out;
        }
        d = Py_BuildValue("{s:I,s:d,s:d}", "type", m->type, "time",
                m->start_time, param_name, m->param);
        if (d == NULL) {
            goto out;
        }
        PyList_SET_ITEM(l, j, d);
        j++;
        m = m ->next;
    }
    ret = l;
    l = NULL;
out:
    Py_XDECREF(l);
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
    if (status < 0) {
        handle_library_error(status);
        goto out;
    }
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
    {"get_time", (PyCFunction) Simulator_get_time, METH_NOARGS,
            "Returns the current simulation time" },
    {"get_tree_file_name",
            (PyCFunction) Simulator_get_tree_file_name, METH_NOARGS,
            "Returns the name of the tree file." },
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
    {"get_num_breakpoints", (PyCFunction) Simulator_get_num_breakpoints,
            METH_NOARGS, "Returns the number of trees" },
    {"get_ancestors", (PyCFunction) Simulator_get_ancestors, METH_NOARGS,
            "Returns the ancestors" },
    {"get_population_models", (PyCFunction) Simulator_get_population_models,
            METH_VARARGS, "Returns the population models"},
    {"run", (PyCFunction) Simulator_run, METH_VARARGS,
            "Simulates until at most the specified time. Returns True\
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
 * TreeFile
 *===================================================================
 */
static int
TreeFile_check_tree_file(TreeFile *self)
{
    int ret = 0;
    if (self->tree_file == NULL) {
        PyErr_SetString(PyExc_SystemError, "tree_file not initialised");
        ret = -1;
    }
    return ret;
}


static void
TreeFile_dealloc(TreeFile* self)
{
    if (self->tree_file != NULL) {
        tree_file_close(self->tree_file);
        PyMem_Free(self->tree_file);
        self->tree_file = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
TreeFile_init(TreeFile *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int tr_ret;
    static char *kwlist[] = {"tree_file_name", NULL};
    char *tree_file_name;
    char mode = 'r';

    self->tree_file = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist,
                &tree_file_name)) {
        goto out;
    }
    self->tree_file = PyMem_Malloc(sizeof(tree_file_t));
    if (self->tree_file == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    tr_ret = tree_file_open(self->tree_file, tree_file_name, mode);
    if (tr_ret != 0) {
        handle_library_error(tr_ret);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyMemberDef TreeFile_members[] = {
    {NULL}  /* Sentinel */
};

static PyObject *
TreeFile_get_num_loci(TreeFile *self)
{
    PyObject *ret = NULL;
    if (TreeFile_check_tree_file(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->tree_file->num_loci);
out:
    return ret;
}

static PyObject *
TreeFile_get_sample_size(TreeFile *self)
{
    PyObject *ret = NULL;
    if (TreeFile_check_tree_file(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->tree_file->sample_size);
out:
    return ret;
}

static PyObject *
TreeFile_get_metadata(TreeFile *self)
{
    PyObject *ret = NULL;
    if (TreeFile_check_tree_file(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("s", self->tree_file->metadata);
out:
    return ret;
}

static PyObject *
TreeFile_iscomplete(TreeFile *self)
{
    PyObject *ret = NULL;
    if (TreeFile_check_tree_file(self) != 0) {
        goto out;
    }
    ret = PyBool_FromLong((long) tree_file_iscomplete(self->tree_file));
out:
    return ret;
}

static PyObject *
TreeFile_issorted(TreeFile *self)
{
    PyObject *ret = NULL;
    if (TreeFile_check_tree_file(self) != 0) {
        goto out;
    }
    ret = PyBool_FromLong((long) tree_file_issorted(self->tree_file));
out:
    return ret;
}

static PyObject *
TreeFile_next(TreeFile *self)
{
    PyObject *ret = NULL;
    coalescence_record_t cr;
    int v;

    v = tree_file_next_record(self->tree_file, &cr);
    if (v < 0) {
        handle_library_error(v);
        goto out;
    }
    if (v == 1) {
        ret = Py_BuildValue("iiiid", cr.left, cr.children[0],
                cr.children[1], cr.parent, cr.time);
    }
out:
    return ret;
}


static PyMethodDef TreeFile_methods[] = {
    {"get_num_loci", (PyCFunction) TreeFile_get_num_loci, METH_NOARGS,
            "Returns the number of loci"},
    {"get_sample_size", (PyCFunction) TreeFile_get_sample_size, METH_NOARGS,
            "Returns the sample size"},
    {"get_metadata", (PyCFunction) TreeFile_get_metadata, METH_NOARGS,
            "Returns the simulation metadata"},
    {"iscomplete", (PyCFunction) TreeFile_iscomplete, METH_NOARGS,
            "Returns True if this file corresponds to a complete simulation."},
    {"issorted", (PyCFunction) TreeFile_issorted, METH_NOARGS,
            "Returns True if this file is sorted"},
    {NULL}  /* Sentinel */
};

static PyTypeObject TreeFileType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.TreeFile",             /* tp_name */
    sizeof(TreeFile),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)TreeFile_dealloc, /* tp_dealloc */
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
    "TreeFile objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    PyObject_SelfIter,            /* tp_iter */
    (iternextfunc) TreeFile_next, /* tp_iternext */
    TreeFile_methods,             /* tp_methods */
    TreeFile_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)TreeFile_init,      /* tp_init */
};

/*===================================================================
 * HaplotypeGenerator
 *===================================================================
 */

static int
HaplotypeGenerator_check_hapgen(HaplotypeGenerator *self)
{
    int ret = 0;
    if (self->hapgen == NULL) {
        PyErr_SetString(PyExc_SystemError, "hapgen not initialised");
        ret = -1;
    }
    return ret;
}


static void
HaplotypeGenerator_dealloc(HaplotypeGenerator* self)
{
    if (self->hapgen != NULL) {
        hapgen_free(self->hapgen);
        PyMem_Free(self->hapgen);
        self->hapgen = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
HaplotypeGenerator_init(HaplotypeGenerator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int hg_ret;
    static char *kwlist[] = {"tree_file_name", "mutation_rate", "random_seed",
        "max_haplotype_length", NULL};
    char *tree_file_name;
    Py_ssize_t max_haplotype_length;
    long random_seed;
    double mutation_rate;

    self->hapgen = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "sdln", kwlist,
                &tree_file_name, &mutation_rate, &random_seed,
                &max_haplotype_length)) {
        goto out;
    }
    self->hapgen = PyMem_Malloc(sizeof(hapgen_t));
    if (self->hapgen == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    hg_ret = hapgen_alloc(self->hapgen, mutation_rate, tree_file_name,
            random_seed, (size_t) max_haplotype_length);
    if (hg_ret != 0) {
        handle_library_error(hg_ret);
        goto out;
    }
    hg_ret = hapgen_generate(self->hapgen);
    if (hg_ret != 0) {
        handle_library_error(hg_ret);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyMemberDef HaplotypeGenerator_members[] = {
    {NULL}  /* Sentinel */
};

static PyObject *
HaplotypeGenerator_get_num_loci(HaplotypeGenerator *self)
{
    PyObject *ret = NULL;
    if (HaplotypeGenerator_check_hapgen(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->hapgen->num_loci);
out:
    return ret;
}

static PyObject *
HaplotypeGenerator_get_sample_size(HaplotypeGenerator *self)
{
    PyObject *ret = NULL;
    if (HaplotypeGenerator_check_hapgen(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->hapgen->sample_size);
out:
    return ret;
}

static PyObject *
HaplotypeGenerator_get_num_trees(HaplotypeGenerator *self)
{
    PyObject *ret = NULL;
    if (HaplotypeGenerator_check_hapgen(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->hapgen->num_trees);
out:
    return ret;
}

static PyObject *
HaplotypeGenerator_get_haplotype_length(HaplotypeGenerator *self)
{
    PyObject *ret = NULL;
    if (HaplotypeGenerator_check_hapgen(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->hapgen->haplotype_length);
out:
    return ret;
}

static PyObject *
HaplotypeGenerator_get_haplotypes(HaplotypeGenerator *self)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *h = NULL;
    int hg_ret;
    uint32_t n = self->hapgen->sample_size;
    uint64_t j;
    char **haplotypes;
    size_t s;

    if (HaplotypeGenerator_check_hapgen(self) != 0) {
        goto out;
    }
    hg_ret = hapgen_get_haplotypes(self->hapgen, &haplotypes, &s);
    if (hg_ret != 0) {
        handle_library_error(hg_ret);
        goto out;
    }
    l = PyList_New(n + 1);
    if (l == NULL) {
        goto out;
    }
    Py_INCREF(Py_None);
    PyList_SET_ITEM(l, 0, Py_None);
    for (j = 1; j <= n; j++) {
        h = PyBytes_FromStringAndSize(haplotypes[j], s);
        if (h == NULL) {
            goto out;
        }
        PyList_SET_ITEM(l, j, h);
    }
    ret = l;
    l = NULL;
out:
    Py_XDECREF(l);
    return ret;
}


static PyMethodDef HaplotypeGenerator_methods[] = {
    {"get_num_loci", (PyCFunction) HaplotypeGenerator_get_num_loci, METH_NOARGS,
            "Returns the number of loci"},
    {"get_sample_size", (PyCFunction) HaplotypeGenerator_get_sample_size, METH_NOARGS,
            "Returns the sample size"},
    {"get_num_trees", (PyCFunction) HaplotypeGenerator_get_num_trees, METH_NOARGS,
            "Returns the sample size"},
    {"get_haplotype_length",
            (PyCFunction) HaplotypeGenerator_get_haplotype_length, METH_NOARGS,
            "Returns the number of segregating sites in the generated haplotypes."},
    {"get_haplotypes", (PyCFunction) HaplotypeGenerator_get_haplotypes, METH_NOARGS,
            "Returns the generated haplotypes."},
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


/*==========================================================
 * Module level functions
 *==========================================================
 */

PyDoc_STRVAR(msprime_sort_tree_file_doc,
"Sorts the specified tree file..\n");

static PyObject *
msprime_sort_tree_file(PyObject *self, PyObject *args)
{
    PyObject *ret = NULL;
    char *filename = NULL;
    tree_file_t *tf = NULL;
    int tf_ret;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        goto out;
    }
    tf = PyMem_Malloc(sizeof(tree_file_t));
    if (tf == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    tf_ret = tree_file_open(tf, filename, 'u');
    if (tf_ret != 0) {
        handle_library_error(tf_ret);
        goto out;
    }
    if (tree_file_issorted(tf)) {
        PyErr_SetString(PyExc_ValueError, "Tree file is already sorted");
        goto out;
    }
    tf_ret = tree_file_sort(tf);
    if (tf_ret != 0) {
        handle_library_error(tf_ret);
        goto out;
    }
    tf_ret = tree_file_close(tf);
    if (tf_ret != 0) {
        handle_library_error(tf_ret);
        goto out;
    }
    ret = Py_None;
    Py_INCREF(Py_None);
out:
    if (tf != NULL) {
        if (tree_file_isopen(tf)) {
            /* we are already in an error condition, so ignore new errors */
            tree_file_close(tf);
        }
        PyMem_Free(tf);
    }
    return ret;
}

static PyMethodDef msprime_methods[] = {
    {"sort_tree_file", (PyCFunction) msprime_sort_tree_file, METH_VARARGS,
        msprime_sort_tree_file_doc},
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
    SimulatorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&SimulatorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&SimulatorType);
    PyModule_AddObject(module, "Simulator", (PyObject *) &SimulatorType);
    TreeFileType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&TreeFileType) < 0) {
        INITERROR;
    }
    Py_INCREF(&TreeFileType);
    PyModule_AddObject(module, "TreeFile", (PyObject *) &TreeFileType);
    HaplotypeGeneratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&HaplotypeGeneratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&HaplotypeGeneratorType);
    PyModule_AddObject(module, "HaplotypeGenerator",
            (PyObject *) &HaplotypeGeneratorType);
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

#ifdef WORDS_BIGENDIAN
    PyErr_Format(PyExc_RuntimeError, "Big Endian systems not currently supported.");
    INITERROR;
#endif

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}


