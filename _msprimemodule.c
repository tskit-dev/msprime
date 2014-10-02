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

#include <Python.h>
#include <structmember.h>
#include <float.h>

#include "msprime.h"

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif
#define PY_SSIZE_T_CLEAN

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
    int iter_state;
} TreeFile;

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
    double time, size, alpha;
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
        value = get_dict_number(item, "time");
        if (value == NULL) {
            goto out;
        }
        time = PyFloat_AsDouble(value);
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
            err = msp_add_constant_population_model(self->sim, time, size);
        } else if (type == POP_MODEL_EXPONENTIAL) {
            value = get_dict_number(item, "alpha");
            if (value == NULL) {
                goto out;
            }
            alpha = PyFloat_AsDouble(value);
            err = msp_add_constant_population_model(self->sim, time, alpha);
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
    if (sim->recombination_rate < 0 || sim->recombination_rate > 1) {
        handle_input_error("must have 0 <= recombination_rate <= 1");
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
        "tree_file_name", "num_loci", "recombination_rate",
        "population_models", "max_memory", "avl_node_block_size",
        "segment_block_size", "node_mapping_block_size", NULL};
    PyObject *population_models = NULL;
    msp_t *sim = PyMem_Malloc(sizeof(msp_t));
    char *cr_filename;
    Py_ssize_t cr_filename_len;

    self->tree_file_name = NULL;
    self->sim = sim;
    if (self->sim == NULL) {
        goto out;
    }
    memset(self->sim, 0, sizeof(msp_t));
    sim->sample_size = 2;
    sim->num_loci = 1;
    sim->random_seed = 1;
    sim->recombination_rate = 0.5;
    sim->max_memory = 10 * 1024 * 1024;
    sim->avl_node_block_size = 10;
    sim->segment_block_size = 10;
    sim->node_mapping_block_size = 10;
    sim->tree_file_name = NULL;
    /* TODO verify these types are compatible! */
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "Ils#|IdO!nnnn", kwlist,
            &sim->sample_size, &sim->random_seed, &cr_filename,
            &cr_filename_len,
            &sim->num_loci, &sim->recombination_rate,
            &PyList_Type, &population_models, &sim->max_memory,
            &sim->avl_node_block_size, &sim->segment_block_size,
            &sim->node_mapping_block_size)) {
        goto out;
    }
    self->tree_file_name = PyMem_Malloc(cr_filename_len + 1);
    if (self->tree_file_name == NULL) {
        goto out;
    }
    strcpy(self->tree_file_name, cr_filename);
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
Simulator_get_tree_file(Simulator  *self)
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
Simulator_get_recombination_rate(Simulator *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sim->recombination_rate);
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
    status = msp_finalise_tree_file(self->sim);
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
    {"get_random_seed", (PyCFunction) Simulator_get_random_seed, METH_NOARGS,
            "Returns the random seed" },
    {"get_time", (PyCFunction) Simulator_get_time, METH_NOARGS,
            "Returns the current simulation time" },
    {"get_tree_file", (PyCFunction) Simulator_get_tree_file, METH_NOARGS,
            "Returns the name of the tree file." },
    {"get_recombination_rate",
            (PyCFunction) Simulator_get_recombination_rate, METH_NOARGS,
            "Returns the rate of recombination between adjacent loci" },
    {"get_num_ancestors", (PyCFunction) Simulator_get_num_ancestors, METH_NOARGS,
            "Returns the number of ancestors" },
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
        tree_file_free(self->tree_file);
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
    tree_file_t *tree_file = PyMem_Malloc(sizeof(tree_file_t));

    self->iter_state = 0;
    self->tree_file = tree_file;
    if (self->tree_file == NULL) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist,
                &tree_file_name)) {
        goto out;
    }
    tr_ret = tree_file_alloc(self->tree_file, tree_file_name);
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
TreeFile_get_num_trees(TreeFile *self)
{
    PyObject *ret = NULL;
    if (TreeFile_check_tree_file(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->tree_file->num_trees);
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
TreeFile_sort(TreeFile *self)
{
    PyObject *ret = NULL;
    int err;

    if (TreeFile_check_tree_file(self) != 0) {
        goto out;
    }
    err = tree_file_sort(self->tree_file);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_None;
    Py_INCREF(ret);
out:
    return ret;
}

static PyObject *
TreeFile_next(TreeFile *self)
{
    PyObject *ret = NULL;
    coalescence_record_t cr;
    int v;

    if (self->iter_state == 0) {
        tree_file_record_iter_init(self->tree_file);
        self->iter_state = 1;
    }
    if (self->iter_state == 1) {
        v = tree_file_record_iter_next(self->tree_file, &cr);
        if (v < 0) {
            handle_library_error(v);
            goto out;
        }
        if (v == 0) {
            /* last record has been read, iteration is done */
            self->iter_state = 2;
        }
        ret = Py_BuildValue("iiiiid", cr.left, cr.right, cr.children[0],
                cr.children[1], cr.parent, cr.time);
    }
out:
    return ret;
}



#if 0
static PyObject *
TreeFile_get_tree(TreeFile *self, PyObject *args)
{
    PyObject *ret = NULL;
    PyObject *py_pi = NULL;
    PyObject *py_tau = NULL;
    int err;
    unsigned int j, n;
    uint32_t l;
    int32_t *pi;
    float *tau;

    if (TreeFile_check_tree_file(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &j)) {
        goto out;
    }
    if (j >= self->tree_file->num_trees) {
        handle_input_error("tree out of bounds");
        goto out;
    }
    err = tree_file_get_tree(self->tree_file, (uint32_t) j, &l, &pi, &tau);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    n = 2 * self->tree_file->sample_size;
    py_pi = PyList_New(n);
    if (py_pi == NULL) {
        goto out;
    }
    py_tau = PyList_New(n);
    if (py_tau == NULL) {
        goto out;
    }
    for (j = 0; j < n; j++) {
        err = PyList_SetItem(py_pi, j, PyLong_FromLong((long) pi[j]));
        if (err < 0) {
            goto out;
        }
        err = PyList_SetItem(py_tau, j, PyFloat_FromDouble((double) tau[j]));
        if (err < 0) {
            goto out;
        }
    }
    ret = Py_BuildValue("(I, O, O)", (unsigned int) l, py_pi, py_tau);
out:
    Py_XDECREF(py_pi);
    Py_XDECREF(py_tau);
    return ret;
}

#endif


static PyMethodDef TreeFile_methods[] = {
    {"get_num_loci", (PyCFunction) TreeFile_get_num_loci, METH_NOARGS,
            "Returns the number of loci"},
    {"get_sample_size", (PyCFunction) TreeFile_get_sample_size, METH_NOARGS,
            "Returns the sample size"},
    {"get_num_trees", (PyCFunction) TreeFile_get_num_trees, METH_NOARGS,
            "Returns the number of trees"},
    {"get_metadata", (PyCFunction) TreeFile_get_metadata, METH_NOARGS,
            "Returns the simulation metadata"},
    {"sort", (PyCFunction) TreeFile_sort, METH_NOARGS,
            "Sorts the coalescence records in the file."},
    /* {"get_tree", (PyCFunction) TreeFile_get_tree, METH_VARARGS, */
    /*         "Returns the tree corresponding to the jth breakpoint."}, */
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
    NULL, NULL, NULL, NULL, NULL
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
    PyObject *module = Py_InitModule3("_msprime", NULL, MODULE_DOC);
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

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}


