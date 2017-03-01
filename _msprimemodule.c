/*
** Copyright (C) 2014-2017 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <structmember.h>

#include <float.h>

#include <hdf5.h>
#include <gsl/gsl_version.h>
#include <gsl/gsl_math.h>

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
    unsigned long seed;
    gsl_rng* rng;
} RandomGenerator;

typedef struct {
    PyObject_HEAD
    mutgen_t *mutgen;
    RandomGenerator *random_generator;
} MutationGenerator;

typedef struct {
    PyObject_HEAD
    node_table_t *node_table;
} NodeTable;

typedef struct {
    PyObject_HEAD
    edgeset_table_t *edgeset_table;
} EdgesetTable;

typedef struct {
    PyObject_HEAD
    mutation_type_table_t *mutation_type_table;
} MutationTypeTable;

typedef struct {
    PyObject_HEAD
    mutation_table_t *mutation_table;
} MutationTable;

typedef struct {
    PyObject_HEAD
    migration_table_t *migration_table;
} MigrationTable;

typedef struct {
    PyObject_HEAD
    msp_t *sim;
    RandomGenerator *random_generator;
} Simulator;

typedef struct {
    PyObject_HEAD
    recomb_map_t *recomb_map;
} RecombinationMap;

typedef struct {
    PyObject_HEAD
    tree_sequence_t *tree_sequence;
} TreeSequence;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    sparse_tree_t *sparse_tree;
} SparseTree;

typedef struct {
    PyObject_HEAD
    SparseTree *sparse_tree;
    int first;
} SparseTreeIterator;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    tree_diff_iterator_t *tree_diff_iterator;
} TreeDiffIterator;

typedef struct {
    PyObject_HEAD
    SparseTree *sparse_tree;
    leaf_list_node_t *head;
    leaf_list_node_t *tail;
    leaf_list_node_t *next;
} LeafListIterator;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    newick_converter_t *newick_converter;
} NewickConverter;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    vcf_converter_t *vcf_converter;
} VcfConverter;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    hapgen_t *haplotype_generator;
} HaplotypeGenerator;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    PyObject *genotypes_buffer;
    vargen_t *variant_generator;
    Py_buffer buffer;
    int buffer_acquired;
} VariantGenerator;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    ld_calc_t *ld_calc;
} LdCalculator;

static void
handle_library_error(int err)
{
    if (err == MSP_ERR_OUT_OF_BOUNDS) {
        PyErr_SetString(PyExc_IndexError, msp_strerror(err));
    } else{
        PyErr_SetString(MsprimeLibraryError, msp_strerror(err));
    }
}

static void
handle_input_error(int err)
{
    PyErr_SetString(MsprimeInputError, msp_strerror(err));
}

static int
parse_sample_ids(PyObject *py_samples, tree_sequence_t *ts, size_t *num_samples,
        node_id_t **samples)
{
    int ret = -1;
    PyObject *item;
    size_t j;
    Py_ssize_t num_samples_local;
    node_id_t *samples_local = NULL;
    uint32_t n = tree_sequence_get_sample_size(ts);

    num_samples_local = PyList_Size(py_samples);
    if (num_samples_local < 2) {
        PyErr_SetString(PyExc_ValueError, "Must provide at least 2 samples");
        goto out;
    }
    samples_local = PyMem_Malloc(num_samples_local * sizeof(node_id_t));
    if (samples_local == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    for (j = 0; j < num_samples_local; j++) {
        item = PyList_GetItem(py_samples, j);
        if (!PyNumber_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "sample id must be a number");
            goto out;
        }
        samples_local[j] = (node_id_t) PyLong_AsLong(item);
        if (samples_local[j] < 0) {
            PyErr_SetString(PyExc_ValueError, "sample IDs must be >= 0");
            goto out;
        }
        if (samples_local[j] >= n) {
            PyErr_SetString(PyExc_ValueError, "sample ids must be < sample_size");
            goto out;
        }
    }
    *num_samples = (size_t) num_samples_local;
    *samples = samples_local;
    samples_local = NULL;
    ret = 0;
out:
    if (samples_local != NULL) {
        PyMem_Free(samples_local);
    }
    return ret;
}

static int
parse_samples(PyObject *py_samples, Py_ssize_t *sample_size, sample_t **samples)
{
    int ret = -1;
    long tmp_long;
    Py_ssize_t j, n;
    PyObject *sample, *value;
    sample_t *ret_samples = NULL;

    n = PyList_Size(py_samples);
    ret_samples = PyMem_Malloc(n * sizeof(sample_t));
    if (ret_samples == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    for (j = 0; j < n; j++) {
        sample = PyList_GetItem(py_samples, j);
        if (!PyTuple_Check(sample)) {
            PyErr_SetString(PyExc_TypeError, "not a tuple");
            goto out;
        }
        if (PyTuple_Size(sample) != 2) {
            PyErr_SetString(PyExc_ValueError,
                    "sample must be (population,time) tuple");
            goto out;
        }
        value = PyTuple_GetItem(sample, 0);
        if (!PyNumber_Check(value)) {
            PyErr_Format(PyExc_TypeError, "'population' is not number");
            goto out;
        }
        tmp_long = PyLong_AsLong(value);
        if (tmp_long < 0) {
            PyErr_SetString(PyExc_ValueError, "negative population IDs not valid");
            goto out;
        }
        ret_samples[j].population_id = (population_id_t) tmp_long;
        value = PyTuple_GetItem(sample, 1);
        if (!PyNumber_Check(value)) {
            PyErr_Format(PyExc_TypeError, "'time' is not number");
            goto out;
        }
        ret_samples[j].time = PyFloat_AsDouble(value);
        if (ret_samples[j].time < 0) {
            PyErr_SetString(PyExc_ValueError, "negative times not valid");
            goto out;
        }
    }
    *samples = ret_samples;
    *sample_size = n;
    ret = 0;
    ret_samples = NULL;
out:
    if (ret_samples != NULL) {
        PyMem_Free(ret_samples);
    }
    return ret;
}

static int
parse_provenance_strings(PyObject *py_provenance_strings, Py_ssize_t *num_provenance_strings,
        char ***provenance_strings)
{

    int ret = -1;
    Py_ssize_t j, n;
    PyObject *item;
    char *s;
    char **ret_provenance_strings = NULL;

    n = PyList_Size(py_provenance_strings);
    ret_provenance_strings = PyMem_Malloc(n * sizeof(char *));
    if (ret_provenance_strings == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    for (j = 0; j < n; j++) {
        item = PyList_GetItem(py_provenance_strings, j);
        assert(item != NULL);
        s = PyBytes_AsString(item);
        if (s == NULL) {
            goto out;
        }
        ret_provenance_strings[j] = s;
    }
    *num_provenance_strings = n;
    *provenance_strings = ret_provenance_strings;
    ret = 0;
out:
    return ret;
}

/*
 * Retrieves the PyObject* corresponding the specified key in the
 * specified dictionary.
 */
static PyObject *
get_dict_value(PyObject *dict, const char *key_str)
{
    PyObject *ret = NULL;
    PyObject *key = Py_BuildValue("s", key_str);

    if (key == NULL) {
        goto out;
    }
    if (!PyDict_Contains(dict, key)) {
        PyErr_Format(PyExc_ValueError, "'%s' not specified", key_str);
        goto out;
    }
    ret = PyDict_GetItem(dict, key);
    assert(ret != NULL);
out:
    Py_DECREF(key);
    return ret;
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

    value = get_dict_value(dict, key_str);
    if (value == NULL) {
        goto out;
    }
    if (!PyNumber_Check(value)) {
        PyErr_Format(PyExc_TypeError, "'%s' is not number", key_str);
        goto out;
    }
    ret = value;
out:
    return ret;
}


static PyObject *
convert_integer_list(size_t *list, size_t size)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_int = NULL;
    size_t j;

    l = PyList_New(size);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < size; j++) {
        py_int = Py_BuildValue("n", (Py_ssize_t) list[j]);
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

static PyObject *
convert_string_list(char **list, size_t size)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_str = NULL;
    size_t j;

    l = PyList_New(size);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < size; j++) {
        assert(list[j] != NULL);
        py_str = PyBytes_FromString(list[j]);
        if (py_str == NULL) {
            Py_DECREF(l);
            goto out;
        }
        PyList_SET_ITEM(l, j, py_str);
    }
    ret = l;
out:
    return ret;
}

static PyObject *
convert_node_id_list(node_id_t *children, size_t num_children)
{
    PyObject *ret = NULL;
    PyObject *t;
    PyObject *py_int;
    size_t j;

    t = PyTuple_New(num_children);
    if (t == NULL) {
        goto out;
    }
    for (j = 0; j < num_children; j++) {
        py_int = Py_BuildValue("i", (int) children[j]);
        if (py_int == NULL) {
            Py_DECREF(children);
            goto out;
        }
        PyTuple_SET_ITEM(t, j, py_int);
    }
    ret = t;
out:
    return ret;
}

static PyObject *
convert_float_list(double *list, size_t size)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_float = NULL;
    size_t j;

    l = PyList_New(size);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < size; j++) {
        py_float = Py_BuildValue("d", list[j]);
        if (py_float == NULL) {
            Py_DECREF(l);
            goto out;
        }
        PyList_SET_ITEM(l, j, py_float);
    }
    ret = l;
out:
    return ret;
}

static PyObject *
make_coalescence_record(coalescence_record_t *cr)
{
    PyObject *children = NULL;
    PyObject *ret = NULL;

    children = convert_node_id_list(cr->children, cr->num_children);
    if (children == NULL) {
        goto out;
    }
    ret = Py_BuildValue("ddiOdi",
            cr->left, cr->right, (int) cr->node, children, cr->time,
            (int) cr->population_id);
out:
    Py_XDECREF(children);
    return ret;
}

static PyObject *
make_coalescence_record_tmp(node_t *node, edgeset_t *edgeset)
{
    PyObject *children = NULL;
    PyObject *ret = NULL;

    children = convert_node_id_list(edgeset->children, edgeset->num_children);
    if (children == NULL) {
        goto out;
    }
    ret = Py_BuildValue("ddIOdi",
            edgeset->left, edgeset->right, (int) edgeset->parent,
            children, node->time, (int) node->population);
out:
    Py_XDECREF(children);
    return ret;
}

static PyObject *
make_node(node_t *r)
{
    PyObject *ret = NULL;

    ret = Py_BuildValue("Idis",
        (unsigned int) r->flags, r->time, (int) r->population, r->name);
    return ret;
}

static PyObject *
make_edgeset(edgeset_t *edgeset)
{
    PyObject *children = NULL;
    PyObject *ret = NULL;

    children = convert_node_id_list(edgeset->children, edgeset->num_children);
    if (children == NULL) {
        goto out;
    }
    ret = Py_BuildValue("ddiO",
            edgeset->left, edgeset->right, (int) edgeset->parent, children);
out:
    Py_XDECREF(children);
    return ret;
}

static PyObject *
make_migration(migration_t *r)
{
    int source = r->source == MSP_NULL_POPULATION_ID ? -1: r->source;
    int dest = r->dest == MSP_NULL_POPULATION_ID ? -1: r->dest;
    PyObject *ret = NULL;

    ret = Py_BuildValue("ddiiid",
            r->left, r->right, (int) r->node, source, dest, r->time);
    return ret;
}

static PyObject *
make_mutation_type(mutation_type_t *mutation_type)
{
    PyObject *ret = NULL;

    ret = Py_BuildValue("ss", mutation_type->ancestral_state,
            mutation_type->derived_state);
    return ret;
}

static PyObject *
make_mutation(mutation_t *mutation)
{
    PyObject *nodes = NULL;
    PyObject *ret = NULL;

    nodes = convert_node_id_list(mutation->nodes, mutation->num_nodes);
    if (nodes == NULL) {
        goto out;
    }
    ret = Py_BuildValue("dOin", mutation->position, nodes,
           (int) mutation->type, (Py_ssize_t) mutation->index);
out:
    Py_XDECREF(nodes);
    return ret;
}

static PyObject *
convert_mutations(mutation_t *mutations, size_t num_mutations)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_mutation = NULL;
    size_t j;

    l = PyList_New(num_mutations);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < num_mutations; j++) {
        py_mutation = make_mutation(&mutations[j]);
        if (py_mutation == NULL) {
            Py_DECREF(l);
            goto out;
        }
        PyList_SET_ITEM(l, j, py_mutation);
    }
    ret = l;
out:
    return ret;
}

static int
parse_node_tuple(PyObject *py_nodes, size_t *size, node_id_t **nodes)
{
    int ret = -1;
    PyObject *item;
    size_t j;
    Py_ssize_t num_nodes_local;
    node_id_t *nodes_local = NULL;

    num_nodes_local = PyTuple_Size(py_nodes);
    nodes_local = PyMem_Malloc(num_nodes_local * sizeof(node_id_t));
    if (nodes_local == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    for (j = 0; j < num_nodes_local; j++) {
        item = PyTuple_GetItem(py_nodes, j);
        if (!PyNumber_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "node id must be a number");
            goto out;
        }
        nodes_local[j] = (node_id_t) PyLong_AsLong(item);
    }
    *size = (size_t) num_nodes_local;
    *nodes = nodes_local;
    nodes_local = NULL;
    ret = 0;
out:
    if (nodes_local != NULL) {
        PyMem_Free(nodes_local);
    }
    return ret;
}

/*===================================================================
 * RandomGenerator
 *===================================================================
 */

static int
RandomGenerator_check_state(RandomGenerator *self)
{
    int ret = 0;
    if (self->rng == NULL) {
        PyErr_SetString(PyExc_SystemError, "RandomGenerator not initialised");
        ret = -1;
    }
    return ret;
}

static void
RandomGenerator_dealloc(RandomGenerator* self)
{
    if (self->rng != NULL) {
        gsl_rng_free(self->rng);
        self->rng = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
RandomGenerator_init(RandomGenerator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    static char *kwlist[] = {"seed", NULL};
    unsigned long long seed = 0;

    self->rng  = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "K", kwlist, &seed)) {
        goto out;
    }
    if (seed == 0 || seed >= (1ULL<<32)) {
        PyErr_Format(PyExc_ValueError,
            "seeds must be greater than 0 and less than 2^32");
        goto out;
    }
    self->seed = seed;
    self->rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(self->rng, self->seed);
    ret = 0;
out:
    return ret;
}

static PyObject *
RandomGenerator_get_seed(RandomGenerator *self)
{
    PyObject *ret = NULL;

    if (RandomGenerator_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("k", self->seed);
out:
    return ret;
}

static PyMemberDef RandomGenerator_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef RandomGenerator_methods[] = {
    {"get_seed", (PyCFunction) RandomGenerator_get_seed,
        METH_NOARGS, "Returns the random seed for this generator."},
    {NULL}  /* Sentinel */
};

static PyTypeObject RandomGeneratorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.RandomGenerator",             /* tp_name */
    sizeof(RandomGenerator),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)RandomGenerator_dealloc, /* tp_dealloc */
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
    "RandomGenerator objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    RandomGenerator_methods,             /* tp_methods */
    RandomGenerator_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)RandomGenerator_init,      /* tp_init */
};

/*===================================================================
 * General table code.
 *===================================================================
 */

static PyObject *
table_get_column_array(size_t num_rows, void *data, int npy_type,
        size_t element_size)
{
    PyObject *ret = NULL;
    PyArrayObject *array;
    npy_intp dims = (npy_intp) num_rows;

    array = (PyArrayObject *) PyArray_EMPTY(1, &dims, npy_type, 0);
    if (array == NULL) {
        goto out;
    }
    memcpy(PyArray_DATA(array), data, num_rows * element_size);
    ret = (PyObject *) array;
out:
    return ret;
}

static PyArrayObject *
table_read_column_array(PyObject *input, int npy_type, size_t *num_rows, bool check_num_rows)
{
    PyArrayObject *ret = NULL;
    PyArrayObject *array = NULL;
    npy_intp *shape;

    array = (PyArrayObject *) PyArray_FROM_OTF(input, npy_type, NPY_ARRAY_IN_ARRAY);
    if (array == NULL) {
        goto out;
    }
    if (PyArray_NDIM(array) != 1) {
        PyErr_SetString(PyExc_ValueError, "Dim != 1");
        goto out;
    }
    shape = PyArray_DIMS(array);
    if (check_num_rows) {
        if (*num_rows != shape[0]) {
            PyErr_SetString(PyExc_ValueError, "Input array dimensions must be equal.");
            goto out;
        }
    } else {
        *num_rows = shape[0];
    }
    ret = array;
out:
    if (ret == NULL) {
        Py_XDECREF(array);
    }
    return ret;
}

/*===================================================================
 * NodeTable
 *===================================================================
 */

static int
NodeTable_check_state(NodeTable *self)
{
    int ret = 0;
    if (self->node_table == NULL) {
        PyErr_SetString(PyExc_SystemError, "NodeTable not initialised");
        ret = -1;
        goto out;
    }
out:
    return ret;
}

static void
NodeTable_dealloc(NodeTable* self)
{
    if (self->node_table != NULL) {
        node_table_free(self->node_table);
        PyMem_Free(self->node_table);
        self->node_table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
NodeTable_init(NodeTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"max_rows_increment", NULL};
    Py_ssize_t max_rows_increment = 1024;
    Py_ssize_t max_name_length_increment = 1;

    self->node_table = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist,
                &max_rows_increment)) {
        goto out;
    }
    if (max_rows_increment <= 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    self->node_table = PyMem_Malloc(sizeof(node_table_t));
    if (self->node_table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = node_table_alloc(self->node_table, (size_t) max_rows_increment,
            (size_t) max_name_length_increment);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
NodeTable_add_row(NodeTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    unsigned int flags = 0;
    double time = 0;
    int population = -1;
    char *name = "";
    Py_ssize_t name_length = 0;
    static char *kwlist[] = {"flags", "time", "population", "name", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|idis#", kwlist,
                &flags, &time, &population, &name, &name_length)) {
        goto out;
    }
    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    err = node_table_add_row(self->node_table, (uint32_t) flags, time,
            (population_id_t) population, name);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
NodeTable_set_columns(NodeTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    size_t num_rows, name_length;
    char *name_data;
    void *population_data;
    PyObject *time_input = NULL;
    PyArrayObject *time_array = NULL;
    PyObject *flags_input = NULL;
    PyArrayObject *flags_array = NULL;
    PyObject *population_input = NULL;
    PyArrayObject *population_array = NULL;
    PyObject *name_input = NULL;
    PyArrayObject *name_array = NULL;
    static char *kwlist[] = {"flags", "time", "population", "name", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|OO", kwlist,
                &flags_input, &time_input, &population_input, &name_input)) {
        goto out;
    }
    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    flags_array = table_read_column_array(flags_input, NPY_UINT32, &num_rows, false);
    if (flags_array == NULL) {
        goto out;
    }
    time_array = table_read_column_array(time_input, NPY_FLOAT64, &num_rows, true);
    if (time_array == NULL) {
        goto out;
    }
    population_data = NULL;
    if (population_input != NULL) {
        population_array = table_read_column_array(population_input, NPY_INT32,
                &num_rows, true);
        if (population_array == NULL) {
            goto out;
        }
        population_data = PyArray_DATA(population_array);
    }
    name_length = 0;
    name_data = NULL;
    if (name_input != NULL) {
        name_array = table_read_column_array(name_input, NPY_INT8, &name_length, false);
        if (name_array == NULL) {
            goto out;
        }
        name_data = PyArray_DATA(name_array);
    }
    err = node_table_set_columns(self->node_table, num_rows,
            PyArray_DATA(flags_array), PyArray_DATA(time_array), population_data,
            name_length, name_data);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    Py_XDECREF(flags_array);
    Py_XDECREF(time_array);
    return ret;
}

static PyObject *
NodeTable_get_max_rows_increment(NodeTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->node_table->max_rows_increment);
out:
    return ret;
}

static PyObject *
NodeTable_get_num_rows(NodeTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->node_table->num_rows);
out:
    return ret;
}

static PyObject *
NodeTable_get_time(NodeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->node_table->num_rows, self->node_table->time,
            NPY_FLOAT64, sizeof(double));
out:
    return ret;
}

static PyObject *
NodeTable_get_flags(NodeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->node_table->num_rows, self->node_table->flags,
            NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyObject *
NodeTable_get_population(NodeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->node_table->num_rows, self->node_table->population,
            NPY_INT32, sizeof(int32_t));
out:
    return ret;
}

static PyObject *
NodeTable_get_name(NodeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->node_table->name_length,
            self->node_table->name, NPY_INT8, sizeof(char));
out:
    return ret;
}



static PyGetSetDef NodeTable_getsetters[] = {
    {"max_rows_increment",
        (getter) NodeTable_get_max_rows_increment, NULL, "The size increment"},
    {"num_rows", (getter) NodeTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"time", (getter) NodeTable_get_time, NULL, "The time array"},
    {"flags", (getter) NodeTable_get_flags, NULL, "The flags array"},
    {"population", (getter) NodeTable_get_population, NULL, "The population array"},
    {"name", (getter) NodeTable_get_name, NULL, "The name array"},
    {NULL}  /* Sentinel */
};

static PyMethodDef NodeTable_methods[] = {
    {"add_row", (PyCFunction) NodeTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    {"set_columns", (PyCFunction) NodeTable_set_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the speficied arrays into the columns."},
    {NULL}  /* Sentinel */
};

static PyTypeObject NodeTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.NodeTable",             /* tp_name */
    sizeof(NodeTable),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)NodeTable_dealloc, /* tp_dealloc */
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
    "NodeTable objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    NodeTable_methods,             /* tp_methods */
    0,                             /* tp_members */
    NodeTable_getsetters,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)NodeTable_init,      /* tp_init */
};

/*===================================================================
 * EdgesetTable
 *===================================================================
 */

static int
EdgesetTable_check_state(EdgesetTable *self)
{
    int ret = 0;
    if (self->edgeset_table == NULL) {
        PyErr_SetString(PyExc_SystemError, "EdgesetTable not initialised");
        ret = -1;
        goto out;
    }
out:
    return ret;
}

static void
EdgesetTable_dealloc(EdgesetTable* self)
{
    if (self->edgeset_table != NULL) {
        edgeset_table_free(self->edgeset_table);
        PyMem_Free(self->edgeset_table);
        self->edgeset_table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
EdgesetTable_init(EdgesetTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {
        "max_rows_increment", "max_children_length_increment", NULL};
    Py_ssize_t max_rows_increment = 1024;
    Py_ssize_t max_children_length_increment = 1024;

    self->edgeset_table = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|nn", kwlist,
                &max_rows_increment, &max_children_length_increment)) {
        goto out;
    }
    if (max_rows_increment <= 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    if (max_children_length_increment <= 0) {
        PyErr_SetString(PyExc_ValueError, "max_children_length_increment must be positive");
        goto out;
    }
    self->edgeset_table = PyMem_Malloc(sizeof(edgeset_table_t));
    if (self->edgeset_table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = edgeset_table_alloc(self->edgeset_table, max_rows_increment,
            max_children_length_increment);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
EdgesetTable_add_row(EdgesetTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    double left = 0.0;
    double right = 1.0;
    int parent;
    PyObject * py_children;
    node_id_t *children = NULL;
    size_t num_children;
    static char *kwlist[] = {"left", "right", "parent", "children", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddiO!", kwlist,
                &left, &right, &parent, &PyTuple_Type, &py_children)) {
        goto out;
    }
    if (EdgesetTable_check_state(self) != 0) {
        goto out;
    }
    err = parse_node_tuple(py_children, &num_children, &children);
    if (err != 0) {
        goto out;
    }
    err = edgeset_table_add_row(self->edgeset_table, left, right, parent,
        num_children, children);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    if (children != NULL) {
        PyMem_Free(children);
    }
    return ret;
}


static PyObject *
EdgesetTable_set_columns(EdgesetTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    size_t num_rows = 0;
    size_t children_length = 0;
    PyObject *left_input = NULL;
    PyArrayObject *left_array = NULL;
    PyObject *right_input = NULL;
    PyArrayObject *right_array = NULL;
    PyObject *parent_input = NULL;
    PyArrayObject *parent_array = NULL;
    PyObject *children_input = NULL;
    PyArrayObject *children_array = NULL;
    static char *kwlist[] = {"left", "right", "parent", "children", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOO", kwlist,
                &left_input, &right_input, &parent_input, &children_input)) {
        goto out;
    }

    left_array = table_read_column_array(left_input, NPY_FLOAT64, &num_rows, false);
    if (left_array == NULL) {
        goto out;
    }
    right_array = table_read_column_array(right_input, NPY_FLOAT64, &num_rows, true);
    if (right_array == NULL) {
        goto out;
    }
    parent_array = table_read_column_array(parent_input, NPY_INT32, &num_rows, true);
    if (parent_array == NULL) {
        goto out;
    }
    children_array = table_read_column_array(children_input, NPY_INT32,
            &children_length, false);
    if (children_array == NULL) {
        goto out;
    }
    err = edgeset_table_set_columns(self->edgeset_table, num_rows,
            PyArray_DATA(left_array), PyArray_DATA(right_array),
            PyArray_DATA(parent_array),
            children_length, PyArray_DATA(children_array));
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    Py_XDECREF(left_array);
    Py_XDECREF(right_array);
    Py_XDECREF(parent_array);
    Py_XDECREF(children_array);
    return ret;
}

static PyObject *
EdgesetTable_get_max_rows_increment(EdgesetTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (EdgesetTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->edgeset_table->max_rows_increment);
out:
    return ret;
}

static PyObject *
EdgesetTable_get_max_children_length_increment(EdgesetTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (EdgesetTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->edgeset_table->max_children_length_increment);
out:
    return ret;
}

static PyObject *
EdgesetTable_get_num_rows(EdgesetTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (EdgesetTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->edgeset_table->num_rows);
out:
    return ret;
}

static PyObject *
EdgesetTable_get_left(EdgesetTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (EdgesetTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->edgeset_table->num_rows, self->edgeset_table->left, NPY_FLOAT64,
            sizeof(double));
out:
    return ret;
}

static PyObject *
EdgesetTable_get_right(EdgesetTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (EdgesetTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->edgeset_table->num_rows, self->edgeset_table->right, NPY_FLOAT64,
            sizeof(double));
out:
    return ret;
}

static PyObject *
EdgesetTable_get_parent(EdgesetTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (EdgesetTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->edgeset_table->num_rows, self->edgeset_table->parent, NPY_INT32,
            sizeof(int32_t));
out:
    return ret;
}

static PyObject *
EdgesetTable_get_children(EdgesetTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (EdgesetTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->edgeset_table->children_length, self->edgeset_table->children, NPY_INT32,
            sizeof(int32_t));
out:
    return ret;
}

static PyGetSetDef EdgesetTable_getsetters[] = {
    {"max_rows_increment",
        (getter) EdgesetTable_get_max_rows_increment, NULL,
        "The size increment"},
    {"max_children_length_increment",
        (getter) EdgesetTable_get_max_children_length_increment, NULL,
        "The total children increment"},
    {"num_rows", (getter) EdgesetTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"left", (getter) EdgesetTable_get_left, NULL, "The left array"},
    {"right", (getter) EdgesetTable_get_right, NULL, "The right array"},
    {"parent", (getter) EdgesetTable_get_parent, NULL, "The parent array"},
    {"children", (getter) EdgesetTable_get_children, NULL, "The children array"},
    {NULL}  /* Sentinel */
};

static PyMethodDef EdgesetTable_methods[] = {
    {"add_row", (PyCFunction) EdgesetTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    {"set_columns", (PyCFunction) EdgesetTable_set_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the speficied arrays into the columns."},
    {NULL}  /* Sentinel */
};

static PyTypeObject EdgesetTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.EdgesetTable",             /* tp_name */
    sizeof(EdgesetTable),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)EdgesetTable_dealloc, /* tp_dealloc */
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
    "EdgesetTable objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    EdgesetTable_methods,             /* tp_methods */
    0,                             /* tp_members */
    EdgesetTable_getsetters,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)EdgesetTable_init,      /* tp_init */
};

/*===================================================================
 * MigrationTable
 *===================================================================
 */

static int
MigrationTable_check_state(MigrationTable *self)
{
    int ret = 0;
    if (self->migration_table == NULL) {
        PyErr_SetString(PyExc_SystemError, "MigrationTable not initialised");
        ret = -1;
        goto out;
    }
out:
    return ret;
}

static void
MigrationTable_dealloc(MigrationTable* self)
{
    if (self->migration_table != NULL) {
        migration_table_free(self->migration_table);
        PyMem_Free(self->migration_table);
        self->migration_table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
MigrationTable_init(MigrationTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"max_rows_increment", NULL};
    Py_ssize_t max_rows_increment = 1024;

    self->migration_table = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist,
                &max_rows_increment)) {
        goto out;
    }
    if (max_rows_increment <= 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    self->migration_table = PyMem_Malloc(sizeof(migration_table_t));
    if (self->migration_table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = migration_table_alloc(self->migration_table, (size_t) max_rows_increment);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
MigrationTable_set_columns(MigrationTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    size_t num_rows;
    PyObject *left_input = NULL;
    PyArrayObject *left_array = NULL;
    PyObject *right_input = NULL;
    PyArrayObject *right_array = NULL;
    PyObject *node_input = NULL;
    PyArrayObject *node_array = NULL;
    PyObject *source_input = NULL;
    PyArrayObject *source_array = NULL;
    PyObject *dest_input = NULL;
    PyArrayObject *dest_array = NULL;
    PyObject *time_input = NULL;
    PyArrayObject *time_array = NULL;
    static char *kwlist[] = {"left", "right", "node", "source", "dest", "time", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOOOO", kwlist,
                &left_input, &right_input, &node_input, &source_input, &dest_input,
                &time_input)) {
        goto out;
    }
    left_array = table_read_column_array(left_input, NPY_FLOAT64, &num_rows, false);
    if (left_array == NULL) {
        goto out;
    }
    right_array = table_read_column_array(right_input, NPY_FLOAT64, &num_rows, true);
    if (right_array == NULL) {
        goto out;
    }
    node_array = table_read_column_array(node_input, NPY_INT32, &num_rows, true);
    if (node_array == NULL) {
        goto out;
    }
    source_array = table_read_column_array(source_input, NPY_INT32, &num_rows, true);
    if (source_array == NULL) {
        goto out;
    }
    dest_array = table_read_column_array(dest_input, NPY_INT32, &num_rows, true);
    if (dest_array == NULL) {
        goto out;
    }
    time_array = table_read_column_array(time_input, NPY_FLOAT64, &num_rows, true);
    if (time_array == NULL) {
        goto out;
    }
    err = migration_table_set_columns(self->migration_table, num_rows,
            PyArray_DATA(left_array), PyArray_DATA(right_array), PyArray_DATA(node_array),
            PyArray_DATA(source_array), PyArray_DATA(dest_array), PyArray_DATA(time_array));
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    Py_XDECREF(left_array);
    Py_XDECREF(right_array);
    Py_XDECREF(node_array);
    Py_XDECREF(source_array);
    Py_XDECREF(dest_array);
    Py_XDECREF(time_array);
    return ret;
}

static PyObject *
MigrationTable_get_max_rows_increment(MigrationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->migration_table->max_rows_increment);
out:
    return ret;
}

static PyObject *
MigrationTable_get_num_rows(MigrationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->migration_table->num_rows);
out:
    return ret;
}

static PyObject *
MigrationTable_get_left(MigrationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->migration_table->num_rows, self->migration_table->left,
            NPY_FLOAT64, sizeof(double));
out:
    return ret;
}

static PyObject *
MigrationTable_get_right(MigrationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->migration_table->num_rows, self->migration_table->right,
            NPY_FLOAT64, sizeof(double));
out:
    return ret;
}

static PyObject *
MigrationTable_get_time(MigrationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->migration_table->num_rows, self->migration_table->time,
            NPY_FLOAT64, sizeof(double));
out:
    return ret;
}

static PyObject *
MigrationTable_get_node(MigrationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->migration_table->num_rows, self->migration_table->node,
            NPY_INT32, sizeof(uint32_t));
out:
    return ret;
}

static PyObject *
MigrationTable_get_source(MigrationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->migration_table->num_rows, self->migration_table->source,
            NPY_INT32, sizeof(uint32_t));
out:
    return ret;
}

static PyObject *
MigrationTable_get_dest(MigrationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->migration_table->num_rows, self->migration_table->dest,
            NPY_INT32, sizeof(uint32_t));
out:
    return ret;
}

static PyGetSetDef MigrationTable_getsetters[] = {
    {"max_rows_increment",
        (getter) MigrationTable_get_max_rows_increment, NULL, "The size increment"},
    {"num_rows", (getter) MigrationTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"left", (getter) MigrationTable_get_left, NULL, "The left array"},
    {"right", (getter) MigrationTable_get_right, NULL, "The right array"},
    {"node", (getter) MigrationTable_get_node, NULL, "The node array"},
    {"source", (getter) MigrationTable_get_source, NULL, "The source array"},
    {"dest", (getter) MigrationTable_get_dest, NULL, "The dest array"},
    {"time", (getter) MigrationTable_get_time, NULL, "The time array"},
    {NULL}  /* Sentinel */
};

static PyMethodDef MigrationTable_methods[] = {
    {"set_columns", (PyCFunction) MigrationTable_set_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the speficied arrays into the columns."},
    {NULL}  /* Sentinel */
};

static PyTypeObject MigrationTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.MigrationTable",             /* tp_name */
    sizeof(MigrationTable),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)MigrationTable_dealloc, /* tp_dealloc */
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
    "MigrationTable objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    MigrationTable_methods,             /* tp_methods */
    0,                             /* tp_members */
    MigrationTable_getsetters,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)MigrationTable_init,      /* tp_init */
};


/*===================================================================
 * MutationTypeTable
 *===================================================================
 */

static int
MutationTypeTable_check_state(MutationTypeTable *self)
{
    int ret = 0;
    if (self->mutation_type_table == NULL) {
        PyErr_SetString(PyExc_SystemError, "MutationTypeTable not initialised");
        ret = -1;
        goto out;
    }
out:
    return ret;
}

static void
MutationTypeTable_dealloc(MutationTypeTable* self)
{
    if (self->mutation_type_table != NULL) {
        mutation_type_table_free(self->mutation_type_table);
        PyMem_Free(self->mutation_type_table);
        self->mutation_type_table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
MutationTypeTable_init(MutationTypeTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"max_rows_increment", NULL};
    Py_ssize_t max_rows_increment = 1;

    self->mutation_type_table = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist, &max_rows_increment)) {
        goto out;
    }
    if (max_rows_increment <= 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    self->mutation_type_table = PyMem_Malloc(sizeof(mutation_type_table_t));
    if (self->mutation_type_table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = mutation_type_table_alloc(self->mutation_type_table, max_rows_increment);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
MutationTypeTable_add_row(MutationTypeTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    char *ancestral_state = NULL;
    char *derived_state = NULL;
    static char *kwlist[] = {"ancestral_state", "derived_state", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ss", kwlist,
                &ancestral_state, &derived_state)) {
        goto out;
    }
    if (MutationTypeTable_check_state(self) != 0) {
        goto out;
    }
    err = mutation_type_table_add_row(self->mutation_type_table, ancestral_state,
            derived_state);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

/* TODO implement */
/* static PyObject * */
/* MutationTypeTable_set_columns(MutationTypeTable *self, PyObject *args, PyObject *kwds) */
/* { */
/*     PyObject *ret = NULL; */
/*     int err; */
/*     size_t num_rows = 0; */
/*     size_t nodes_length = 0; */
/*     PyObject *position_input = NULL; */
/*     PyArrayObject *position_array = NULL; */
/*     PyObject *nodes_input = NULL; */
/*     PyArrayObject *nodes_array = NULL; */

/*     static char *kwlist[] = {"position", "nodes", NULL}; */

/*     if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO", kwlist, */
/*                 &position_input, &nodes_input)) { */
/*         goto out; */
/*     } */
/*     position_array = table_read_column_array(position_input, NPY_FLOAT64, */
/*             &num_rows, false); */
/*     if (position_array == NULL) { */
/*         goto out; */
/*     } */
/*     nodes_array = table_read_column_array(nodes_input, NPY_INT32, &nodes_length, false); */
/*     if (nodes_array == NULL) { */
/*         goto out; */
/*     } */
/*     err = mutation_type_table_set_columns(self->mutation_type_table, num_rows, */
/*             PyArray_DATA(position_array), nodes_length, PyArray_DATA(nodes_array)); */
/*     if (err != 0) { */
/*         handle_library_error(err); */
/*         goto out; */
/*     } */
/*     ret = Py_BuildValue(""); */
/* out: */
/*     Py_XDECREF(position_array); */
/*     Py_XDECREF(nodes_array); */
/*     return ret; */
/* } */

static PyObject *
MutationTypeTable_get_max_rows_increment(MutationTypeTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (MutationTypeTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->mutation_type_table->max_rows_increment);
out:
    return ret;
}

static PyObject *
MutationTypeTable_get_num_rows(MutationTypeTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (MutationTypeTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->mutation_type_table->num_rows);
out:
    return ret;
}

static PyObject *
MutationTypeTable_get_ancestral_state(MutationTypeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTypeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->mutation_type_table->num_rows, self->mutation_type_table->ancestral_state,
            NPY_INT8, sizeof(char));
out:
    return ret;
}

static PyObject *
MutationTypeTable_get_derived_state(MutationTypeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTypeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->mutation_type_table->num_rows, self->mutation_type_table->derived_state,
            NPY_INT8, sizeof(char));
out:
    return ret;
}

static PyGetSetDef MutationTypeTable_getsetters[] = {
    {"max_rows_increment",
        (getter) MutationTypeTable_get_max_rows_increment, NULL,
        "The size increment"},
    {"num_rows",
        (getter) MutationTypeTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"ancestral_state", (getter) MutationTypeTable_get_ancestral_state, NULL,
        "The ancestral state array."},
    {"derived_state", (getter) MutationTypeTable_get_derived_state, NULL,
        "The derived state array."},
    {NULL}  /* Sentinel */
};

static PyMethodDef MutationTypeTable_methods[] = {
    {"add_row", (PyCFunction) MutationTypeTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    /* {"set_columns", (PyCFunction) MutationTypeTable_set_columns, METH_VARARGS|METH_KEYWORDS, */
    /*     "Copies the data in the speficied arrays into the columns."}, */
    {NULL}  /* Sentinel */
};

static PyTypeObject MutationTypeTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.MutationTypeTable",             /* tp_name */
    sizeof(MutationTypeTable),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)MutationTypeTable_dealloc, /* tp_dealloc */
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
    "MutationTypeTable objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    MutationTypeTable_methods,             /* tp_methods */
    0,                             /* tp_members */
    MutationTypeTable_getsetters,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)MutationTypeTable_init,      /* tp_init */
};


/*===================================================================
 * MutationTable
 *===================================================================
 */

static int
MutationTable_check_state(MutationTable *self)
{
    int ret = 0;
    if (self->mutation_table == NULL) {
        PyErr_SetString(PyExc_SystemError, "MutationTable not initialised");
        ret = -1;
        goto out;
    }
out:
    return ret;
}

static void
MutationTable_dealloc(MutationTable* self)
{
    if (self->mutation_table != NULL) {
        mutation_table_free(self->mutation_table);
        PyMem_Free(self->mutation_table);
        self->mutation_table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
MutationTable_init(MutationTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {
        "max_rows_increment", "max_nodes_length_increment", NULL};
    Py_ssize_t max_rows_increment = 1024;
    Py_ssize_t max_nodes_length_increment = 1024;

    self->mutation_table = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|nn", kwlist,
                &max_rows_increment, &max_nodes_length_increment)) {
        goto out;
    }
    if (max_rows_increment <= 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    if (max_nodes_length_increment <= 0) {
        PyErr_SetString(PyExc_ValueError, "max_nodes_length_increment must be positive");
        goto out;
    }
    self->mutation_table = PyMem_Malloc(sizeof(mutation_table_t));
    if (self->mutation_table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = mutation_table_alloc(self->mutation_table, max_rows_increment,
            max_nodes_length_increment);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
MutationTable_add_row(MutationTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    double position = 0.0;
    PyObject * py_nodes;
    node_id_t *nodes = NULL;
    int type = 0;
    size_t num_nodes;
    static char *kwlist[] = {"position", "nodes", "type", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dO!|i", kwlist,
                &position, &PyTuple_Type, &py_nodes, &type)) {
        goto out;
    }
    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    err = parse_node_tuple(py_nodes, &num_nodes, &nodes);
    if (err != 0) {
        goto out;
    }
    err = mutation_table_add_row(self->mutation_table, position,
            (mutation_type_id_t) type, num_nodes, nodes);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    if (nodes != NULL) {
        PyMem_Free(nodes);
    }
    return ret;
}


static PyObject *
MutationTable_set_columns(MutationTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    size_t num_rows = 0;
    size_t nodes_length = 0;
    PyObject *position_input = NULL;
    PyArrayObject *position_array = NULL;
    PyObject *nodes_input = NULL;
    PyArrayObject *nodes_array = NULL;
    PyObject *type_input = NULL;
    PyArrayObject *type_array = NULL;

    static char *kwlist[] = {"position", "nodes", "type", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOO", kwlist,
                &position_input, &nodes_input, &type_input)) {
        goto out;
    }
    position_array = table_read_column_array(position_input, NPY_FLOAT64,
            &num_rows, false);
    if (position_array == NULL) {
        goto out;
    }
    nodes_array = table_read_column_array(nodes_input, NPY_INT32, &nodes_length, false);
    if (nodes_array == NULL) {
        goto out;
    }
    type_array = table_read_column_array(type_input, NPY_UINT8, &num_rows, true);
    if (type_array == NULL) {
        goto out;
    }
    err = mutation_table_set_columns(self->mutation_table, num_rows,
            PyArray_DATA(position_array), PyArray_DATA(type_array),
            nodes_length, PyArray_DATA(nodes_array));
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    Py_XDECREF(position_array);
    Py_XDECREF(nodes_array);
    return ret;
}

static PyObject *
MutationTable_get_max_rows_increment(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->mutation_table->max_rows_increment);
out:
    return ret;
}

static PyObject *
MutationTable_get_max_nodes_length_increment(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->mutation_table->max_nodes_length_increment);
out:
    return ret;
}

static PyObject *
MutationTable_get_num_rows(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->mutation_table->num_rows);
out:
    return ret;
}

static PyObject *
MutationTable_get_position(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->mutation_table->num_rows, self->mutation_table->position, NPY_FLOAT64,
            sizeof(double));
out:
    return ret;
}

static PyObject *
MutationTable_get_type(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->mutation_table->num_rows, self->mutation_table->type, NPY_UINT8,
            sizeof(uint8_t));
out:
    return ret;
}

static PyObject *
MutationTable_get_nodes(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->mutation_table->nodes_length, self->mutation_table->nodes, NPY_INT32,
            sizeof(uint32_t));
out:
    return ret;
}

static PyGetSetDef MutationTable_getsetters[] = {
    {"max_rows_increment",
        (getter) MutationTable_get_max_rows_increment, NULL,
        "The size increment"},
    {"max_nodes_length_increment",
        (getter) MutationTable_get_max_nodes_length_increment, NULL,
        "The total nodes increment"},
    {"num_rows",
        (getter) MutationTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"position", (getter) MutationTable_get_position, NULL, "The position array"},
    {"type", (getter) MutationTable_get_type, NULL, "The type array"},
    {"nodes", (getter) MutationTable_get_nodes, NULL, "The nodes array"},
    {NULL}  /* Sentinel */
};

static PyMethodDef MutationTable_methods[] = {
    {"add_row", (PyCFunction) MutationTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    {"set_columns", (PyCFunction) MutationTable_set_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the speficied arrays into the columns."},
    {NULL}  /* Sentinel */
};

static PyTypeObject MutationTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.MutationTable",             /* tp_name */
    sizeof(MutationTable),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)MutationTable_dealloc, /* tp_dealloc */
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
    "MutationTable objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    MutationTable_methods,             /* tp_methods */
    0,                             /* tp_members */
    MutationTable_getsetters,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)MutationTable_init,      /* tp_init */
};

/*===================================================================
 * MutationGenerator
 *===================================================================
 */

static int
MutationGenerator_check_state(MutationGenerator *self)
{
    int ret = 0;
    if (self->mutgen == NULL) {
        PyErr_SetString(PyExc_SystemError, "MutationGenerator not initialised");
        ret = -1;
        goto out;
    }
    ret = RandomGenerator_check_state(self->random_generator);
out:
    return ret;
}

static void
MutationGenerator_dealloc(MutationGenerator* self)
{
    if (self->mutgen != NULL) {
        mutgen_free(self->mutgen);
        PyMem_Free(self->mutgen);
        self->mutgen = NULL;
    }
    Py_XDECREF(self->random_generator);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
MutationGenerator_init(MutationGenerator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"random_generator", "mutation_rate", NULL};
    double mutation_rate = 0;
    RandomGenerator *random_generator = NULL;

    self->mutgen = NULL;
    self->random_generator = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!d", kwlist,
            &RandomGeneratorType, &random_generator, &mutation_rate)) {
        goto out;
    }
    self->random_generator = random_generator;
    Py_INCREF(self->random_generator);
    if (RandomGenerator_check_state(self->random_generator) != 0) {
        goto out;
    }
    if (mutation_rate < 0) {
        PyErr_Format(PyExc_ValueError, "mutation_rate must be >= 0");
        goto out;
    }
    self->mutgen = PyMem_Malloc(sizeof(mutgen_t));
    if (self->mutgen == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = mutgen_alloc(self->mutgen, mutation_rate, random_generator->rng);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
MutationGenerator_get_mutation_rate(MutationGenerator *self)
{
    PyObject *ret = NULL;

    if (MutationGenerator_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->mutgen->mutation_rate);
out:
    return ret;
}


static PyObject *
MutationGenerator_generate(MutationGenerator *self, PyObject *args, PyObject *kwds)
{
    int err;
    PyObject *ret = NULL;
    NodeTable *nodes = NULL;
    EdgesetTable *edgesets = NULL;
    MutationTable *mutations = NULL;
    MutationTypeTable *mutation_types = NULL;
    static char *kwlist[] = {"nodes", "edgesets", "mutation_types", "mutations", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!O!O!", kwlist,
            &NodeTableType, &nodes,
            &EdgesetTableType, &edgesets,
            &MutationTypeTableType, &mutation_types,
            &MutationTableType, &mutations)) {
        goto out;
    }
    if (MutationGenerator_check_state(self) != 0) {
        goto out;
    }
    if (NodeTable_check_state(nodes) != 0) {
        goto out;
    }
    if (EdgesetTable_check_state(edgesets) != 0) {
        goto out;
    }
    if (MutationTypeTable_check_state(mutation_types) != 0) {
        goto out;
    }
    if (MutationTable_check_state(mutations) != 0) {
        goto out;
    }
    err = mutgen_generate_tables_tmp(self->mutgen, nodes->node_table,
            edgesets->edgeset_table);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    err = mutgen_populate_tables(self->mutgen, mutation_types->mutation_type_table,
            mutations->mutation_table);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyMemberDef MutationGenerator_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef MutationGenerator_methods[] = {
    {"get_mutation_rate", (PyCFunction) MutationGenerator_get_mutation_rate,
        METH_NOARGS, "Returns the mutation rate for this mutation generator."},
    {"generate", (PyCFunction) MutationGenerator_generate,
        METH_VARARGS|METH_KEYWORDS,
        "Generate mutations and write to the specified table."},
    {NULL}  /* Sentinel */
};

static PyTypeObject MutationGeneratorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.MutationGenerator",             /* tp_name */
    sizeof(MutationGenerator),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)MutationGenerator_dealloc, /* tp_dealloc */
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
    "MutationGenerator objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    MutationGenerator_methods,             /* tp_methods */
    MutationGenerator_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)MutationGenerator_init,      /* tp_init */
};
/*===================================================================
 * RecombinationMap
 *===================================================================
 */

static int
RecombinationMap_check_recomb_map(RecombinationMap *self)
{
    int ret = 0;
    if (self->recomb_map == NULL) {
        PyErr_SetString(PyExc_ValueError, "recomb_map not initialised");
        ret = -1;
    }
    return ret;
}

static void
RecombinationMap_dealloc(RecombinationMap* self)
{
    if (self->recomb_map != NULL) {
        recomb_map_free(self->recomb_map);
        PyMem_Free(self->recomb_map);
        self->recomb_map = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
RecombinationMap_init(RecombinationMap *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"num_loci", "positions", "rates", NULL};
    Py_ssize_t size, j;
    PyObject *py_positions = NULL;
    PyObject *py_rates = NULL;
    double *positions = NULL;
    double *rates = NULL;
    unsigned int num_loci = 0;
    PyObject *item;

    self->recomb_map = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "IO!O!", kwlist,
            &num_loci, &PyList_Type, &py_positions, &PyList_Type,
            &py_rates)) {
        goto out;
    }
    if (PyList_Size(py_positions) != PyList_Size(py_rates)) {
        PyErr_SetString(PyExc_ValueError,
            "positions and rates list must be the same length");
        goto out;
    }
    size = PyList_Size(py_positions);
    positions = PyMem_Malloc(size * sizeof(double));
    rates = PyMem_Malloc(size * sizeof(double));
    if (positions == NULL || rates == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    for (j = 0; j < size; j++) {
        item = PyList_GetItem(py_positions, j);
        if (!PyNumber_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "position must be a number");
            goto out;
        }
        positions[j] = PyFloat_AsDouble(item);
        item = PyList_GetItem(py_rates, j);
        if (!PyNumber_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "rates must be a number");
            goto out;
        }
        rates[j] = PyFloat_AsDouble(item);
    }
    self->recomb_map = PyMem_Malloc(sizeof(recomb_map_t));
    if (self->recomb_map == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = recomb_map_alloc(self->recomb_map, (uint32_t) num_loci,
            positions[size - 1], positions, rates, size);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    if (positions != NULL) {
        PyMem_Free(positions);
    }
    if (rates != NULL) {
        PyMem_Free(rates);
    }
    return ret;
}

static PyObject *
RecombinationMap_genetic_to_physical(RecombinationMap *self, PyObject *args)
{
    PyObject *ret = NULL;
    double genetic_x, physical_x;
    uint32_t num_loci;

    if (RecombinationMap_check_recomb_map(self) != 0) {
        goto out;
    }
    num_loci = recomb_map_get_num_loci(self->recomb_map);
    if (!PyArg_ParseTuple(args, "d", &genetic_x)) {
        goto out;
    }
    if (genetic_x < 0 || genetic_x > num_loci) {
        PyErr_SetString(PyExc_ValueError,
                "coordinates must be 0 <= x <= num_loci");
        goto out;
    }
    physical_x = recomb_map_genetic_to_phys(self->recomb_map, genetic_x);
    ret = Py_BuildValue("d", physical_x);
out:
    return ret;
}

static PyObject *
RecombinationMap_physical_to_genetic(RecombinationMap *self, PyObject *args)
{
    PyObject *ret = NULL;
    double genetic_x, physical_x, sequence_length;

    if (RecombinationMap_check_recomb_map(self) != 0) {
        goto out;
    }
    sequence_length = recomb_map_get_sequence_length(self->recomb_map);
    if (!PyArg_ParseTuple(args, "d", &physical_x)) {
        goto out;
    }
    if (physical_x < 0 || physical_x > sequence_length) {
        PyErr_SetString(PyExc_ValueError,
            "coordinates must be 0 <= x <= sequence_length");
        goto out;
    }
    genetic_x = recomb_map_phys_to_genetic(self->recomb_map, physical_x);
    ret = Py_BuildValue("d", genetic_x);
out:
    return ret;
}

static PyObject *
RecombinationMap_get_per_locus_recombination_rate(RecombinationMap *self)
{
    PyObject *ret = NULL;

    if (RecombinationMap_check_recomb_map(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d",
        recomb_map_get_per_locus_recombination_rate(self->recomb_map));
out:
    return ret;
}

static PyObject *
RecombinationMap_get_total_recombination_rate(RecombinationMap *self)
{
    PyObject *ret = NULL;

    if (RecombinationMap_check_recomb_map(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d",
        recomb_map_get_total_recombination_rate(self->recomb_map));
out:
    return ret;
}

static PyObject *
RecombinationMap_get_num_loci(RecombinationMap *self)
{
    PyObject *ret = NULL;

    if (RecombinationMap_check_recomb_map(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n",
        (Py_ssize_t) recomb_map_get_num_loci(self->recomb_map));
out:
    return ret;
}

static PyObject *
RecombinationMap_get_size(RecombinationMap *self)
{
    PyObject *ret = NULL;

    if (RecombinationMap_check_recomb_map(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n",
        (Py_ssize_t) recomb_map_get_size(self->recomb_map));
out:
    return ret;
}


static PyObject *
RecombinationMap_get_positions(RecombinationMap *self)
{
    PyObject *ret = NULL;
    double *positions = NULL;
    size_t size;
    int err;

    if (RecombinationMap_check_recomb_map(self) != 0) {
        goto out;
    }
    size = recomb_map_get_size(self->recomb_map);
    positions = PyMem_Malloc(size * sizeof(double));
    if (positions == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = recomb_map_get_positions(self->recomb_map, positions);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = convert_float_list(positions, size);
out:
    if (positions != NULL) {
        PyMem_Free(positions);
    }
    return ret;
}

static PyObject *
RecombinationMap_get_rates(RecombinationMap *self)
{
    PyObject *ret = NULL;
    double *rates = NULL;
    size_t size;
    int err;

    if (RecombinationMap_check_recomb_map(self) != 0) {
        goto out;
    }
    size = recomb_map_get_size(self->recomb_map);
    rates = PyMem_Malloc(size * sizeof(double));
    if (rates == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = recomb_map_get_rates(self->recomb_map, rates);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = convert_float_list(rates, size);
out:
    if (rates != NULL) {
        PyMem_Free(rates);
    }
    return ret;
}

static PyMemberDef RecombinationMap_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef RecombinationMap_methods[] = {
    {"genetic_to_physical", (PyCFunction) RecombinationMap_genetic_to_physical,
        METH_VARARGS, "Converts the specified value into physical coordinates."},
    {"physical_to_genetic", (PyCFunction) RecombinationMap_physical_to_genetic,
        METH_VARARGS, "Converts the specified value into genetic coordinates."},
    {"get_total_recombination_rate",
        (PyCFunction) RecombinationMap_get_total_recombination_rate, METH_NOARGS,
        "Returns the total product of physical distance times recombination rate"},
    {"get_per_locus_recombination_rate",
        (PyCFunction) RecombinationMap_get_per_locus_recombination_rate,
        METH_NOARGS,
        "Returns the recombination rate between loci implied by this map"},
    {"get_num_loci", (PyCFunction) RecombinationMap_get_num_loci, METH_NOARGS,
        "Returns the number discrete loci in the genetic map."},
    {"get_size", (PyCFunction) RecombinationMap_get_size, METH_NOARGS,
        "Returns the number of physical  positions in this map."},
    {"get_positions",
        (PyCFunction) RecombinationMap_get_positions, METH_NOARGS,
        "Returns the positions in this recombination map."},
    {"get_rates",
        (PyCFunction) RecombinationMap_get_rates, METH_NOARGS,
        "Returns the rates in this recombination map."},
    {NULL}  /* Sentinel */
};

static PyTypeObject RecombinationMapType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.RecombinationMap",             /* tp_name */
    sizeof(RecombinationMap),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)RecombinationMap_dealloc, /* tp_dealloc */
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
    "RecombinationMap objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    RecombinationMap_methods,             /* tp_methods */
    RecombinationMap_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)RecombinationMap_init,      /* tp_init */
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
    int ret = -1;
    int err;

    self->tree_sequence = NULL;
    self->tree_sequence = PyMem_Malloc(sizeof(tree_sequence_t));
    if (self->tree_sequence == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tree_sequence_initialise(self->tree_sequence);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
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
    int flags = 0;
    static char *kwlist[] = {"path", "zlib_compression", NULL};

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|i", kwlist,
                &path, &zlib_compression)) {
        goto out;
    }
    if (zlib_compression) {
        flags = MSP_DUMP_ZLIB_COMPRESSION;
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
TreeSequence_load_tables(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    int err;
    PyObject *ret = NULL;
    NodeTable *py_nodes = NULL;
    EdgesetTable *py_edgesets = NULL;
    MigrationTable *py_migrations = NULL;
    MutationTypeTable *py_mutation_types = NULL;
    MutationTable *py_mutations = NULL;
    PyObject *py_provenance_strings = NULL;
    Py_ssize_t num_provenance_strings = 0;
    char **provenance_strings = NULL;
    node_table_t *nodes = NULL;
    edgeset_table_t *edgesets = NULL;
    migration_table_t *migrations = NULL;
    mutation_table_t *mutations = NULL;
    mutation_type_table_t *mutation_types = NULL;

    static char *kwlist[] = {"nodes", "edgesets", "migrations",
        "mutation_types", "mutations", "provenance_strings", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|O!O!O!O!", kwlist,
            &NodeTableType, &py_nodes,
            &EdgesetTableType, &py_edgesets,
            &MigrationTableType, &py_migrations,
            &MutationTypeTableType, &py_mutation_types,
            &MutationTableType, &py_mutations,
            &PyList_Type, &py_provenance_strings)) {
        goto out;
    }
    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (NodeTable_check_state(py_nodes) != 0) {
        goto out;
    }
    nodes = py_nodes->node_table;
    if (EdgesetTable_check_state(py_edgesets) != 0) {
        goto out;
    }
    edgesets = py_edgesets->edgeset_table;
    if (py_migrations != NULL) {
        if (MigrationTable_check_state(py_migrations) != 0) {
            goto out;
        }
        migrations = py_migrations->migration_table;
    }
    if (py_mutation_types != NULL) {
        if (MutationTypeTable_check_state(py_mutation_types) != 0) {
            goto out;
        }
        mutation_types = py_mutation_types->mutation_type_table;
    }
    if (py_mutations != NULL) {
        if (MutationTable_check_state(py_mutations) != 0) {
            goto out;
        }
        mutations = py_mutations->mutation_table;
    }
    num_provenance_strings = 0;
    if (py_provenance_strings != NULL) {
        if (parse_provenance_strings(py_provenance_strings, &num_provenance_strings,
                    &provenance_strings) != 0) {
            goto out;
        }
    }
    if ((mutations == NULL) != (mutation_types == NULL)) {
        PyErr_SetString(PyExc_TypeError, "Must specify both mutations and mutation types");
        goto out;
    }
    err = tree_sequence_load_tables_tmp(self->tree_sequence,
        nodes, edgesets, migrations, mutation_types, mutations,
        num_provenance_strings, provenance_strings);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
TreeSequence_dump_tables(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    int err;
    PyObject *ret = NULL;
    NodeTable *py_nodes = NULL;
    EdgesetTable *py_edgesets = NULL;
    MigrationTable *py_migrations = NULL;
    MutationTypeTable *py_mutation_types = NULL;
    MutationTable *py_mutations = NULL;
    node_table_t *nodes = NULL;
    edgeset_table_t *edgesets = NULL;
    migration_table_t *migrations = NULL;
    mutation_type_table_t *mutation_types = NULL;
    mutation_table_t *mutations = NULL;
    size_t num_provenance_strings = 0;
    char **provenance_strings = NULL;
    static char *kwlist[] = {"nodes", "edgesets", "migrations",
        "mutation_types", "mutations", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|O!O!O!", kwlist,
            &NodeTableType, &py_nodes,
            &EdgesetTableType, &py_edgesets,
            &MigrationTableType, &py_migrations,
            &MutationTypeTableType, &py_mutation_types,
            &MutationTableType, &py_mutations)) {
        goto out;
    }
    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (NodeTable_check_state(py_nodes) != 0) {
        goto out;
    }
    nodes = py_nodes->node_table;
    if (EdgesetTable_check_state(py_edgesets) != 0) {
        goto out;
    }
    edgesets = py_edgesets->edgeset_table;
    if (py_migrations != NULL) {
        if (MigrationTable_check_state(py_migrations) != 0) {
            goto out;
        }
        migrations = py_migrations->migration_table;
    }
    if (py_mutation_types != NULL) {
        if (MutationTypeTable_check_state(py_mutation_types) != 0) {
            goto out;
        }
        mutation_types = py_mutation_types->mutation_type_table;
    }
    if (py_mutations != NULL) {
        if (MutationTable_check_state(py_mutations) != 0) {
            goto out;
        }
        mutations = py_mutations->mutation_table;
    }
    if ((mutations == NULL) != (mutation_types == NULL)) {
        PyErr_SetString(PyExc_TypeError, "Must specify both mutations and mutation types");
        goto out;
    }
    err = tree_sequence_dump_tables_tmp(self->tree_sequence,
        nodes, edgesets, migrations, mutation_types, mutations,
        &num_provenance_strings, &provenance_strings);
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
    PyObject *ret = NULL;
    static char *kwlist[] = {"path", NULL};

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &path)) {
        goto out;
    }
    /* Silence the low-level error reporting HDF5 */
    if (H5Eset_auto(H5E_DEFAULT, NULL, NULL) < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Error silencing HDF5 errors");
        goto out;
    }
    err = tree_sequence_load(self->tree_sequence, path, flags);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
TreeSequence_get_provenance_strings(TreeSequence *self)
{
    int err;
    PyObject *ret = NULL;
    size_t num_provenance_strings;
    char **provenance_strings;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    err = tree_sequence_get_provenance_strings(self->tree_sequence,
            &num_provenance_strings, &provenance_strings);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = convert_string_list(provenance_strings, num_provenance_strings);
out:
    return ret;
}

static PyObject *
TreeSequence_get_record(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    edgeset_t edgeset;
    node_t node;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tree_sequence_get_num_edgesets(
        self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tree_sequence_get_edgeset(self->tree_sequence, (size_t) record_index,
            &edgeset);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    err = tree_sequence_get_node(self->tree_sequence, edgeset.parent, &node);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_coalescence_record_tmp(&node, &edgeset);
out:
    return ret;
}

static PyObject *
TreeSequence_get_node(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    node_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tree_sequence_get_num_nodes(self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tree_sequence_get_node(self->tree_sequence, (size_t) record_index, &record);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_node(&record);
out:
    return ret;
}

static PyObject *
TreeSequence_get_edgeset(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    edgeset_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tree_sequence_get_num_edgesets(self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tree_sequence_get_edgeset(self->tree_sequence, (size_t) record_index, &record);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_edgeset(&record);
out:
    return ret;
}

static PyObject *
TreeSequence_get_migration(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    migration_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tree_sequence_get_num_migrations(
        self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tree_sequence_get_migration(self->tree_sequence,
            (size_t) record_index, &record);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_migration(&record);
out:
    return ret;
}

static PyObject *
TreeSequence_get_mutation_type(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    mutation_type_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tree_sequence_get_num_mutation_types(
        self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tree_sequence_get_mutation_type(self->tree_sequence,
            (size_t) record_index, &record);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_mutation_type(&record);
out:
    return ret;
}

static PyObject *
TreeSequence_get_mutation(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    mutation_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tree_sequence_get_num_mutations(
        self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tree_sequence_get_mutation(self->tree_sequence,
            (size_t) record_index, &record);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_mutation(&record);
out:
    return ret;
}


static PyObject *
TreeSequence_get_num_edgesets(TreeSequence *self, PyObject *args)
{
    PyObject *ret = NULL;
    size_t num_records;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_records = tree_sequence_get_num_edgesets(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_records);
out:
    return ret;
}

static PyObject *
TreeSequence_get_num_migrations(TreeSequence *self, PyObject *args)
{
    PyObject *ret = NULL;
    size_t num_records;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_records = tree_sequence_get_num_migrations(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_records);
out:
    return ret;
}



static PyObject *
TreeSequence_get_num_trees(TreeSequence *self, PyObject *args)
{
    PyObject *ret = NULL;
    size_t num_trees;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_trees = tree_sequence_get_num_trees(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_trees);
out:
    return ret;
}

static PyObject *
TreeSequence_get_sequence_length(TreeSequence  *self)
{
    PyObject *ret = NULL;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d",
        tree_sequence_get_sequence_length(self->tree_sequence));
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

/* TODO refactor this to be get_node() */
static PyObject *
TreeSequence_get_sample(TreeSequence *self, PyObject *args)
{
    PyObject *ret = NULL;
    unsigned int u;
    node_t node;
    int err;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &u)) {
        goto out;
    }
    if (u >= tree_sequence_get_sample_size(self->tree_sequence)) {
        PyErr_SetString(PyExc_IndexError, "out of bounds");
        goto out;
    }
    err = tree_sequence_get_node(self->tree_sequence, u, &node);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("id", (int) node.population, node.time);
out:
    return ret;
}

static PyObject *
TreeSequence_get_pairwise_diversity(TreeSequence *self, PyObject *args,
        PyObject *kwds)
{
    PyObject *ret = NULL;
    PyObject *py_samples = NULL;
    static char *kwlist[] = {"samples", NULL};
    node_id_t *samples = NULL;
    size_t num_samples = 0;
    double pi;
    int err;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
            &PyList_Type, &py_samples)) {
        goto out;
    }
    if (parse_sample_ids(py_samples, self->tree_sequence, &num_samples, &samples) != 0) {
        goto out;
    }
    err = tree_sequence_get_pairwise_diversity(
        self->tree_sequence, samples, (uint32_t) num_samples, &pi);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("d", pi);
out:
    if (samples != NULL) {
        PyMem_Free(samples);
    }

    return ret;
}

/* Forward declaration */
static PyTypeObject TreeSequenceType;
static PyObject *
TreeSequence_simplify(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    PyObject *py_samples = NULL;
    static char *kwlist[] = {"output", "samples", "filter_root_mutations", NULL};
    node_id_t *samples = NULL;
    size_t num_samples = 0;
    TreeSequence *output = NULL;
    int filter_root_mutations = 1;
    int flags = 0;
    int err;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|i", kwlist,
            &TreeSequenceType, &output,
            &PyList_Type, &py_samples,
            &filter_root_mutations)) {
        goto out;
    }
    if (TreeSequence_check_tree_sequence(output) != 0) {
        goto out;
    }
    if (parse_sample_ids(py_samples, self->tree_sequence, &num_samples, &samples) != 0) {
        goto out;
    }
    if (filter_root_mutations) {
        flags |= MSP_FILTER_ROOT_MUTATIONS;
    }
    err = tree_sequence_simplify(
        self->tree_sequence, samples, (uint32_t) num_samples, flags, output->tree_sequence);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    if (samples != NULL) {
        PyMem_Free(samples);
    }
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
TreeSequence_get_num_mutation_types(TreeSequence  *self)
{
    PyObject *ret = NULL;
    size_t num_mutation_types;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_mutation_types = tree_sequence_get_num_mutation_types(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_mutation_types);
out:
    return ret;
}

static PyMemberDef TreeSequence_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef TreeSequence_methods[] = {
    {"dump", (PyCFunction) TreeSequence_dump,
        METH_VARARGS|METH_KEYWORDS,
        "Writes the tree sequence out to the specified path."},
    {"load", (PyCFunction) TreeSequence_load,
        METH_VARARGS|METH_KEYWORDS,
        "Loads a tree sequence from the specified path."},
    {"load_tables", (PyCFunction) TreeSequence_load_tables,
        METH_VARARGS|METH_KEYWORDS,
        "Loads a tree sequence from the specified set of tables"},
    {"dump_tables", (PyCFunction) TreeSequence_dump_tables,
        METH_VARARGS|METH_KEYWORDS,
        "Dumps the tree sequence to the specified set of tables"},
    {"get_provenance_strings", (PyCFunction) TreeSequence_get_provenance_strings,
        METH_NOARGS, "Returns the list of provenance strings."},
    {"get_record", (PyCFunction) TreeSequence_get_record, METH_VARARGS,
        "Returns the record at the specified index."},
    {"get_node",
        (PyCFunction) TreeSequence_get_node, METH_VARARGS,
        "Returns the node record at the specified index."},
    {"get_edgeset",
        (PyCFunction) TreeSequence_get_edgeset, METH_VARARGS,
        "Returns the edgeset record at the specified index."},
    {"get_migration",
        (PyCFunction) TreeSequence_get_migration, METH_VARARGS,
        "Returns the migration record at the specified index."},
    {"get_mutation_type",
        (PyCFunction) TreeSequence_get_mutation_type, METH_VARARGS,
        "Returns the mutation type record at the specified index."},
    {"get_mutation",
        (PyCFunction) TreeSequence_get_mutation, METH_VARARGS,
        "Returns the mutation record at the specified index."},
    {"get_num_edgesets", (PyCFunction) TreeSequence_get_num_edgesets,
        METH_NOARGS, "Returns the number of coalescence records." },
    {"get_num_migrations", (PyCFunction) TreeSequence_get_num_migrations,
        METH_NOARGS, "Returns the number of migration records." },
    {"get_num_trees", (PyCFunction) TreeSequence_get_num_trees,
        METH_NOARGS, "Returns the number of trees in the tree sequence." },
    {"get_sequence_length", (PyCFunction) TreeSequence_get_sequence_length,
        METH_NOARGS, "Returns the sequence length in bases." },
    {"get_num_mutation_types", (PyCFunction) TreeSequence_get_num_mutation_types,
        METH_NOARGS, "Returns the number of mutation types" },
    {"get_num_mutations", (PyCFunction) TreeSequence_get_num_mutations, METH_NOARGS,
        "Returns the number of mutations" },
    {"get_num_nodes", (PyCFunction) TreeSequence_get_num_nodes, METH_NOARGS,
        "Returns the number of unique nodes in the tree sequence." },
    {"get_sample_size", (PyCFunction) TreeSequence_get_sample_size, METH_NOARGS,
        "Returns the sample size" },
    {"get_sample", (PyCFunction) TreeSequence_get_sample, METH_VARARGS,
        "Returns a dictionary describing the specified sample." },
    {"get_pairwise_diversity",
        (PyCFunction) TreeSequence_get_pairwise_diversity,
        METH_VARARGS|METH_KEYWORDS, "Returns the average pairwise diversity." },
    {"simplify", (PyCFunction) TreeSequence_simplify,
        METH_VARARGS|METH_KEYWORDS,
        "Returns a simplified version of this tree sequence."},
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
        PyErr_SetString(PyExc_RuntimeError, "sparse_tree not initialised");
        ret = -1;
    }
    return ret;
}

static int
SparseTree_check_bounds(SparseTree *self, int node)
{
    int ret = 0;
    if (node < 0 || node >= self->sparse_tree->num_nodes) {
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
    Py_XDECREF(self->tree_sequence);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
SparseTree_init(SparseTree *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"tree_sequence", "flags", "tracked_leaves",
        NULL};
    PyObject *py_tracked_leaves = NULL;
    TreeSequence *tree_sequence = NULL;
    node_id_t *tracked_leaves = NULL;
    int flags = 0;
    uint32_t j, n, num_tracked_leaves;
    PyObject *item;

    self->sparse_tree = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|iO!", kwlist,
            &TreeSequenceType, &tree_sequence,
            &flags, &PyList_Type, &py_tracked_leaves)) {
        goto out;
    }
    self->tree_sequence = tree_sequence;
    Py_INCREF(self->tree_sequence);
    if (TreeSequence_check_tree_sequence(tree_sequence) != 0) {
        goto out;
    }
    n = tree_sequence_get_sample_size(tree_sequence->tree_sequence);
    num_tracked_leaves = 0;
    if (py_tracked_leaves != NULL) {
        if (!flags & MSP_LEAF_COUNTS) {
            PyErr_SetString(PyExc_ValueError,
                "Cannot specified tracked_leaves without count_leaves flag");
            goto out;
        }
        num_tracked_leaves = PyList_Size(py_tracked_leaves);
    }
    tracked_leaves = PyMem_Malloc(num_tracked_leaves * sizeof(node_id_t));
    if (tracked_leaves == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    for (j = 0; j < num_tracked_leaves; j++) {
        item = PyList_GetItem(py_tracked_leaves, j);
        if (!PyNumber_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "leaf must be a number");
            goto out;
        }
        tracked_leaves[j] = (node_id_t) PyLong_AsLong(item);
        if (tracked_leaves[j] >= n) {
            PyErr_SetString(PyExc_ValueError, "leaves must be < sample_size");
            goto out;
        }
    }
    self->sparse_tree = PyMem_Malloc(sizeof(sparse_tree_t));
    if (self->sparse_tree == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = sparse_tree_alloc(self->sparse_tree, tree_sequence->tree_sequence,
           flags);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    if (flags & MSP_LEAF_COUNTS) {
        err = sparse_tree_set_tracked_leaves(self->sparse_tree, num_tracked_leaves,
                tracked_leaves);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    ret = 0;
out:
    if (tracked_leaves != NULL) {
        PyMem_Free(tracked_leaves);
    }
    return ret;
}

static PyObject *
SparseTree_free(SparseTree *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    /* This method is need because we have dangling references to
     * trees after a for loop and we can't run set_mutations.
     */
    sparse_tree_free(self->sparse_tree);
    PyMem_Free(self->sparse_tree);
    self->sparse_tree = NULL;
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
SparseTree_get_sample_size(SparseTree *self)
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
SparseTree_get_num_nodes(SparseTree *self)
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
SparseTree_get_index(SparseTree *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sparse_tree->index);
out:
    return ret;
}

static PyObject *
SparseTree_get_root(SparseTree *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("i", (int) self->sparse_tree->root);
out:
    return ret;
}

static PyObject *
SparseTree_get_left(SparseTree *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sparse_tree->left);
out:
    return ret;
}

static PyObject *
SparseTree_get_right(SparseTree *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sparse_tree->right);
out:
    return ret;
}

static PyObject *
SparseTree_get_flags(SparseTree *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("i", self->sparse_tree->flags);
out:
    return ret;
}

static int
SparseTree_get_node_argument(SparseTree *self, PyObject *args, int *node)
{
    int ret = -1;
    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", node)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, *node)) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
SparseTree_get_parent(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    node_id_t parent;
    int node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    parent = self->sparse_tree->parent[node];
    ret = Py_BuildValue("i", (int) parent);
out:
    return ret;
}

static PyObject *
SparseTree_get_population(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    population_id_t population;
    int node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    population = self->sparse_tree->population[node];
    ret = Py_BuildValue("i", (int) population);
out:
    return ret;
}

static PyObject *
SparseTree_get_time(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    double time;
    int node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
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
    node_id_t *children;
    size_t num_children;
    int err, node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    err = sparse_tree_get_children(self->sparse_tree,
            (node_id_t) node, &num_children, &children);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = convert_node_id_list(children, num_children);
out:
    return ret;
}

static PyObject *
SparseTree_get_mrca(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    int err;
    node_id_t mrca;
    int u, v;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "ii", &u, &v)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, u)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, v)) {
        goto out;
    }
    err = sparse_tree_get_mrca(self->sparse_tree, (node_id_t) u,
            (node_id_t) v, &mrca);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("i", (int) mrca);
out:
    return ret;
}

static PyObject *
SparseTree_get_num_leaves(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    size_t num_leaves;
    int err, node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    err = sparse_tree_get_num_leaves(self->sparse_tree, (node_id_t) node,
            &num_leaves);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("I", (unsigned int) num_leaves);
out:
    return ret;
}

static PyObject *
SparseTree_get_num_tracked_leaves(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    size_t num_tracked_leaves;
    int err, node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    err = sparse_tree_get_num_tracked_leaves(self->sparse_tree, (node_id_t) node,
            &num_tracked_leaves);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("I", (unsigned int) num_tracked_leaves);
out:
    return ret;
}

static PyObject *
SparseTree_get_mutations(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = convert_mutations(self->sparse_tree->mutations,
            self->sparse_tree->num_mutations);
out:
    return ret;
}

static PyObject *
SparseTree_get_num_mutations(SparseTree  *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sparse_tree->num_mutations);
out:
    return ret;
}


static PyMemberDef SparseTree_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef SparseTree_methods[] = {
    {"free", (PyCFunction) SparseTree_free, METH_NOARGS,
            "Frees the underlying tree object." },
    {"get_num_nodes", (PyCFunction) SparseTree_get_num_nodes, METH_NOARGS,
            "Returns the number of nodes in the sparse tree." },
    {"get_sample_size", (PyCFunction) SparseTree_get_sample_size, METH_NOARGS,
            "Returns the sample size" },
    {"get_index", (PyCFunction) SparseTree_get_index, METH_NOARGS,
            "Returns the index this tree occupies within the tree sequence." },
    {"get_root", (PyCFunction) SparseTree_get_root, METH_NOARGS,
            "Returns the root of the tree." },
    {"get_left", (PyCFunction) SparseTree_get_left, METH_NOARGS,
            "Returns the left-most coordinate (inclusive)." },
    {"get_right", (PyCFunction) SparseTree_get_right, METH_NOARGS,
            "Returns the right-most coordinate (exclusive)." },
    {"get_mutations", (PyCFunction) SparseTree_get_mutations, METH_NOARGS,
            "Returns the list of mutations on this tree." },
    {"get_flags", (PyCFunction) SparseTree_get_flags, METH_NOARGS,
            "Returns the value of the flags variable." },
    {"get_num_mutations", (PyCFunction) SparseTree_get_num_mutations, METH_NOARGS,
            "Returns the number of mutations on this tree." },
    {"get_parent", (PyCFunction) SparseTree_get_parent, METH_VARARGS,
            "Returns the parent of node u" },
    {"get_time", (PyCFunction) SparseTree_get_time, METH_VARARGS,
            "Returns the time of node u" },
    {"get_population", (PyCFunction) SparseTree_get_population, METH_VARARGS,
            "Returns the population of node u" },
    {"get_children", (PyCFunction) SparseTree_get_children, METH_VARARGS,
            "Returns the children of node u" },
    {"get_mrca", (PyCFunction) SparseTree_get_mrca, METH_VARARGS,
            "Returns the MRCA of nodes u and v" },
    {"get_num_leaves", (PyCFunction) SparseTree_get_num_leaves, METH_VARARGS,
            "Returns the number of leaves below node u." },
    {"get_num_tracked_leaves", (PyCFunction) SparseTree_get_num_tracked_leaves,
            METH_VARARGS,
            "Returns the number of tracked leaves below node u." },
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
    PyObject *children = NULL;
    int err;
    double length;
    size_t list_size, j;
    node_record_t *records_out, *records_in, *record;

    if (TreeDiffIterator_check_state(self) != 0) {
        goto out;
    }
    err = tree_diff_iterator_next(self->tree_diff_iterator, &length,
            &records_out, &records_in);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    if (err == 1) {
        /* out records */
        record = records_out;
        list_size = 0;
        while (record != NULL) {
            list_size++;
            record = record->next;
        }
        out_list = PyList_New(list_size);
        if (out_list == NULL) {
            goto out;
        }
        record = records_out;
        j = 0;
        while (record != NULL) {
            children = convert_node_id_list(record->children, record->num_children);
            if (children == NULL) {
                goto out;
            }
            value = Py_BuildValue("IOd", (unsigned int) record->node,
                    children, record->time);
            Py_DECREF(children);
            if (value == NULL) {
                goto out;
            }
            PyList_SET_ITEM(out_list, j, value);
            record = record->next;
            j++;
        }
        /* in records */
        record = records_in;
        list_size = 0;
        while (record != NULL) {
            list_size++;
            record = record->next;
        }
        in_list = PyList_New(list_size);
        if (in_list == NULL) {
            goto out;
        }
        record = records_in;
        j = 0;
        while (record != NULL) {
            children = convert_node_id_list(record->children, record->num_children);
            if (children == NULL) {
                goto out;
            }
            value = Py_BuildValue("IOd", (unsigned int) record->node,
                    children, record->time);
            Py_DECREF(children);
            if (value == NULL) {
                goto out;
            }
            PyList_SET_ITEM(in_list, j, value);
            record = record->next;
            j++;
        }
        ret = Py_BuildValue("dOO", length, out_list, in_list);
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
 * LeafListIterator
 *===================================================================
 */

static int
LeafListIterator_check_state(LeafListIterator *self)
{
    int ret = 0;
    if (self->sparse_tree == NULL) {
        PyErr_SetString(PyExc_SystemError, "iterator not initialised");
        ret = -1;
    }
    return ret;
}

static void
LeafListIterator_dealloc(LeafListIterator* self)
{
    Py_XDECREF(self->sparse_tree);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
LeafListIterator_init(LeafListIterator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"sparse_tree", "node", NULL};
    unsigned int node = 0;
    SparseTree *sparse_tree = NULL;

    self->sparse_tree = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!I", kwlist,
            &SparseTreeType, &sparse_tree, &node)) {
        goto out;
    }
    self->sparse_tree = sparse_tree;
    Py_INCREF(self->sparse_tree);
    if (SparseTree_check_sparse_tree(sparse_tree) != 0) {
        goto out;
    }
    if (SparseTree_check_bounds(self->sparse_tree, node)) {
        goto out;
    }
    err = sparse_tree_get_leaf_list(self->sparse_tree->sparse_tree,
            (uint32_t) node, &self->head, &self->tail);
    self->next = self->head;
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
LeafListIterator_next(LeafListIterator  *self)
{
    PyObject *ret = NULL;

    if (LeafListIterator_check_state(self) != 0) {
        goto out;
    }
    if (self->next != NULL) {
        ret = Py_BuildValue("I", (unsigned int) self->next->node);
        if (ret == NULL) {
            goto out;
        }
        /* Get the next value */
        if (self->next == self->tail) {
            self->next = NULL;
        } else {
            self->next = self->next->next;
        }
    }
out:
    return ret;
}

static PyMemberDef LeafListIterator_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef LeafListIterator_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject LeafListIteratorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.LeafListIterator",             /* tp_name */
    sizeof(LeafListIterator),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)LeafListIterator_dealloc, /* tp_dealloc */
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
    "LeafListIterator objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    PyObject_SelfIter,                    /* tp_iter */
    (iternextfunc) LeafListIterator_next, /* tp_iternext */
    LeafListIterator_methods,             /* tp_methods */
    LeafListIterator_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LeafListIterator_init,      /* tp_init */
};


/*===================================================================
 * SparseTreeIterator
 *===================================================================
 */

static int
SparseTreeIterator_check_state(SparseTreeIterator *self)
{
    int ret = 0;
    return ret;
}

static void
SparseTreeIterator_dealloc(SparseTreeIterator* self)
{
    Py_XDECREF(self->sparse_tree);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
SparseTreeIterator_init(SparseTreeIterator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    static char *kwlist[] = {"sparse_tree", NULL};
    SparseTree *sparse_tree;

    self->first = 1;
    self->sparse_tree = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
            &SparseTreeType, &sparse_tree)) {
        goto out;
    }
    self->sparse_tree = sparse_tree;
    Py_INCREF(self->sparse_tree);
    if (SparseTree_check_sparse_tree(self->sparse_tree) != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
SparseTreeIterator_next(SparseTreeIterator  *self)
{
    PyObject *ret = NULL;
    int err;

    if (SparseTreeIterator_check_state(self) != 0) {
        goto out;
    }

    if (self->first) {
        err = sparse_tree_first(self->sparse_tree->sparse_tree);
        self->first = 0;
    } else {
        err = sparse_tree_next(self->sparse_tree->sparse_tree);
    }
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    if (err == 1) {
        ret = (PyObject *) self->sparse_tree;
        Py_INCREF(ret);
    }
out:
    return ret;
}

static PyMemberDef SparseTreeIterator_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef SparseTreeIterator_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject SparseTreeIteratorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.SparseTreeIterator",             /* tp_name */
    sizeof(SparseTreeIterator),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)SparseTreeIterator_dealloc, /* tp_dealloc */
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
    "SparseTreeIterator objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    PyObject_SelfIter,                    /* tp_iter */
    (iternextfunc) SparseTreeIterator_next, /* tp_iternext */
    SparseTreeIterator_methods,             /* tp_methods */
    SparseTreeIterator_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)SparseTreeIterator_init,      /* tp_init */
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
    static char *kwlist[] = {"tree_sequence", "precision", "Ne", NULL};
    int precision = 3;
    double Ne = 0.25; /* default to 1/4 for coalescent time units. */
    TreeSequence *tree_sequence;

    self->newick_converter = NULL;
    self->tree_sequence = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|id", kwlist,
            &TreeSequenceType, &tree_sequence, &precision, &Ne)) {
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
            self->tree_sequence->tree_sequence, (size_t) precision, Ne);
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
    double length;
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
        ret = Py_BuildValue("ds", length, tree);
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
 * VcfConverter
 *===================================================================
 */

static int
VcfConverter_check_state(VcfConverter *self)
{
    int ret = 0;
    if (self->vcf_converter == NULL) {
        PyErr_SetString(PyExc_SystemError, "converter not initialised");
        ret = -1;
    }
    return ret;
}

static void
VcfConverter_dealloc(VcfConverter* self)
{
    if (self->vcf_converter != NULL) {
        vcf_converter_free(self->vcf_converter);
        PyMem_Free(self->vcf_converter);
        self->vcf_converter = NULL;
    }
    Py_XDECREF(self->tree_sequence);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
VcfConverter_init(VcfConverter *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"tree_sequence", "ploidy", NULL};
    unsigned int ploidy = 1;
    TreeSequence *tree_sequence;

    self->vcf_converter = NULL;
    self->tree_sequence = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|I", kwlist,
            &TreeSequenceType, &tree_sequence, &ploidy)) {
        goto out;
    }
    self->tree_sequence = tree_sequence;
    Py_INCREF(self->tree_sequence);
    if (TreeSequence_check_tree_sequence(self->tree_sequence) != 0) {
        goto out;
    }
    if (ploidy < 1) {
        PyErr_SetString(PyExc_ValueError, "Ploidy must be >= 1");
        goto out;
    }
    self->vcf_converter = PyMem_Malloc(sizeof(vcf_converter_t));
    if (self->vcf_converter == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = vcf_converter_alloc(self->vcf_converter,
            self->tree_sequence->tree_sequence, ploidy);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
VcfConverter_next(VcfConverter  *self)
{
    PyObject *ret = NULL;
    char *record;
    int err;

    if (VcfConverter_check_state(self) != 0) {
        goto out;
    }
    err = vcf_converter_next(self->vcf_converter, &record);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    if (err == 1) {
        ret = Py_BuildValue("s", record);
    }
out:
    return ret;
}

static PyObject *
VcfConverter_get_header(VcfConverter *self)
{
    PyObject *ret = NULL;
    int err;
    char *header;

    if (VcfConverter_check_state(self) != 0) {
        goto out;
    }
    err = vcf_converter_get_header(self->vcf_converter, &header);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("s", header);
out:
    return ret;
}

static PyMemberDef VcfConverter_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef VcfConverter_methods[] = {
    {"get_header", (PyCFunction) VcfConverter_get_header, METH_NOARGS,
            "Returns the VCF header as plain text." },
    {NULL}  /* Sentinel */
};

static PyTypeObject VcfConverterType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.VcfConverter",             /* tp_name */
    sizeof(VcfConverter),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)VcfConverter_dealloc, /* tp_dealloc */
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
    "VcfConverter objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    PyObject_SelfIter,                    /* tp_iter */
    (iternextfunc) VcfConverter_next, /* tp_iternext */
    VcfConverter_methods,             /* tp_methods */
    VcfConverter_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)VcfConverter_init,      /* tp_init */
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
 * VariantGenerator
 *===================================================================
 */

static int
VariantGenerator_check_state(VariantGenerator *self)
{
    int ret = 0;
    if (self->variant_generator == NULL) {
        PyErr_SetString(PyExc_SystemError, "converter not initialised");
        ret = -1;
    }
    return ret;
}

static void
VariantGenerator_dealloc(VariantGenerator* self)
{
    if (self->variant_generator != NULL) {
        vargen_free(self->variant_generator);
        PyMem_Free(self->variant_generator);
        self->variant_generator = NULL;
    }
    Py_XDECREF(self->tree_sequence);
    Py_XDECREF(self->genotypes_buffer);
    if (self->buffer_acquired) {
        PyBuffer_Release(&self->buffer);
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
VariantGenerator_init(VariantGenerator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"tree_sequence", "genotypes_buffer", "as_char", NULL};
    TreeSequence *tree_sequence = NULL;
    PyObject *genotypes_buffer = NULL;
    int as_char = 0;
    int flags = 0;
    size_t sample_size;

    self->variant_generator = NULL;
    self->tree_sequence = NULL;
    self->genotypes_buffer = NULL;
    self->buffer_acquired = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O|i", kwlist,
            &TreeSequenceType, &tree_sequence, &genotypes_buffer,
            &as_char)) {
        goto out;
    }
    self->tree_sequence = tree_sequence;
    Py_INCREF(self->tree_sequence);
    self->genotypes_buffer = genotypes_buffer;
    Py_INCREF(self->genotypes_buffer);
    if (TreeSequence_check_tree_sequence(self->tree_sequence) != 0) {
        goto out;
    }
    sample_size = tree_sequence_get_sample_size(
            self->tree_sequence->tree_sequence);
    if (!PyObject_CheckBuffer(genotypes_buffer)) {
        PyErr_SetString(PyExc_TypeError,
            "genotypes buffer must support the Python buffer protocol.");
        goto out;
    }
    if (PyObject_GetBuffer(genotypes_buffer, &self->buffer,
                PyBUF_SIMPLE|PyBUF_WRITABLE) != 0) {
        goto out;
    }
    self->buffer_acquired = 1;
    if (sample_size * sizeof(uint8_t) > self->buffer.len) {
        PyErr_SetString(PyExc_BufferError, "genotypes buffer is too small");
        goto out;
    }
    self->variant_generator = PyMem_Malloc(sizeof(vargen_t));
    if (self->variant_generator == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    flags = as_char? MSP_GENOTYPES_AS_CHAR: 0;
    err = vargen_alloc(self->variant_generator,
            self->tree_sequence->tree_sequence, flags);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
VariantGenerator_next(VariantGenerator *self)
{
    PyObject *ret = NULL;
    mutation_t *mutation;
    int err;
    char *genotypes = (char *) self->buffer.buf;

    if (VariantGenerator_check_state(self) != 0) {
        goto out;
    }
    err = vargen_next(self->variant_generator, &mutation, genotypes);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    if (err == 1) {
        ret = make_mutation(mutation);
    }
out:
    return ret;
}

static PyMemberDef VariantGenerator_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef VariantGenerator_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject VariantGeneratorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.VariantGenerator",             /* tp_name */
    sizeof(VariantGenerator),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)VariantGenerator_dealloc, /* tp_dealloc */
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
    "VariantGenerator objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    PyObject_SelfIter,                    /* tp_iter */
    (iternextfunc) VariantGenerator_next, /* tp_iternext */
    VariantGenerator_methods,             /* tp_methods */
    VariantGenerator_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)VariantGenerator_init,      /* tp_init */
};

/*===================================================================
 * LdCalculator
 *===================================================================
 */

static int
LdCalculator_check_state(LdCalculator *self)
{
    int ret = 0;
    if (self->ld_calc == NULL) {
        PyErr_SetString(PyExc_SystemError, "converter not initialised");
        ret = -1;
    }
    return ret;
}

static void
LdCalculator_dealloc(LdCalculator* self)
{
    if (self->ld_calc != NULL) {
        ld_calc_free(self->ld_calc);
        PyMem_Free(self->ld_calc);
        self->ld_calc = NULL;
    }
    Py_XDECREF(self->tree_sequence);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
LdCalculator_init(LdCalculator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"tree_sequence", NULL};
    TreeSequence *tree_sequence;

    self->ld_calc = NULL;
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
    self->ld_calc = PyMem_Malloc(sizeof(ld_calc_t));
    if (self->ld_calc == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->ld_calc, 0, sizeof(ld_calc_t));
    err = ld_calc_alloc(self->ld_calc, self->tree_sequence->tree_sequence);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
LdCalculator_get_r2(LdCalculator *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t a, b;
    double r2;

    if (LdCalculator_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "nn", &a, &b)) {
        goto out;
    }
    Py_BEGIN_ALLOW_THREADS
    err = ld_calc_get_r2(self->ld_calc, (size_t) a, (size_t) b, &r2);
    Py_END_ALLOW_THREADS
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("d", r2);
out:
    return ret;
}

static PyObject *
LdCalculator_get_r2_array(LdCalculator *self, PyObject *args, PyObject *kwds)
{
    int err;
    PyObject *ret = NULL;
    static char *kwlist[] = {
        "dest", "source_index", "direction", "max_mutations",
        "max_distance", NULL};
    PyObject *dest = NULL;
    Py_buffer buffer;
    Py_ssize_t source_index;
    Py_ssize_t max_mutations = -1;
    double max_distance = DBL_MAX;
    int direction = MSP_DIR_FORWARD;
    size_t num_r2_values = 0;
    int buffer_acquired = 0;

    if (LdCalculator_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "On|ind", kwlist,
            &dest, &source_index, &direction, &max_mutations, &max_distance)) {
        goto out;
    }
    if (direction != MSP_DIR_FORWARD && direction != MSP_DIR_REVERSE) {
        PyErr_SetString(PyExc_ValueError,
            "direction must be FORWARD or REVERSE");
        goto out;
    }
    if (max_distance < 0) {
        PyErr_SetString(PyExc_ValueError, "max_distance must be >= 0");
        goto out;
    }
    if (!PyObject_CheckBuffer(dest)) {
        PyErr_SetString(PyExc_TypeError,
            "dest buffer must support the Python buffer protocol.");
        goto out;
    }
    if (PyObject_GetBuffer(dest, &buffer, PyBUF_SIMPLE|PyBUF_WRITABLE) != 0) {
        goto out;
    }
    buffer_acquired = 1;
    if (max_mutations == -1) {
        max_mutations = buffer.len / sizeof(double);
    } else if (max_mutations * sizeof(double) > buffer.len) {
        PyErr_SetString(PyExc_BufferError,
            "dest buffer is too small for the results");
        goto out;
    }

    Py_BEGIN_ALLOW_THREADS
    err = ld_calc_get_r2_array(
        self->ld_calc, (size_t) source_index, direction,
        (size_t) max_mutations, max_distance,
        (double *) buffer.buf, &num_r2_values);
    Py_END_ALLOW_THREADS
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) num_r2_values);
out:
    if (buffer_acquired) {
        PyBuffer_Release(&buffer);
    }
    return ret;
}

static PyMemberDef LdCalculator_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef LdCalculator_methods[] = {
    {"get_r2", (PyCFunction) LdCalculator_get_r2, METH_VARARGS,
        "Returns the value of the r2 statistic between the specified pair of "
        "mutation indexes"},
    {"get_r2_array", (PyCFunction) LdCalculator_get_r2_array,
        METH_VARARGS|METH_KEYWORDS,
        "Returns r2 statistic for a given mutation over specified range"},
    {NULL}  /* Sentinel */
};

static PyTypeObject LdCalculatorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.LdCalculator",             /* tp_name */
    sizeof(LdCalculator),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)LdCalculator_dealloc, /* tp_dealloc */
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
    "LdCalculator objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    LdCalculator_methods,             /* tp_methods */
    LdCalculator_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LdCalculator_init,      /* tp_init */
};


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
Simulator_parse_population_configuration(Simulator *self,
        PyObject *py_pop_config)
{
    int ret = -1;
    Py_ssize_t j, num_populations;
    double initial_size, growth_rate;
    int err;
    PyObject *item, *value;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    num_populations = PyList_Size(py_pop_config);
    if (num_populations == 0) {
        PyErr_SetString(PyExc_ValueError, "Empty population configuration");
        goto out;
    }
    err = msp_set_num_populations(self->sim, (size_t) num_populations);
    if (err != 0) {
        handle_input_error(err);
        goto out;
    }
    for (j = 0; j < PyList_Size(py_pop_config); j++) {
        item = PyList_GetItem(py_pop_config, j);
        if (!PyDict_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "not a dictionary");
            goto out;
        }
        value = get_dict_number(item, "initial_size");
        if (value == NULL) {
            goto out;
        }
        initial_size = PyFloat_AsDouble(value);
        value = get_dict_number(item, "growth_rate");
        if (value == NULL) {
            goto out;
        }
        growth_rate = PyFloat_AsDouble(value);
        err = msp_set_population_configuration(self->sim, j,
                initial_size, growth_rate);
        if (err != 0) {
            handle_input_error(err);
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

static int
Simulator_parse_migration_matrix(Simulator *self,
        PyObject *py_migration_matrix)
{
    int ret = -1;
    int err;
    size_t num_populations, j, size;
    PyObject *value;
    double *migration_matrix = NULL;

    size = PyList_Size(py_migration_matrix);
    migration_matrix = PyMem_Malloc(size * sizeof(double));
    if (migration_matrix == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    num_populations = msp_get_num_populations(self->sim);
    if (num_populations * num_populations != size) {
        PyErr_Format(PyExc_ValueError,
            "Migration matrix must be a flattened "
            "num_populations*num_populations square array");
        goto out;
    }
    for (j = 0; j < size; j++) {
        value = PyList_GetItem(py_migration_matrix, j);
        if (!PyNumber_Check(value)) {
            PyErr_Format(PyExc_TypeError, "Migration rate not a number");
            goto out;
        }
        migration_matrix[j] = PyFloat_AsDouble(value);
        if (migration_matrix[j] < 0.0) {
            PyErr_Format(PyExc_ValueError, "Negative values not permitted");
            goto out;
        }
    }
    err = msp_set_migration_matrix(self->sim, size, migration_matrix);
    if (err != 0) {
        handle_input_error(err);
        goto out;
    }
    ret = 0;
out:
    if (migration_matrix != NULL) {
        PyMem_Free(migration_matrix);
    }
    return ret;
}

static int
Simulator_parse_simulation_model(Simulator *self, PyObject *py_model)
{
    int ret = -1;
    int err = -1;
    PyObject *py_name = NULL;
    PyObject *hudson_s = NULL;
    PyObject *smc_s = NULL;
    PyObject *smc_prime_s = NULL;
    PyObject *dirac_s = NULL;
    PyObject *beta_s = NULL;
    PyObject *value;
    int is_hudson, is_smc, is_smc_prime, is_dirac, is_beta;
    double psi, alpha, truncation_point;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    hudson_s = Py_BuildValue("s", "hudson");
    if (hudson_s == NULL) {
        goto out;
    }
    smc_s = Py_BuildValue("s", "smc");
    if (smc_s == NULL) {
        goto out;
    }
    smc_prime_s = Py_BuildValue("s", "smc_prime");
    if (smc_prime_s == NULL) {
        goto out;
    }
    dirac_s = Py_BuildValue("s", "dirac");
    if (dirac_s == NULL) {
        goto out;
    }
    beta_s = Py_BuildValue("s", "beta");
    if (beta_s == NULL) {
        goto out;
    }
    py_name = get_dict_value(py_model, "name");
    if (py_name == NULL) {
        goto out;
    }

    /* We need to go through this tedious rigmarole because of string
     * handling in Python 3. By pushing the comparison up into Python
     * we don't need to worry about encodings, etc, etc.
     */
    is_hudson = PyObject_RichCompareBool(py_name, hudson_s, Py_EQ);
    if (is_hudson == -1) {
        goto out;
    }
    if (is_hudson) {
        err = msp_set_simulation_model_non_parametric(self->sim, MSP_MODEL_HUDSON);
    }

    is_smc = PyObject_RichCompareBool(py_name, smc_s, Py_EQ);
    if (is_smc == -1) {
        goto out;
    }
    if (is_smc) {
        err = msp_set_simulation_model_non_parametric(self->sim, MSP_MODEL_SMC);
    }

    is_smc_prime = PyObject_RichCompareBool(py_name, smc_prime_s, Py_EQ);
    if (is_smc_prime == -1) {
        goto out;
    }
    if (is_smc_prime) {
        err = msp_set_simulation_model_non_parametric(self->sim, MSP_MODEL_SMC_PRIME);
    }

    is_dirac = PyObject_RichCompareBool(py_name, dirac_s, Py_EQ);
    if (is_dirac == -1) {
        goto out;
    }
    if (is_dirac) {
        value = get_dict_number(py_model, "psi");
        if (value == NULL) {
            goto out;
        }
        psi = PyFloat_AsDouble(value);
        /* TODO range checking on psi */
        err = msp_set_simulation_model_dirac(self->sim, psi);
    }

    is_beta = PyObject_RichCompareBool(py_name, beta_s, Py_EQ);
    if (is_beta == -1) {
        goto out;
    }
    if (is_beta) {
        value = get_dict_number(py_model, "alpha");
        if (value == NULL) {
            goto out;
        }
        alpha = PyFloat_AsDouble(value);
        value = get_dict_number(py_model, "truncation_point");
        if (value == NULL) {
            goto out;
        }
        truncation_point = PyFloat_AsDouble(value);
        /* TODO range checking on alpha and truncation_point */
        err = msp_set_simulation_model_beta(self->sim, alpha, truncation_point);
    }

    if (! (is_hudson || is_smc || is_smc_prime || is_dirac || is_beta)) {
        PyErr_SetString(PyExc_ValueError, "Unknown simulation model");
        goto out;
    }
    if (err != 0) {
        handle_input_error(err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(hudson_s);
    Py_XDECREF(smc_s);
    Py_XDECREF(smc_prime_s);
    Py_XDECREF(beta_s);
    Py_XDECREF(dirac_s);
    return ret;
}


static int
Simulator_parse_demographic_events(Simulator *self, PyObject *py_events)
{
    int ret = -1;
    Py_ssize_t j;
    double time, initial_size, growth_rate, migration_rate, proportion,
           strength;
    int err, population_id, matrix_index, source, destination;
    int is_population_parameter_change, is_migration_rate_change,
        is_mass_migration, is_simple_bottleneck, is_instantaneous_bottleneck;
    PyObject *item, *value, *type;
    PyObject *population_parameter_change_s = NULL;
    PyObject *migration_rate_change_s = NULL;
    PyObject *mass_migration_s = NULL;
    PyObject *simple_bottleneck_s = NULL;
    PyObject *instantaneous_bottleneck_s = NULL;
    PyObject *initial_size_s = NULL;
    PyObject *growth_rate_s = NULL;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    /* Create the Python strings for comparison with the types and
     * dictionary lookups */
    population_parameter_change_s = Py_BuildValue("s",
            "population_parameters_change");
    if (population_parameter_change_s == NULL) {
        goto out;
    }
    migration_rate_change_s = Py_BuildValue("s", "migration_rate_change");
    if (migration_rate_change_s == NULL) {
        goto out;
    }
    mass_migration_s = Py_BuildValue("s", "mass_migration");
    if (mass_migration_s == NULL) {
        goto out;
    }
    simple_bottleneck_s = Py_BuildValue("s", "simple_bottleneck");
    if (simple_bottleneck_s == NULL) {
        goto out;
    }
    instantaneous_bottleneck_s = Py_BuildValue("s", "instantaneous_bottleneck");
    if (instantaneous_bottleneck_s == NULL) {
        goto out;
    }
    initial_size_s = Py_BuildValue("s", "initial_size");
    if (initial_size_s == NULL) {
        goto out;
    }
    growth_rate_s = Py_BuildValue("s", "growth_rate");
    if (growth_rate_s == NULL) {
        goto out;
    }

    for (j = 0; j < PyList_Size(py_events); j++) {
        item = PyList_GetItem(py_events, j);
        if (!PyDict_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "not a dictionary");
            goto out;
        }
        value = get_dict_number(item, "time");
        if (value == NULL) {
            goto out;
        }
        time = PyFloat_AsDouble(value);
        if (time < 0) {
            PyErr_SetString(PyExc_ValueError, "negative times not valid");
            goto out;
        }
        type = get_dict_value(item, "type");
        if (type == NULL) {
            goto out;
        }
        /* We need to go through this tedious rigmarole because of string
         * handling in Python 3. By pushing the comparison up into Python
         * we don't need to worry about encodings, etc, etc.
         */
        is_population_parameter_change = PyObject_RichCompareBool(type,
                population_parameter_change_s, Py_EQ);
        if (is_population_parameter_change == -1) {
            goto out;
        }
        is_migration_rate_change = PyObject_RichCompareBool(type,
                migration_rate_change_s, Py_EQ);
        if (is_migration_rate_change == -1) {
            goto out;
        }
        is_mass_migration = PyObject_RichCompareBool(type, mass_migration_s,
                Py_EQ);
        if (is_mass_migration == -1) {
            goto out;
        }
        is_simple_bottleneck = PyObject_RichCompareBool(
                type, simple_bottleneck_s, Py_EQ);
        if (is_simple_bottleneck == -1) {
            goto out;
        }
        is_instantaneous_bottleneck = PyObject_RichCompareBool(
                type, instantaneous_bottleneck_s, Py_EQ);
        if (is_instantaneous_bottleneck == -1) {
            goto out;
        }
        if (is_population_parameter_change) {
            initial_size = GSL_NAN;
            if (PyDict_Contains(item, initial_size_s)) {
                value = get_dict_number(item, "initial_size");
                if (value == NULL) {
                    goto out;
                }
                initial_size = PyFloat_AsDouble(value);
            }
            growth_rate = GSL_NAN;
            if (PyDict_Contains(item, growth_rate_s)) {
                value = get_dict_number(item, "growth_rate");
                if (value == NULL) {
                    goto out;
                }
                growth_rate = PyFloat_AsDouble(value);
            }
            value = get_dict_number(item, "population_id");
            if (value == NULL) {
                goto out;
            }
            population_id = (int) PyLong_AsLong(value);
            err = msp_add_population_parameters_change(self->sim, time,
                    population_id, initial_size, growth_rate);
        } else if (is_migration_rate_change) {
            value = get_dict_number(item, "migration_rate");
            if (value == NULL) {
                goto out;
            }
            migration_rate = PyFloat_AsDouble(value);
            value = get_dict_number(item, "matrix_index");
            if (value == NULL) {
                goto out;
            }
            matrix_index = (int) PyLong_AsLong(value);
            err = msp_add_migration_rate_change(self->sim, time, matrix_index,
                    migration_rate);
        } else if (is_mass_migration) {
            value = get_dict_number(item, "proportion");
            if (value == NULL) {
                goto out;
            }
            proportion = PyFloat_AsDouble(value);
            value = get_dict_number(item, "source");
            if (value == NULL) {
                goto out;
            }
            source = (int) PyLong_AsLong(value);
            value = get_dict_number(item, "destination");
            if (value == NULL) {
                goto out;
            }
            destination = (int) PyLong_AsLong(value);
            err = msp_add_mass_migration(self->sim, time, source, destination,
                    proportion);
        } else if (is_simple_bottleneck) {
            value = get_dict_number(item, "proportion");
            if (value == NULL) {
                goto out;
            }
            proportion = PyFloat_AsDouble(value);
            value = get_dict_number(item, "population_id");
            if (value == NULL) {
                goto out;
            }
            population_id = (int) PyLong_AsLong(value);
            err = msp_add_simple_bottleneck(self->sim, time, population_id,
                    proportion);
        } else if (is_instantaneous_bottleneck) {
            value = get_dict_number(item, "strength");
            if (value == NULL) {
                goto out;
            }
            strength = PyFloat_AsDouble(value);
            value = get_dict_number(item, "population_id");
            if (value == NULL) {
                goto out;
            }
            population_id = (int) PyLong_AsLong(value);
            err = msp_add_instantaneous_bottleneck(self->sim, time, population_id,
                    strength);
        } else {
            PyErr_Format(PyExc_ValueError, "Unknown demographic event type");
            goto out;
        }
        if (err != 0) {
            handle_input_error(err);
            goto out;
        }
    }
    ret = 0;
out:
    Py_XDECREF(population_parameter_change_s);
    Py_XDECREF(migration_rate_change_s);
    Py_XDECREF(mass_migration_s);
    Py_XDECREF(simple_bottleneck_s);
    Py_XDECREF(instantaneous_bottleneck_s);
    Py_XDECREF(initial_size_s);
    Py_XDECREF(growth_rate_s);
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
    Py_XDECREF(self->random_generator);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
Simulator_init(Simulator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int sim_ret;
    static char *kwlist[] = {"samples", "random_generator",
        "num_loci", "scaled_recombination_rate",
        "population_configuration", "migration_matrix", "demographic_events",
        "model", "max_memory", "avl_node_block_size", "segment_block_size",
        "node_mapping_block_size", "coalescence_record_block_size",
        "migration_block_size", "store_migrations", NULL};
    PyObject *py_samples = NULL;
    PyObject *migration_matrix = NULL;
    PyObject *population_configuration = NULL;
    PyObject *demographic_events = NULL;
    PyObject *py_model = NULL;
    RandomGenerator *random_generator = NULL;
    sample_t *samples = NULL;
    /* parameter defaults */
    Py_ssize_t sample_size = 2;
    Py_ssize_t num_loci = 1;
    double scaled_recombination_rate = 0.0;
    Py_ssize_t max_memory = 10 * 1024 * 1024;
    Py_ssize_t avl_node_block_size = 10;
    Py_ssize_t segment_block_size = 10;
    Py_ssize_t node_mapping_block_size = 10;
    Py_ssize_t coalescence_record_block_size = 10;
    Py_ssize_t migration_block_size = 10;
    int store_migrations = 0;

    self->sim = NULL;
    self->random_generator = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|ndO!O!O!O!nnnnnni", kwlist,
            &PyList_Type, &py_samples,
            &RandomGeneratorType, &random_generator,
            &num_loci, &scaled_recombination_rate,
            &PyList_Type, &population_configuration,
            &PyList_Type, &migration_matrix,
            &PyList_Type, &demographic_events,
            &PyDict_Type, &py_model,
            &max_memory, &avl_node_block_size, &segment_block_size,
            &node_mapping_block_size, &coalescence_record_block_size,
            &migration_block_size, &store_migrations)) {
        goto out;
    }
    self->random_generator = random_generator;
    Py_INCREF(self->random_generator);
    if (RandomGenerator_check_state(self->random_generator) != 0) {
        goto out;
    }
    if (parse_samples(py_samples, &sample_size, &samples) != 0) {
        goto out;
    }
    self->sim = PyMem_Malloc(sizeof(msp_t));
    if (self->sim == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    sim_ret = msp_alloc(self->sim, (size_t) sample_size, samples,
            self->random_generator->rng);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    if (py_model != NULL) {
        if (Simulator_parse_simulation_model(self, py_model) != 0) {
            goto out;
        }
    }
    sim_ret = msp_set_store_migrations(self->sim, (bool) store_migrations);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    sim_ret = msp_set_num_loci(self->sim, (size_t) num_loci);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    sim_ret = msp_set_scaled_recombination_rate(self->sim,
            scaled_recombination_rate);
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
    sim_ret = msp_set_migration_block_size(self->sim,
            (size_t) migration_block_size);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    if (population_configuration != NULL) {
        if (Simulator_parse_population_configuration(self,
                population_configuration) != 0) {
            goto out;
        }
        if (migration_matrix == NULL) {
            PyErr_SetString(PyExc_ValueError,
                "A migration matrix must be provided when a non-default "
                "population configuration is used.");
            goto out;
        }
        if (Simulator_parse_migration_matrix(self,
                migration_matrix) != 0) {
            goto out;
        }
    } else if (migration_matrix != NULL) {
        PyErr_SetString(PyExc_ValueError,
            "Cannot supply migration_matrix without "
            "population_configuration.");
        goto out;
    }
    if (demographic_events != NULL) {
        if (Simulator_parse_demographic_events(self, demographic_events) != 0) {
            goto out;
        }
    }
    sim_ret = msp_initialise(self->sim);
    if (sim_ret != 0) {
        handle_input_error(sim_ret);
        goto out;
    }
    ret = 0;
out:
    if (samples != NULL) {
        PyMem_Free(samples);
    }
    return ret;
}

static PyMemberDef Simulator_members[] = {
    {NULL}  /* Sentinel */
};


static PyObject *
Simulator_get_model(Simulator *self)
{
    PyObject *ret = NULL;
    PyObject *d = NULL;
    PyObject *value = NULL;
    simulation_model_t *model;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    model = msp_get_model(self->sim);
    d = Py_BuildValue("{ss}", "name", msp_get_model_name(self->sim));
    if (model->type == MSP_MODEL_DIRAC) {
        value = Py_BuildValue("d", model->params.dirac_coalescent.psi);
        if (value == NULL) {
            goto out;
        }
        if (PyDict_SetItemString(d, "psi", value) != 0) {
            goto out;
        }
        Py_DECREF(value);
        value = NULL;
    } else if (model->type == MSP_MODEL_BETA) {
        value = Py_BuildValue("d", model->params.beta_coalescent.alpha);
        if (value == NULL) {
            goto out;
        }
        if (PyDict_SetItemString(d, "alpha", value) != 0) {
            goto out;
        }
        Py_DECREF(value);
        value = NULL;
        value = Py_BuildValue("d", model->params.beta_coalescent.truncation_point);
        if (value == NULL) {
            goto out;
        }
        if (PyDict_SetItemString(d, "truncation_point", value) != 0) {
            goto out;
        }
        Py_DECREF(value);
        value = NULL;
    }
    ret = d;
    d = NULL;
out:
    Py_XDECREF(d);
    Py_XDECREF(value);
    return ret;
}

static PyObject *
Simulator_get_num_loci(Simulator *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_num_loci(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_store_migrations(Simulator *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("i", (Py_ssize_t) msp_get_store_migrations(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_sample_size(Simulator *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_sample_size(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_populations(Simulator *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_num_populations(self->sim));
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
Simulator_get_max_memory(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", self->sim->max_memory);
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
Simulator_get_migration_block_size(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n",
            (Py_ssize_t) self->sim->migration_block_size);
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
Simulator_get_num_common_ancestor_events(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n",
        (Py_ssize_t) msp_get_num_common_ancestor_events(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_rejected_common_ancestor_events(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n",
        (Py_ssize_t) msp_get_num_rejected_common_ancestor_events(self->sim));
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
    ret = Py_BuildValue("n",
        (Py_ssize_t) msp_get_num_recombination_events(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_migration_events(Simulator  *self)
{
    PyObject *ret = NULL;
    size_t *num_migration_events = NULL;
    size_t num_populations;
    int err;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    num_populations = msp_get_num_populations(self->sim);
    num_migration_events = PyMem_Malloc(
        num_populations * num_populations * sizeof(size_t));
    if (num_migration_events == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = msp_get_num_migration_events(self->sim, num_migration_events);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = convert_integer_list(num_migration_events,
            num_populations * num_populations);
out:
    if (num_migration_events != NULL) {
        PyMem_Free(num_migration_events);
    }
    return ret;
}

static PyObject *
Simulator_get_num_multiple_recombination_events(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sim->num_multiple_re_events);
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
Simulator_get_num_migration_blocks(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n",
            (Py_ssize_t) msp_get_num_migration_blocks(self->sim));
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
Simulator_get_num_migrations(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n",
            (Py_ssize_t) msp_get_num_migrations(self->sim));
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
        t = Py_BuildValue("(I,I,I,I)", u->left, u->right, u->value,
                u->population_id);
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
Simulator_get_breakpoints(Simulator *self)
{
    PyObject *ret = NULL;
    size_t *breakpoints = NULL;
    size_t num_breakpoints;
    int err;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    num_breakpoints = msp_get_num_breakpoints(self->sim);
    breakpoints = PyMem_Malloc(num_breakpoints * sizeof(size_t));
    if (breakpoints == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = msp_get_breakpoints(self->sim, breakpoints);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = convert_integer_list(breakpoints, num_breakpoints);
out:
    if (breakpoints != NULL) {
        PyMem_Free(breakpoints);
    }
    return ret;
}

static PyObject *
Simulator_get_migration_matrix(Simulator *self)
{
    PyObject *ret = NULL;
    double *migration_matrix = NULL;
    size_t N;
    int err;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    N = msp_get_num_populations(self->sim);
    migration_matrix = PyMem_Malloc(N * N * sizeof(double));
    if (migration_matrix == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = msp_get_migration_matrix(self->sim, migration_matrix);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = convert_float_list(migration_matrix, N * N);
out:
    if (migration_matrix != NULL) {
        PyMem_Free(migration_matrix);
    }
    return ret;
}

static PyObject *
Simulator_get_coalescence_records(Simulator *self)
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
    err = msp_get_coalescence_records(self->sim, &coalescence_records);
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
        py_cr = make_coalescence_record(cr);
        if (py_cr == NULL) {
            Py_DECREF(l);
            goto out;
        }
        PyList_SET_ITEM(l, j, py_cr);
    }
    ret = l;
out:
    return ret;
}

static PyObject *
Simulator_get_migrations(Simulator *self)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_mr = NULL;
    migration_t *migrations = NULL;
    migration_t *mr;
    size_t num_migrations, j;
    int err;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    num_migrations = msp_get_num_migrations(self->sim);
    err = msp_get_migrations(self->sim, &migrations);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    l = PyList_New(num_migrations);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < num_migrations; j++) {
        mr = &migrations[j];
        py_mr = make_migration(mr);
        if (py_mr == NULL) {
            Py_DECREF(l);
            goto out;
        }
        PyList_SET_ITEM(l, j, py_mr);
    }
    ret = l;
out:
    return ret;
}

static PyObject *
Simulator_get_population_configuration(Simulator *self)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *d = NULL;
    size_t j = 0;
    size_t num_populations;
    int sim_ret = 0;
    double initial_size, growth_rate;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    num_populations = msp_get_num_populations(self->sim);
    l = PyList_New(num_populations);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < num_populations; j++) {
        sim_ret = msp_get_population_configuration(self->sim, j,
            &initial_size, &growth_rate);
        if (sim_ret != 0) {
            handle_library_error(sim_ret);
            goto out;
        }
        d = Py_BuildValue("{s:d,s:d}",
               "initial_size", initial_size,
               "growth_rate", growth_rate);
        if (d == NULL) {
            goto out;
        }
        PyList_SET_ITEM(l, j, d);
    }
    ret = l;
    l = NULL;
out:
    Py_XDECREF(l);
    return ret;
}

static PyObject *
Simulator_get_samples(Simulator *self)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *t = NULL;
    sample_t *samples = NULL;
    size_t j = 0;
    size_t sample_size;
    int population;
    int sim_ret = 0;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    sample_size = msp_get_sample_size(self->sim);
    sim_ret = msp_get_samples(self->sim, &samples);
    if (sim_ret != 0) {
        handle_library_error(sim_ret);
        goto out;
    }
    l = PyList_New(sample_size);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < sample_size; j++) {
        population = samples[j].population_id == MSP_NULL_POPULATION_ID? -1:
            samples[j].population_id;
        t = Py_BuildValue("id", population, samples[j].time);
        if (t == NULL) {
            goto out;
        }
        PyList_SET_ITEM(l, j, t);
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
    /* return True if complete coalescence has occured */
    ret = coalesced ? Py_True : Py_False;
    Py_INCREF(ret);
out:
    return ret;
}

static PyObject *
Simulator_run_event(Simulator *self)
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

static PyObject *
Simulator_populate_tables(Simulator *self, PyObject *args, PyObject *kwds)
{
    int err;
    PyObject *ret = NULL;
    NodeTable *nodes = NULL;
    EdgesetTable *edgesets = NULL;
    MigrationTable *migrations = NULL;
    RecombinationMap *recombination_map = NULL;
    recomb_map_t *recomb_map = NULL;
    double Ne = 0.25; /* default to coalescent time */
    static char *kwlist[] = {"nodes", "edgesets", "migrations",
        "Ne", "recombination_map", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!O!|dO!", kwlist,
            &NodeTableType, &nodes,
            &EdgesetTableType, &edgesets,
            &MigrationTableType, &migrations,
            &Ne,
            &RecombinationMapType, &recombination_map)) {
        goto out;
    }
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (NodeTable_check_state(nodes) != 0) {
        goto out;
    }
    if (EdgesetTable_check_state(edgesets) != 0) {
        goto out;
    }
    if (MigrationTable_check_state(migrations) != 0) {
        goto out;
    }
    if (recombination_map != NULL) {
        if (RecombinationMap_check_recomb_map(recombination_map) != 0) {
            goto out;
        }
        recomb_map = recombination_map->recomb_map;
    }
    err = msp_populate_tables(self->sim, Ne, recomb_map,
        nodes->node_table, edgesets->edgeset_table,
        migrations->migration_table);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
Simulator_reset(Simulator *self)
{
    PyObject *ret = NULL;
    int status;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    status = msp_reset(self->sim);
    if (status < 0) {
        handle_library_error(status);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
Simulator_debug_demography(Simulator *self)
{
    PyObject *ret = NULL;
    int status;
    double end_time;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    status = msp_debug_demography(self->sim, &end_time);
    if (status < 0) {
        handle_library_error(status);
        goto out;
    }
    ret = Py_BuildValue("d", end_time);
out:
    return ret;
}


static PyMethodDef Simulator_methods[] = {
    {"get_model", (PyCFunction) Simulator_get_model, METH_NOARGS,
            "Returns the simulation model" },
    {"get_num_loci", (PyCFunction) Simulator_get_num_loci, METH_NOARGS,
            "Returns the number of loci" },
    {"get_store_migrations",
            (PyCFunction) Simulator_get_store_migrations, METH_NOARGS,
            "Returns True if the simulator should store migration records." },
    {"get_sample_size", (PyCFunction) Simulator_get_sample_size, METH_NOARGS,
            "Returns the sample size" },
    {"get_num_populations", (PyCFunction) Simulator_get_num_populations, METH_NOARGS,
            "Returns the number of populations." },
    {"get_scaled_recombination_rate",
            (PyCFunction) Simulator_get_scaled_recombination_rate, METH_NOARGS,
            "Returns the scaled recombination rate." },
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
    {"get_migration_block_size",
            (PyCFunction) Simulator_get_migration_block_size,
            METH_NOARGS, "Returns the migration record block size" },
    {"get_time", (PyCFunction) Simulator_get_time, METH_NOARGS,
            "Returns the current simulation time" },
    {"get_num_ancestors", (PyCFunction) Simulator_get_num_ancestors, METH_NOARGS,
            "Returns the number of ancestors" },
    {"get_num_common_ancestor_events",
            (PyCFunction) Simulator_get_num_common_ancestor_events, METH_NOARGS,
            "Returns the number of common_ancestor_events" },
    {"get_num_rejected_common_ancestor_events",
            (PyCFunction) Simulator_get_num_rejected_common_ancestor_events,
            METH_NOARGS, "Returns the number of rejected common_ancestor_events" },
    {"get_num_recombination_events",
            (PyCFunction) Simulator_get_num_recombination_events, METH_NOARGS,
            "Returns the number of recombination_events" },
    {"get_num_migration_events",
            (PyCFunction) Simulator_get_num_migration_events, METH_NOARGS,
            "Returns the number of migration events" },
    {"get_num_multiple_recombination_events",
            (PyCFunction) Simulator_get_num_multiple_recombination_events,
            METH_NOARGS,
            "Returns the number of recombination_events that occur at an "
            "existing breakpoint" },
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
    {"get_num_migration_blocks",
            (PyCFunction) Simulator_get_num_migration_blocks, METH_NOARGS,
            "Returns the number of coalescence record memory blocks"},
    {"get_num_breakpoints", (PyCFunction) Simulator_get_num_breakpoints,
            METH_NOARGS, "Returns the number of recombination breakpoints" },
    {"get_num_coalescence_records",
            (PyCFunction) Simulator_get_num_coalescence_records,
            METH_NOARGS, "Returns the number of coalescence records" },
    {"get_num_migrations",
            (PyCFunction) Simulator_get_num_migrations,
            METH_NOARGS, "Returns the number of migration records" },
    {"get_used_memory", (PyCFunction) Simulator_get_used_memory,
            METH_NOARGS, "Returns the approximate amount of memory used." },
    {"get_ancestors", (PyCFunction) Simulator_get_ancestors, METH_NOARGS,
            "Returns the ancestors" },
    {"get_breakpoints", (PyCFunction) Simulator_get_breakpoints,
            METH_NOARGS, "Returns the list of breakpoints." },
    {"get_migration_matrix", (PyCFunction) Simulator_get_migration_matrix,
            METH_NOARGS, "Returns the migration matrix." },
    {"get_coalescence_records", (PyCFunction) Simulator_get_coalescence_records,
            METH_NOARGS, "Returns the coalescence records." },
    {"get_migrations", (PyCFunction) Simulator_get_migrations,
            METH_NOARGS, "Returns the migration records." },
    {"get_population_configuration",
            (PyCFunction) Simulator_get_population_configuration, METH_NOARGS,
            "Returns the population configurations"},
    {"get_samples",
            (PyCFunction) Simulator_get_samples, METH_NOARGS,
            "Returns the samples"},
    {"run", (PyCFunction) Simulator_run, METH_VARARGS,
            "Simulates until at most the specified time. Returns True\
            if sample has coalesced and False otherwise." },
    {"reset", (PyCFunction) Simulator_reset, METH_NOARGS,
            "Resets the simulation so it's ready for another replicate."},
    {"populate_tables",
        (PyCFunction) Simulator_populate_tables, METH_VARARGS|METH_KEYWORDS,
        "Updates the specified tables to reflect the state of this simulator"},
    {"run_event", (PyCFunction) Simulator_run_event, METH_NOARGS,
            "Simulates exactly one event. Returns True "
            "if sample has coalesced and False otherwise." },
    {"debug_demography", (PyCFunction) Simulator_debug_demography, METH_NOARGS,
            "Runs the state of the simulator forward for one demographic event."},
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


static PyObject *
msprime_h5close(PyObject *self)
{
    herr_t status;
    PyObject *ret = NULL;

    status = H5close();
    if (status != 0) {
        PyErr_SetString(PyExc_SystemError, "Error calling H5close");
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
msprime_get_library_version_str(PyObject *self)
{
    return Py_BuildValue("s", MSP_LIBRARY_VERSION_STR);
}


static PyMethodDef msprime_methods[] = {
    {"get_gsl_version", (PyCFunction) msprime_get_gsl_version, METH_NOARGS,
            "Returns the version of GSL we are linking against." },
    {"get_hdf5_version", (PyCFunction) msprime_get_hdf5_version, METH_NOARGS,
            "Returns the version of HDF5 we are linking against." },
    {"h5close", (PyCFunction) msprime_h5close, METH_NOARGS,
            "Calls H5close()" },
    {"get_library_version_str", (PyCFunction) msprime_get_library_version_str,
            METH_NOARGS, "Returns the version of the msp C library." },
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
    import_array();

    /* RandomGenerator type */
    RandomGeneratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&RandomGeneratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&RandomGeneratorType);
    PyModule_AddObject(module, "RandomGenerator", (PyObject *) &RandomGeneratorType);

    /* NodeTable type */
    NodeTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&NodeTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&NodeTableType);
    PyModule_AddObject(module, "NodeTable", (PyObject *) &NodeTableType);

    /* NodeTable type */
    NodeTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&NodeTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&NodeTableType);
    PyModule_AddObject(module, "NodeTable", (PyObject *) &NodeTableType);

    /* EdgesetTable type */
    EdgesetTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&EdgesetTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&EdgesetTableType);
    PyModule_AddObject(module, "EdgesetTable", (PyObject *) &EdgesetTableType);

    /* MigrationTable type */
    MigrationTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&MigrationTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&MigrationTableType);
    PyModule_AddObject(module, "MigrationTable", (PyObject *) &MigrationTableType);

    /* MutationTypeTable type */
    MutationTypeTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&MutationTypeTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&MutationTypeTableType);
    PyModule_AddObject(module, "MutationTypeTable", (PyObject *) &MutationTypeTableType);

    /* MutationTable type */
    MutationTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&MutationTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&MutationTableType);
    PyModule_AddObject(module, "MutationTable", (PyObject *) &MutationTableType);

    /* MutationGenerator type */
    MutationGeneratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&MutationGeneratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&MutationGeneratorType);
    PyModule_AddObject(module, "MutationGenerator", (PyObject *) &MutationGeneratorType);

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

    /* RecombinationMap type */
    RecombinationMapType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&RecombinationMapType) < 0) {
        INITERROR;
    }
    Py_INCREF(&RecombinationMapType);
    PyModule_AddObject(module, "RecombinationMap", (PyObject *) &RecombinationMapType);

    /* SparseTree type */
    SparseTreeType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&SparseTreeType) < 0) {
        INITERROR;
    }
    Py_INCREF(&SparseTreeType);
    PyModule_AddObject(module, "SparseTree", (PyObject *) &SparseTreeType);

    /* SparseTreeIterator type */
    SparseTreeIteratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&SparseTreeIteratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&SparseTreeIteratorType);
    PyModule_AddObject(module, "SparseTreeIterator",
            (PyObject *) &SparseTreeIteratorType);

    /* TreeDiffIterator type */
    TreeDiffIteratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&TreeDiffIteratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&TreeDiffIteratorType);
    PyModule_AddObject(module, "TreeDiffIterator", (PyObject *) &TreeDiffIteratorType);

    /* LeafListIterator type */
    LeafListIteratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&LeafListIteratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&LeafListIteratorType);
    PyModule_AddObject(module, "LeafListIterator", (PyObject *) &LeafListIteratorType);

    /* NewickConverter type */
    NewickConverterType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&NewickConverterType) < 0) {
        INITERROR;
    }
    Py_INCREF(&NewickConverterType);
    PyModule_AddObject(module, "NewickConverter", (PyObject *) &NewickConverterType);

    /* VcfConverter type */
    VcfConverterType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&VcfConverterType) < 0) {
        INITERROR;
    }
    Py_INCREF(&VcfConverterType);
    PyModule_AddObject(module, "VcfConverter", (PyObject *) &VcfConverterType);

    /* HaplotypeGenerator type */
    HaplotypeGeneratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&HaplotypeGeneratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&HaplotypeGeneratorType);
    PyModule_AddObject(module, "HaplotypeGenerator",
            (PyObject *) &HaplotypeGeneratorType);

    /* VariantGenerator type */
    VariantGeneratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&VariantGeneratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&VariantGeneratorType);
    PyModule_AddObject(module, "VariantGenerator", (PyObject *) &VariantGeneratorType);

    /* LdCalculator type */
    LdCalculatorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&LdCalculatorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&LdCalculatorType);
    PyModule_AddObject(module, "LdCalculator", (PyObject *) &LdCalculatorType);

    /* Errors and constants */
    MsprimeInputError = PyErr_NewException("_msprime.InputError", NULL, NULL);
    Py_INCREF(MsprimeInputError);
    PyModule_AddObject(module, "InputError", MsprimeInputError);
    MsprimeLibraryError = PyErr_NewException("_msprime.LibraryError", NULL, NULL);
    Py_INCREF(MsprimeLibraryError);
    PyModule_AddObject(module, "LibraryError", MsprimeLibraryError);

    /* Node flags */
    PyModule_AddIntConstant(module, "NODE_IS_SAMPLE", MSP_NODE_IS_SAMPLE);

    /* Tree flags */
    PyModule_AddIntConstant(module, "LEAF_COUNTS", MSP_LEAF_COUNTS);
    PyModule_AddIntConstant(module, "LEAF_LISTS", MSP_LEAF_LISTS);
    /* Directions */
    PyModule_AddIntConstant(module, "FORWARD", MSP_DIR_FORWARD);
    PyModule_AddIntConstant(module, "REVERSE", MSP_DIR_REVERSE);

    /* turn off GSL error handler so we don't abort on memory error */
    gsl_set_error_handler_off();

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}


