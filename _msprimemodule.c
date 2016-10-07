/*
** Copyright (C) 2014-2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
parse_coalescence_record(PyObject *tuple, coalescence_record_t *cr)
{
    int ret = -1;
    size_t size, j;
    long v;
    PyObject *item;
    PyObject *children;
    const char *err = "Coalescence records must be tuples of the form "
        "(left, right, node, (children), time, population_id)";

    if (!PyTuple_Check(tuple)) {
        PyErr_SetString(
            PyExc_TypeError, "Coalescence records must be a tuple.");
        goto out;
    }
    size = PyTuple_Size(tuple);
    if (size != 6) {
        PyErr_SetString(PyExc_ValueError, err);
        goto out;
    }
    /* left */
    item = PyTuple_GetItem(tuple, 0);
    if (!PyNumber_Check(item)) {
        PyErr_Format(PyExc_TypeError, "'left' is not number");
        goto out;
    }
    cr->left = PyFloat_AsDouble(item);
    if (PyErr_Occurred() != NULL) {
        goto out;
    }
    /* right */
    item = PyTuple_GetItem(tuple, 1);
    if (!PyNumber_Check(item)) {
        PyErr_Format(PyExc_TypeError, "'right' is not number");
        goto out;
    }
    cr->right = PyFloat_AsDouble(item);
    if (PyErr_Occurred() != NULL) {
        goto out;
    }
    /* node */
    item = PyTuple_GetItem(tuple, 2);
    if (!PyNumber_Check(item)) {
        PyErr_Format(PyExc_TypeError, "'node' is not number");
        goto out;
    }
    cr->node = (uint32_t) PyLong_AsLong(item);
    if (PyErr_Occurred() != NULL) {
        goto out;
    }
    /* children */
    children = PyTuple_GetItem(tuple, 3);
    if (!PyTuple_Check(children)) {
        PyErr_SetString(PyExc_TypeError, "children must be a tuple.");
        goto out;
    }
    if (PyTuple_Size(children) < 2) {
        PyErr_SetString(PyExc_ValueError, "Must be >= 2 children");
        goto out;
    }
    cr->num_children = PyTuple_Size(children);
    cr->children = PyMem_Malloc(cr->num_children * sizeof(uint32_t));
    if (cr->children == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    for (j = 0; j < cr->num_children; j++) {
        item = PyTuple_GetItem(children, j);
        if (!PyNumber_Check(item)) {
            PyErr_Format(PyExc_TypeError, "children[%d]' is not number",
                (int)j);
            goto out;
        }
        cr->children[j] = (uint32_t) PyLong_AsLong(item);
        if (PyErr_Occurred() != NULL) {
            goto out;
        }
    }
    /* time */
    item = PyTuple_GetItem(tuple, 4);
    if (!PyNumber_Check(item)) {
        PyErr_Format(PyExc_TypeError, "'time' is not number");
        goto out;
    }
    cr->time = PyFloat_AsDouble(item);
    if (PyErr_Occurred() != NULL) {
        goto out;
    }
    /* population */
    item = PyTuple_GetItem(tuple, 5);
    if (!PyNumber_Check(item)) {
        PyErr_Format(PyExc_TypeError, "'population_id' is not number");
        goto out;
    }
    v = PyLong_AsLong(item);
    if (PyErr_Occurred() != NULL) {
        goto out;
    }
    if (v == -1) {
        cr->population_id = MSP_NULL_POPULATION_ID;
    } else {
        cr->population_id = (uint8_t) v;
    }
    ret = 0;
out:
    return ret;
}

static int
parse_samples(PyObject *py_samples, Py_ssize_t *sample_size,
        sample_t **samples)
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
        ret_samples[j].population_id = (uint8_t) tmp_long;
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
        py_str = Py_BuildValue("s", list[j]);
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
convert_children(uint32_t *children, uint32_t num_children)
{
    PyObject *ret = NULL;
    PyObject *t;
    PyObject *py_int;
    uint32_t j;

    t = PyTuple_New(num_children);
    if (t == NULL) {
        goto out;
    }
    for (j = 0; j < num_children; j++) {
        py_int = Py_BuildValue("I", (unsigned int) children[j]);
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

static inline PyObject *
convert_mutation(mutation_t *mutation)
{
    return Py_BuildValue("dIn", mutation->position,
        (unsigned int) mutation->node, (Py_ssize_t) mutation->index);
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
        py_mutation = convert_mutation(&mutations[j]);
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

static PyObject *
make_coalescence_record(coalescence_record_t *cr)
{
    int population_id =
        cr->population_id == MSP_NULL_POPULATION_ID ? -1: cr->population_id;
    PyObject *children = NULL;
    PyObject *ret = NULL;

    children = convert_children(cr->children, cr->num_children);
    if (children == NULL) {
        goto out;
    }
    ret = Py_BuildValue("ddIOdi",
            cr->left, cr->right, (unsigned int) cr->node, children,
            cr->time, population_id);
out:
    Py_XDECREF(children);
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
Simulator_parse_demographic_events(Simulator *self, PyObject *py_events)
{
    int ret = -1;
    Py_ssize_t j;
    double time, initial_size, growth_rate, migration_rate, proportion;
    int err, population_id, matrix_index, source, destination;
    int is_population_parameter_change, is_migration_rate_change,
        is_mass_migration, is_bottleneck;
    PyObject *item, *value, *type;
    PyObject *population_parameter_change_s = NULL;
    PyObject *migration_rate_change_s = NULL;
    PyObject *mass_migration_s = NULL;
    PyObject *bottleneck_s = NULL;
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
    bottleneck_s = Py_BuildValue("s", "bottleneck");
    if (bottleneck_s == NULL) {
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
        is_bottleneck = PyObject_RichCompareBool(type, bottleneck_s, Py_EQ);
        if (is_bottleneck == -1) {
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
        } else if (is_bottleneck) {
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
            err = msp_add_bottleneck(self->sim, time, population_id,
                    proportion);
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
    Py_DECREF(population_parameter_change_s);
    Py_DECREF(migration_rate_change_s);
    Py_DECREF(mass_migration_s);
    Py_DECREF(bottleneck_s);
    Py_DECREF(initial_size_s);
    Py_DECREF(growth_rate_s);
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
        "max_memory", "avl_node_block_size", "segment_block_size",
        "node_mapping_block_size", "coalescence_record_block_size", NULL};
    PyObject *py_samples = NULL;
    PyObject *migration_matrix = NULL;
    PyObject *population_configuration = NULL;
    PyObject *demographic_events = NULL;
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

    self->sim = NULL;
    self->random_generator = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|ndO!O!O!nnnnn", kwlist,
            &PyList_Type, &py_samples,
            &RandomGeneratorType, &random_generator,
            &num_loci, &scaled_recombination_rate,
            &PyList_Type, &population_configuration,
            &PyList_Type, &migration_matrix,
            &PyList_Type, &demographic_events,
            &max_memory, &avl_node_block_size, &segment_block_size,
            &node_mapping_block_size, &coalescence_record_block_size)) {
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
        if (Simulator_parse_demographic_events(self,
                    demographic_events) != 0) {
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
Simulator_get_num_loci(Simulator  *self)
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
    {"get_num_loci", (PyCFunction) Simulator_get_num_loci, METH_NOARGS,
            "Returns the number of loci" },
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
    {"get_time", (PyCFunction) Simulator_get_time, METH_NOARGS,
            "Returns the current simulation time" },
    {"get_num_ancestors", (PyCFunction) Simulator_get_num_ancestors, METH_NOARGS,
            "Returns the number of ancestors" },
    {"get_num_common_ancestor_events",
            (PyCFunction) Simulator_get_num_common_ancestor_events, METH_NOARGS,
            "Returns the number of common_ancestor_events" },
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
    {"get_migration_matrix", (PyCFunction) Simulator_get_migration_matrix,
            METH_NOARGS, "Returns the migration matrix." },
    {"get_coalescence_records", (PyCFunction) Simulator_get_coalescence_records,
            METH_NOARGS, "Returns the coalescence records." },
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
    {"run_event", (PyCFunction) Simulator_run_event, METH_NOARGS,
            "Simulates exactly one event. Returns True\
            if sample has coalesced and False otherwise." },
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
    self->tree_sequence = NULL;
    return 0;
}

static PyObject *
TreeSequence_create(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Simulator *sim = NULL;
    RecombinationMap *recomb_map = NULL;
    double Ne = 0.25; /* default to 1/4 for coalescent time units. */

    if (self->tree_sequence != NULL) {
        PyErr_SetString(PyExc_ValueError, "tree_sequence already created");
        goto out;
    }
    if (!PyArg_ParseTuple(args, "O!O!|d",
                &SimulatorType, &sim,
                &RecombinationMapType, &recomb_map, &Ne)) {
        goto out;
    }
    if (Simulator_check_sim(sim) != 0) {
        goto out;
    }
    if (!msp_is_completed(sim->sim)) {
        PyErr_SetString(PyExc_ValueError, "Simulation not completed");
        goto out;
    }
    if (RecombinationMap_check_recomb_map(recomb_map) != 0) {
        goto out;
    }
    self->tree_sequence = PyMem_Malloc(sizeof(tree_sequence_t));
    if (self->tree_sequence == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->tree_sequence, 0, sizeof(tree_sequence_t));
    err = tree_sequence_create(self->tree_sequence, sim->sim,
            recomb_map->recomb_map, Ne);
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
        flags = MSP_ZLIB_COMPRESSION;
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
TreeSequence_load_records(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    int err;
    PyObject *ret = NULL;
    PyObject *py_records = NULL;
    PyObject *py_samples = NULL;
    PyObject *item;
    coalescence_record_t *records = NULL;
    sample_t *samples = NULL;
    Py_ssize_t sample_size;
    size_t num_records, j;
    static char *kwlist[] = {"records", "samples", NULL};

    if (self->tree_sequence != NULL) {
        PyErr_SetString(PyExc_ValueError, "TreeSequence already initialised");
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|O!", kwlist,
                &PyList_Type, &py_records, &PyList_Type, &py_samples)) {
        goto out;
    }
    self->tree_sequence = PyMem_Malloc(sizeof(tree_sequence_t));
    if (self->tree_sequence == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->tree_sequence, 0, sizeof(tree_sequence_t));
    num_records = PyList_Size(py_records);
    records = PyMem_Malloc(num_records * sizeof(coalescence_record_t));
    if (records == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(records, 0, num_records * sizeof(coalescence_record_t));
    for (j = 0; j < num_records; j++) {
        item = PyList_GetItem(py_records, j);
        if (parse_coalescence_record(item, &records[j]) != 0) {
            goto out;
        }
    }
    err = tree_sequence_load_records(
            self->tree_sequence, num_records, records);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    if (py_samples != NULL) {
        if (parse_samples(py_samples, &sample_size, &samples) != 0) {
            goto out;
        }
        err = tree_sequence_set_samples(
                self->tree_sequence, sample_size, samples);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    ret = Py_BuildValue("");
out:
    if (records != NULL) {
        for (j = 0; j < num_records; j++) {
            if (records[j].children != NULL) {
                PyMem_Free(records[j].children);
            }
        }
        PyMem_Free(records);
    }
    if (samples != NULL) {
        PyMem_Free(samples);
    }
    if (ret == NULL && self->tree_sequence != NULL) {
        /* Ensure that the state of the tree sequence is consistent */
        tree_sequence_free(self->tree_sequence);
        PyMem_Free(self->tree_sequence);
        self->tree_sequence = NULL;
    }
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

    if (self->tree_sequence != NULL) {
        PyErr_SetString(PyExc_ValueError, "TreeSequence already initialised");
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist,
                &path)) {
        goto out;
    }
    self->tree_sequence = PyMem_Malloc(sizeof(tree_sequence_t));
    if (self->tree_sequence == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->tree_sequence, 0, sizeof(tree_sequence_t));
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
    static char *kwlist[] = {"mutation_rate", "random_generator", NULL};
    mutgen_t *mutgen = NULL;
    double mutation_rate;
    RandomGenerator *random_generator = NULL;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dO!", kwlist,
            &mutation_rate, &RandomGeneratorType, &random_generator)) {
        goto out;
    }
    if (RandomGenerator_check_state(random_generator) != 0) {
        goto out;
    }
    mutgen = PyMem_Malloc(sizeof(mutgen_t));
    if (mutgen == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = mutgen_alloc(mutgen, self->tree_sequence, mutation_rate,
            random_generator->rng);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    err = mutgen_generate(mutgen);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    err = tree_sequence_set_mutations(self->tree_sequence, mutgen->num_mutations,
            mutgen->mutations);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    if (mutgen != NULL) {
        mutgen_free(mutgen);
        PyMem_Free(mutgen);
    }
    return ret;
}

static PyObject *
TreeSequence_set_mutations(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    int err;
    size_t j;
    PyObject *ret = NULL;
    PyObject *item, *node, *pos;
    PyObject *py_mutation_list = NULL;
    static char *kwlist[] = {"mutations", NULL};
    size_t num_mutations = 0;
    mutation_t *mutations = NULL;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
            &PyList_Type, &py_mutation_list)) {
        goto out;
    }
    num_mutations = PyList_Size(py_mutation_list);
    mutations = PyMem_Malloc(num_mutations * sizeof(mutation_t));
    for (j = 0; j < num_mutations; j++) {
        item = PyList_GetItem(py_mutation_list, j);
        if (!PyTuple_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "not a tuple");
            goto out;
        }
        if (PyTuple_Size(item) < 2) {
            PyErr_SetString(PyExc_ValueError,
                    "mutations must (node, pos, ...) tuples");
            goto out;
        }
        pos = PyTuple_GetItem(item, 0);
        node = PyTuple_GetItem(item, 1);
        if (!PyNumber_Check(pos)) {
            PyErr_SetString(PyExc_TypeError, "position must be a number");
            goto out;
        }
        if (!PyNumber_Check(node)) {
            PyErr_SetString(PyExc_TypeError, "node must be a number");
            goto out;
        }
        mutations[j].position = PyFloat_AsDouble(pos);
        mutations[j].node = (uint32_t) PyLong_AsLong(node);
    }
    err = tree_sequence_set_mutations(self->tree_sequence, num_mutations,
            mutations);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    if (mutations != NULL) {
        PyMem_Free(mutations);
    }
    return ret;
}

static PyObject *
TreeSequence_add_provenance_string(TreeSequence *self, PyObject *args)
{
    int err;
    char *s;
    PyObject *ret = NULL;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "s", &s)) {
        goto out;
    }
    if (strlen(s) == 0) {
        PyErr_SetString(PyExc_ValueError,
                "Empty string is not permitted for provenance.");
        goto out;
    }
    err = tree_sequence_add_provenance_string(self->tree_sequence, s);
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
    int order = MSP_ORDER_TIME;
    Py_ssize_t record_index, num_records;
    coalescence_record_t *cr;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n|i", &record_index, &order)) {
        goto out;
    }
    num_records = (Py_ssize_t) tree_sequence_get_num_coalescence_records(
        self->tree_sequence);
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
    ret = make_coalescence_record(cr);
out:
    return ret;
}

static PyObject *
TreeSequence_get_mutations(TreeSequence *self)
{
    PyObject *ret = NULL;
    mutation_t *mutations;
    size_t num_mutations;
    int err;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_mutations = tree_sequence_get_num_mutations(self->tree_sequence);
    err = tree_sequence_get_mutations(self->tree_sequence, &mutations);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = convert_mutations(mutations, num_mutations);
out:
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

static PyObject *
TreeSequence_get_sample(TreeSequence *self, PyObject *args)
{
    PyObject *ret = NULL;
    unsigned int node;
    sample_t sample;
    int population, err;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &node)) {
        goto out;
    }
    err = tree_sequence_get_sample(self->tree_sequence, node, &sample);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    population = sample.population_id;
    if (sample.population_id == MSP_NULL_POPULATION_ID) {
        population = -1;
    }
    ret = Py_BuildValue("id", population, sample.time);
out:
    return ret;
}

static PyObject *
TreeSequence_get_pairwise_diversity(TreeSequence *self, PyObject *args,
        PyObject *kwds)
{
    PyObject *ret = NULL;
    PyObject *py_samples = NULL;
    PyObject *item;
    static char *kwlist[] = {"samples", NULL};
    uint32_t *samples = NULL;
    size_t num_samples = 0;
    uint32_t j, n;
    double pi;
    int err;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
            &PyList_Type, &py_samples)) {
        goto out;
    }
    n = tree_sequence_get_sample_size(self->tree_sequence);
    num_samples = PyList_Size(py_samples);
    if (num_samples < 2) {
        PyErr_SetString(PyExc_ValueError, "Must provide at least 2 samples");
        goto out;
    }
    samples = PyMem_Malloc(num_samples * sizeof(uint32_t));
    if (samples == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    for (j = 0; j < num_samples; j++) {
        item = PyList_GetItem(py_samples, j);
        if (!PyNumber_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "sample id must be a number");
            goto out;
        }
        samples[j] = (uint32_t) PyLong_AsLong(item);
        if (samples[j] >= n) {
            PyErr_SetString(PyExc_ValueError,
                    "sample ids must be < sample_size");
            goto out;
        }
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
    {"load_records", (PyCFunction) TreeSequence_load_records,
        METH_VARARGS|METH_KEYWORDS,
        "Loads a tree sequence from the specified set of records"},
    {"generate_mutations", (PyCFunction) TreeSequence_generate_mutations,
        METH_VARARGS|METH_KEYWORDS,
        "Generates mutations under the infinite sites model"},
    {"set_mutations", (PyCFunction) TreeSequence_set_mutations,
        METH_VARARGS|METH_KEYWORDS,
        "Sets the mutations to the specified list of tuples."},
    {"add_provenance_string", (PyCFunction) TreeSequence_add_provenance_string,
        METH_VARARGS, "Appends a provenance string to the list."},
    {"get_provenance_strings", (PyCFunction) TreeSequence_get_provenance_strings,
        METH_NOARGS, "Returns the list of provenance strings."},
    {"get_mutations", (PyCFunction) TreeSequence_get_mutations,
        METH_NOARGS, "Returns the list of mutations"},
    {"get_record", (PyCFunction) TreeSequence_get_record, METH_VARARGS,
        "Returns the record at the specified index."},
    {"get_num_records", (PyCFunction) TreeSequence_get_num_records,
        METH_NOARGS, "Returns the number of coalescence records." },
    {"get_num_trees", (PyCFunction) TreeSequence_get_num_trees,
        METH_NOARGS, "Returns the number of trees in the tree sequence." },
    {"get_sequence_length", (PyCFunction) TreeSequence_get_sequence_length,
        METH_NOARGS, "Returns the sequence length in bases." },
    {"get_num_mutations", (PyCFunction) TreeSequence_get_num_mutations, METH_NOARGS,
        "Returns the number of loci" },
    {"get_num_nodes", (PyCFunction) TreeSequence_get_num_nodes, METH_NOARGS,
        "Returns the number of unique nodes in the tree sequence." },
    {"get_sample_size", (PyCFunction) TreeSequence_get_sample_size, METH_NOARGS,
        "Returns the sample size" },
    {"get_sample", (PyCFunction) TreeSequence_get_sample, METH_VARARGS,
        "Returns a dictionary describing the specified sample." },
    {"get_pairwise_diversity",
        (PyCFunction) TreeSequence_get_pairwise_diversity,
        METH_VARARGS|METH_KEYWORDS, "Returns the average pairwise diversity." },
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
SparseTree_check_bounds(SparseTree *self, unsigned int node)
{
    int ret = 0;
    if (node >= self->sparse_tree->num_nodes) {
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
    uint32_t *tracked_leaves = NULL;
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
    tracked_leaves = PyMem_Malloc(num_tracked_leaves * sizeof(uint32_t));
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
        tracked_leaves[j] = (uint32_t) PyLong_AsLong(item);
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

static PyObject *
SparseTree_get_parent(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    unsigned int node;
    uint32_t parent;

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
    ret = Py_BuildValue("i", (int) parent);
out:
    return ret;
}

static PyObject *
SparseTree_get_population(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    unsigned int node;
    int population;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &node)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, node)) {
        goto out;
    }
    population = self->sparse_tree->population[node];
    if (population == MSP_NULL_POPULATION_ID) {
        population = -1;
    }
    ret = Py_BuildValue("i", population);
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
    uint32_t *children, num_children;
    unsigned int node;
    int err;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &node)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, node)) {
        goto out;
    }
    err = sparse_tree_get_children(self->sparse_tree,
            (uint32_t) node, &num_children, &children);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    if (num_children == 0) {
        ret = Py_BuildValue("()");
    } else {
        ret = convert_children(children, num_children);
    }
out:
    return ret;
}

static PyObject *
SparseTree_get_mrca(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    int err;
    uint32_t mrca;
    unsigned int u, v;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "II", &u, &v)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, u)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, v)) {
        goto out;
    }
    err = sparse_tree_get_mrca(self->sparse_tree, (uint32_t) u,
            (uint32_t) v, &mrca);
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
    unsigned int node;
    uint32_t num_leaves;
    int err;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &node)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, node)) {
        goto out;
    }
    err = sparse_tree_get_num_leaves(self->sparse_tree, (uint32_t) node,
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
    unsigned int node;
    uint32_t num_tracked_leaves;
    int err;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &node)) {
        goto out;
    }
    if (SparseTree_check_bounds(self, node)) {
        goto out;
    }
    err = sparse_tree_get_num_tracked_leaves(self->sparse_tree, (uint32_t) node,
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
            children = convert_children(record->children, record->num_children);
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
            children = convert_children(record->children, record->num_children);
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
    static char *kwlist[] = {"tree_sequence", "genotypes_buffer", "as_char"};
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
        ret = convert_mutation(mutation);
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
    /* RandomGenerator type */
    RandomGeneratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&RandomGeneratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&RandomGeneratorType);
    PyModule_AddObject(module, "RandomGenerator", (PyObject *) &RandomGeneratorType);
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
    PyModule_AddObject(module, "TreeDiffIterator",
            (PyObject *) &TreeDiffIteratorType);
    /* LeafListIterator type */
    LeafListIteratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&LeafListIteratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&LeafListIteratorType);
    PyModule_AddObject(module, "LeafListIterator",
            (PyObject *) &LeafListIteratorType);
    /* NewickConverter type */
    NewickConverterType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&NewickConverterType) < 0) {
        INITERROR;
    }
    Py_INCREF(&NewickConverterType);
    PyModule_AddObject(module, "NewickConverter",
            (PyObject *) &NewickConverterType);
    /* VcfConverter type */
    VcfConverterType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&VcfConverterType) < 0) {
        INITERROR;
    }
    Py_INCREF(&VcfConverterType);
    PyModule_AddObject(module, "VcfConverter",
            (PyObject *) &VcfConverterType);
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
    PyModule_AddObject(module, "VariantGenerator",
            (PyObject *) &VariantGeneratorType);
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
    MsprimeLibraryError = PyErr_NewException("_msprime.LibraryError", NULL,
            NULL);
    Py_INCREF(MsprimeLibraryError);
    PyModule_AddObject(module, "LibraryError", MsprimeLibraryError);

    PyModule_AddIntConstant(module, "MSP_ORDER_TIME", MSP_ORDER_TIME);
    PyModule_AddIntConstant(module, "MSP_ORDER_LEFT", MSP_ORDER_LEFT);
    PyModule_AddIntConstant(module, "MSP_ORDER_RIGHT", MSP_ORDER_RIGHT);
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


