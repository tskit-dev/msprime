/*
** Copyright (C) 2014-2020 University of Oxford
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
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <float.h>

#include <gsl/gsl_version.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>

#include "msprime.h"
#include "likelihood.h"

/* We keep a reference to the gsl_error_handler so it can be restored if needed */
static gsl_error_handler_t *old_gsl_error_handler;

static PyObject *MsprimeInputError;
static PyObject *MsprimeLibraryError;

/* A lightweight wrapper for a table collection. This serves only as a wrapper
 * around a pointer and a way move to data in-and-out of the low level structures
 * via the canonical dictionary encoding.
 *
 * Copied from _tskitmodule.c 2020-08-07, 2df6c04
 */
typedef struct {
    PyObject_HEAD
    tsk_table_collection_t *tables;
} LightweightTableCollection;

typedef struct {
    PyObject_HEAD
    unsigned long seed;
    gsl_rng* rng;
} RandomGenerator;

/* TODO we should refactor some of the code for dealing with the
 * mutation_model in this base class (which currently does nothing).
 */
typedef struct {
    PyObject_HEAD
} BaseMutationModel;

typedef struct {
    PyObject_HEAD
    mutation_model_t *mutation_model;
} MatrixMutationModel;

typedef struct {
    PyObject_HEAD
    mutation_model_t *mutation_model;
} SLiMMutationModel;

typedef struct {
    PyObject_HEAD
    mutation_model_t *mutation_model;
} InfiniteAllelesMutationModel;

typedef struct {
    PyObject_HEAD
    msp_t *sim;
    RandomGenerator *random_generator;
    LightweightTableCollection *tables;
} Simulator;

static void
handle_library_error(int err)
{
    PyErr_SetString(MsprimeLibraryError, msp_strerror(err));
}

static void
handle_tskit_library_error(int err)
{
    PyErr_SetString(MsprimeLibraryError, tsk_strerror(err));
}

static void
handle_input_error(const char *section, int err)
{
    PyErr_Format(MsprimeInputError, "Input error in %s: %s", section, msp_strerror(err));
}

static int
double_PyArray_converter(PyObject *in, PyObject **converted)
{
    int ret = 0;
    PyObject *array = PyArray_FROMANY(in, NPY_FLOAT64, 1, 1, NPY_ARRAY_IN_ARRAY);

    if (array == NULL) {
        goto out;
    }
    *converted = array;
    ret = 1;
out:
    return ret;;
}

/*
 * Retrieves the PyObject* corresponding the specified key in the
 * specified dictionary.
 *
 * NB This returns a *borrowed reference*, so don't DECREF it!
 */
static PyObject *
get_dict_value(PyObject *dict, const char *key_str)
{
    PyObject *ret = NULL;

    ret = PyDict_GetItemString(dict, key_str);
    if (ret == NULL) {
        PyErr_Format(PyExc_ValueError, "'%s' not specified", key_str);
    }
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

static int
parse_rate_map(PyObject *py_rate_map, size_t *ret_size,
        PyArrayObject **ret_position, PyArrayObject **ret_rate)
{
    int ret = -1;
    PyObject *position = NULL;
    PyObject *rate = NULL;
    PyArrayObject *position_array = NULL;
    PyArrayObject *rate_array = NULL;
    npy_intp *dims, size;

    position = get_dict_value(py_rate_map, "position");
    if (position == NULL) {
        goto out;
    }
    rate = get_dict_value(py_rate_map, "rate");
    if (rate == NULL) {
        goto out;
    }
    position_array = (PyArrayObject *) PyArray_FROMANY(
            position, NPY_FLOAT64, 1, 1, NPY_ARRAY_IN_ARRAY);
    if (position_array == NULL) {
        goto out;
    }
    rate_array = (PyArrayObject *) PyArray_FROMANY(
            rate, NPY_FLOAT64, 1, 1, NPY_ARRAY_IN_ARRAY);
    if (rate_array == NULL) {
        goto out;
    }
    dims = PyArray_DIMS(rate_array);
    /* size is the number of intervals in the rate map, so the number of
     * positions must be 1+ this. */
    size = dims[0];
    dims = PyArray_DIMS(position_array);
    if (dims[0] != size + 1) {
        PyErr_SetString(PyExc_ValueError,
                "The position array must be one larger than rates");
        goto out;
    }
    *ret_size = size;
    *ret_position = position_array;
    *ret_rate = rate_array;
    position_array = NULL;
    rate_array = NULL;
    ret = 0;
out:
    Py_XDECREF(position_array);
    Py_XDECREF(rate_array);
    return ret;
}

/*===================================================================
 * General table code.
 *===================================================================
 */

/***********************************************
 * Start of code copied  _tskitmodule.
 ***********************************************/

/* NOTE: this code was copied from _tskitmodule as the efficient way to
 * import and export TableCollection data. It is unlikely to change
 * much over time, but if updates need to be made it would be better
 * to copy the code wholesale. The tests in ``test_dict_encoding.py``
 * are designed to test this code thoroughly, and also come from
 * tskit.
 *
 * Copied from _tskitmodule.c 2020-08-07, 2df6c04
 *
 * NOTE: we really should find a better mechanism for doing this:
 * perhaps a C header file that can be included by projects that
 * want to inferface between the C API and Python?
 */

static PyObject *
make_Py_Unicode_FromStringAndLength(const char *str, size_t length) {
    PyObject *ret = NULL;

    /* Py_BuildValue returns Py_None for zero length, we would rather
       return a zero-length string */
    if (length == 0) {
        ret = PyUnicode_FromString("");
    } else {
        ret = Py_BuildValue("s#", str, length);
    }
    return ret;
}

/*===================================================================
 * General table code.
 *===================================================================
 */

/*
 * Retrieves the PyObject* corresponding the specified key in the
 * specified dictionary. If required is true, raise a TypeError if the
 * value is None or absent.
 *
 * NB This returns a *borrowed reference*, so don't DECREF it!
 */
static PyObject *
get_table_dict_value(PyObject *dict, const char *key_str, bool required)
{
    PyObject *ret = NULL;

    ret = PyDict_GetItemString(dict, key_str);
    if (ret == NULL) {
        ret = Py_None;
    }
    if (required && ret == Py_None) {
        PyErr_Format(PyExc_TypeError, "'%s' is required", key_str);
        ret = NULL;
    }
    return ret;
}

static PyArrayObject *
table_read_column_array(PyObject *input, int npy_type, size_t *num_rows, bool check_num_rows)
{
    PyArrayObject *ret = NULL;
    PyArrayObject *array = NULL;
    npy_intp *shape;

    array = (PyArrayObject *) PyArray_FROMANY(input, npy_type, 1, 1, NPY_ARRAY_IN_ARRAY);
    if (array == NULL) {
        goto out;
    }
    shape = PyArray_DIMS(array);
    if (check_num_rows) {
        if (*num_rows != (size_t) shape[0]) {
            PyErr_SetString(PyExc_ValueError, "Input array dimensions must be equal.");
            goto out;
        }
    } else {
        *num_rows = (size_t) shape[0];
    }
    ret = array;
    array = NULL;
out:
    Py_XDECREF(array);
    return ret;
}

static PyArrayObject *
table_read_offset_array(PyObject *input, size_t *num_rows, size_t length, bool check_num_rows)
{
    PyArrayObject *ret = NULL;
    PyArrayObject *array = NULL;
    npy_intp *shape;
    uint32_t *data;

    array = (PyArrayObject *) PyArray_FROMANY(input, NPY_UINT32, 1, 1, NPY_ARRAY_IN_ARRAY);
    if (array == NULL) {
        goto out;
    }
    shape = PyArray_DIMS(array);
    if (! check_num_rows) {
        *num_rows = shape[0];
        if (*num_rows == 0) {
            PyErr_SetString(PyExc_ValueError, "Offset arrays must have at least one element");
            goto out;
        }
        *num_rows -= 1;
    }
    if (shape[0] != (npy_intp) (*num_rows + 1)) {
        PyErr_SetString(PyExc_ValueError, "offset columns must have n + 1 rows.");
        goto out;
    }
    data = PyArray_DATA(array);
    if (data[*num_rows] != (uint32_t) length) {
        PyErr_SetString(PyExc_ValueError, "Bad offset column encoding");
        goto out;
    }
    ret = array;
out:
    if (ret == NULL) {
        Py_XDECREF(array);
    }
    return ret;
}

static const char *
parse_metadata_schema_arg(PyObject *arg, Py_ssize_t* metadata_schema_length)
{
    const char *ret = NULL;
    if (arg == NULL) {
        PyErr_Format(
            PyExc_AttributeError,
            "Cannot del metadata_schema, set to empty string (\"\") to clear.");
        goto out;
    }
    ret = PyUnicode_AsUTF8AndSize(arg, metadata_schema_length);
    if (ret == NULL) {
        goto out;
    }
out:
    return ret;
}

static int
parse_individual_table_dict(tsk_individual_table_t *table, PyObject *dict, bool clear_table)
{
    int err;
    int ret = -1;
    size_t num_rows, metadata_length, location_length;
    char *metadata_data = NULL;
    double *location_data = NULL;
    uint32_t *metadata_offset_data = NULL;
    uint32_t *location_offset_data = NULL;
    PyObject *flags_input = NULL;
    PyArrayObject *flags_array = NULL;
    PyObject *location_input = NULL;
    PyArrayObject *location_array = NULL;
    PyObject *location_offset_input = NULL;
    PyArrayObject *location_offset_array = NULL;
    PyObject *metadata_input = NULL;
    PyArrayObject *metadata_array = NULL;
    PyObject *metadata_offset_input = NULL;
    PyArrayObject *metadata_offset_array = NULL;
    PyObject *metadata_schema_input = NULL;
    const char *metadata_schema = NULL;
    Py_ssize_t metadata_schema_length = 0;

    /* Get the input values */
    flags_input = get_table_dict_value(dict, "flags", true);
    if (flags_input == NULL) {
        goto out;
    }
    location_input = get_table_dict_value(dict, "location", false);
    if (location_input == NULL) {
        goto out;
    }
    location_offset_input = get_table_dict_value(dict, "location_offset", false);
    if (location_offset_input == NULL) {
        goto out;
    }
    metadata_input = get_table_dict_value(dict, "metadata", false);
    if (metadata_input == NULL) {
        goto out;
    }
    metadata_offset_input = get_table_dict_value(dict, "metadata_offset", false);
    if (metadata_offset_input == NULL) {
        goto out;
    }
    metadata_schema_input = get_table_dict_value(dict, "metadata_schema", false);
    if (metadata_schema_input == NULL) {
        goto out;
    }

    /* Pull out the arrays */
    flags_array = table_read_column_array(flags_input, NPY_UINT32, &num_rows, false);
    if (flags_array == NULL) {
        goto out;
    }
    if ((location_input == Py_None) != (location_offset_input == Py_None)) {
        PyErr_SetString(PyExc_TypeError,
                "location and location_offset must be specified together");
        goto out;
    }
    if (location_input != Py_None) {
        location_array = table_read_column_array(location_input, NPY_FLOAT64,
                &location_length, false);
        if (location_array == NULL) {
            goto out;
        }
        location_data = PyArray_DATA(location_array);
        location_offset_array = table_read_offset_array(location_offset_input, &num_rows,
                location_length, true);
        if (location_offset_array == NULL) {
            goto out;
        }
        location_offset_data = PyArray_DATA(location_offset_array);
    }
    if ((metadata_input == Py_None) != (metadata_offset_input == Py_None)) {
        PyErr_SetString(PyExc_TypeError,
                "metadata and metadata_offset must be specified together");
        goto out;
    }
    if (metadata_input != Py_None) {
        metadata_array = table_read_column_array(metadata_input, NPY_INT8,
                &metadata_length, false);
        if (metadata_array == NULL) {
            goto out;
        }
        metadata_data = PyArray_DATA(metadata_array);
        metadata_offset_array = table_read_offset_array(metadata_offset_input, &num_rows,
                metadata_length, true);
        if (metadata_offset_array == NULL) {
            goto out;
        }
        metadata_offset_data = PyArray_DATA(metadata_offset_array);
    }

    if (metadata_schema_input != Py_None) {
        metadata_schema = parse_metadata_schema_arg(
            metadata_schema_input, &metadata_schema_length);
        if (metadata_schema == NULL) {
            goto out;
        }
        err = tsk_individual_table_set_metadata_schema(
            table, metadata_schema, metadata_schema_length);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }

    if (clear_table) {
        err = tsk_individual_table_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_individual_table_append_columns(table, num_rows,
            PyArray_DATA(flags_array),
            location_data, location_offset_data,
            metadata_data, metadata_offset_data);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(flags_array);
    Py_XDECREF(location_array);
    Py_XDECREF(location_offset_array);
    Py_XDECREF(metadata_array);
    Py_XDECREF(metadata_offset_array);
    return ret;
}

static int
parse_node_table_dict(tsk_node_table_t *table, PyObject *dict, bool clear_table)
{
    int err;
    int ret = -1;
    size_t num_rows, metadata_length;
    char *metadata_data = NULL;
    uint32_t *metadata_offset_data = NULL;
    void *population_data = NULL;
    void *individual_data = NULL;
    PyObject *time_input = NULL;
    PyArrayObject *time_array = NULL;
    PyObject *flags_input = NULL;
    PyArrayObject *flags_array = NULL;
    PyObject *population_input = NULL;
    PyArrayObject *population_array = NULL;
    PyObject *individual_input = NULL;
    PyArrayObject *individual_array = NULL;
    PyObject *metadata_input = NULL;
    PyArrayObject *metadata_array = NULL;
    PyObject *metadata_offset_input = NULL;
    PyArrayObject *metadata_offset_array = NULL;
    PyObject *metadata_schema_input = NULL;
    const char *metadata_schema = NULL;
    Py_ssize_t metadata_schema_length = 0;

    /* Get the input values */
    flags_input = get_table_dict_value(dict, "flags", true);
    if (flags_input == NULL) {
        goto out;
    }
    time_input = get_table_dict_value(dict, "time", true);
    if (time_input == NULL) {
        goto out;
    }
    population_input = get_table_dict_value(dict, "population", false);
    if (population_input == NULL) {
        goto out;
    }
    individual_input = get_table_dict_value(dict, "individual", false);
    if (individual_input == NULL) {
        goto out;
    }
    metadata_input = get_table_dict_value(dict, "metadata", false);
    if (metadata_input == NULL) {
        goto out;
    }
    metadata_offset_input = get_table_dict_value(dict, "metadata_offset", false);
    if (metadata_offset_input == NULL) {
        goto out;
    }
    metadata_schema_input = get_table_dict_value(dict, "metadata_schema", false);
    if (metadata_schema_input == NULL) {
        goto out;
    }

    /* Create the arrays */
    flags_array = table_read_column_array(flags_input, NPY_UINT32, &num_rows, false);
    if (flags_array == NULL) {
        goto out;
    }
    time_array = table_read_column_array(time_input, NPY_FLOAT64, &num_rows, true);
    if (time_array == NULL) {
        goto out;
    }
    if (population_input != Py_None) {
        population_array = table_read_column_array(population_input, NPY_INT32,
                &num_rows, true);
        if (population_array == NULL) {
            goto out;
        }
        population_data = PyArray_DATA(population_array);
    }
    if (individual_input != Py_None) {
        individual_array = table_read_column_array(individual_input, NPY_INT32,
                &num_rows, true);
        if (individual_array == NULL) {
            goto out;
        }
        individual_data = PyArray_DATA(individual_array);
    }
    if ((metadata_input == Py_None) != (metadata_offset_input == Py_None)) {
        PyErr_SetString(PyExc_TypeError,
                "metadata and metadata_offset must be specified together");
        goto out;
    }
    if (metadata_input != Py_None) {
        metadata_array = table_read_column_array(metadata_input, NPY_INT8,
                &metadata_length, false);
        if (metadata_array == NULL) {
            goto out;
        }
        metadata_data = PyArray_DATA(metadata_array);
        metadata_offset_array = table_read_offset_array(metadata_offset_input, &num_rows,
                metadata_length, true);
        if (metadata_offset_array == NULL) {
            goto out;
        }
        metadata_offset_data = PyArray_DATA(metadata_offset_array);
    }
    if (metadata_schema_input != Py_None) {
        metadata_schema = parse_metadata_schema_arg(
            metadata_schema_input, &metadata_schema_length);
        if (metadata_schema == NULL) {
            goto out;
        }
        err = tsk_node_table_set_metadata_schema(
            table, metadata_schema, metadata_schema_length);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }

    if (clear_table) {
        err = tsk_node_table_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_node_table_append_columns(table, num_rows,
            PyArray_DATA(flags_array), PyArray_DATA(time_array), population_data,
            individual_data, metadata_data, metadata_offset_data);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(flags_array);
    Py_XDECREF(time_array);
    Py_XDECREF(population_array);
    Py_XDECREF(individual_array);
    Py_XDECREF(metadata_array);
    Py_XDECREF(metadata_offset_array);
    return ret;
}

static int
parse_edge_table_dict(tsk_edge_table_t *table, PyObject *dict, bool clear_table)
{
    int ret = -1;
    int err;
    size_t num_rows = 0;
    size_t metadata_length;
    char *metadata_data = NULL;
    uint32_t *metadata_offset_data = NULL;
    PyObject *left_input = NULL;
    PyArrayObject *left_array = NULL;
    PyObject *right_input = NULL;
    PyArrayObject *right_array = NULL;
    PyObject *parent_input = NULL;
    PyArrayObject *parent_array = NULL;
    PyObject *child_input = NULL;
    PyArrayObject *child_array = NULL;
    PyObject *metadata_input = NULL;
    PyArrayObject *metadata_array = NULL;
    PyObject *metadata_offset_input = NULL;
    PyArrayObject *metadata_offset_array = NULL;
    PyObject *metadata_schema_input = NULL;
    const char *metadata_schema = NULL;
    Py_ssize_t metadata_schema_length = 0;

    /* Get the input values */
    left_input = get_table_dict_value(dict, "left", true);
    if (left_input == NULL) {
        goto out;
    }
    right_input = get_table_dict_value(dict, "right", true);
    if (right_input == NULL) {
        goto out;
    }
    parent_input = get_table_dict_value(dict, "parent", true);
    if (parent_input == NULL) {
        goto out;
    }
    child_input = get_table_dict_value(dict, "child", true);
    if (child_input == NULL) {
        goto out;
    }
    metadata_input = get_table_dict_value(dict, "metadata", false);
    if (metadata_input == NULL) {
        goto out;
    }
    metadata_offset_input = get_table_dict_value(dict, "metadata_offset", false);
    if (metadata_offset_input == NULL) {
        goto out;
    }
    metadata_schema_input = get_table_dict_value(dict, "metadata_schema", false);
    if (metadata_schema_input == NULL) {
        goto out;
    }


    /* Create the arrays */
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
    child_array = table_read_column_array(child_input, NPY_INT32, &num_rows, true);
    if (child_array == NULL) {
        goto out;
    }
    if ((metadata_input == Py_None) != (metadata_offset_input == Py_None)) {
        PyErr_SetString(PyExc_TypeError,
                "metadata and metadata_offset must be specified together");
        goto out;
    }
    if (metadata_input != Py_None) {
        metadata_array = table_read_column_array(metadata_input, NPY_INT8,
                &metadata_length, false);
        if (metadata_array == NULL) {
            goto out;
        }
        metadata_data = PyArray_DATA(metadata_array);
        metadata_offset_array = table_read_offset_array(metadata_offset_input, &num_rows,
                metadata_length, true);
        if (metadata_offset_array == NULL) {
            goto out;
        }
        metadata_offset_data = PyArray_DATA(metadata_offset_array);
    }
    if (metadata_schema_input != Py_None) {
        metadata_schema = parse_metadata_schema_arg(
            metadata_schema_input, &metadata_schema_length);
        if (metadata_schema == NULL) {
            goto out;
        }
        err = tsk_edge_table_set_metadata_schema(
            table, metadata_schema, metadata_schema_length);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }

    if (clear_table) {
        err = tsk_edge_table_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_edge_table_append_columns(table, num_rows,
            PyArray_DATA(left_array), PyArray_DATA(right_array),
            PyArray_DATA(parent_array), PyArray_DATA(child_array),
            metadata_data, metadata_offset_data);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(left_array);
    Py_XDECREF(right_array);
    Py_XDECREF(parent_array);
    Py_XDECREF(child_array);
    Py_XDECREF(metadata_array);
    Py_XDECREF(metadata_offset_array);
    return ret;
}

static int
parse_migration_table_dict(tsk_migration_table_t *table, PyObject *dict, bool clear_table)
{
    int err;
    int ret = -1;
    size_t num_rows;
    size_t metadata_length;
    char *metadata_data = NULL;
    uint32_t *metadata_offset_data = NULL;
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
    PyObject *metadata_input = NULL;
    PyArrayObject *metadata_array = NULL;
    PyObject *metadata_offset_input = NULL;
    PyArrayObject *metadata_offset_array = NULL;
    PyObject *metadata_schema_input = NULL;
    const char *metadata_schema = NULL;
    Py_ssize_t metadata_schema_length = 0;

    /* Get the input values */
    left_input = get_table_dict_value(dict, "left", true);
    if (left_input == NULL) {
        goto out;
    }
    right_input = get_table_dict_value(dict, "right", true);
    if (right_input == NULL) {
        goto out;
    }
    node_input = get_table_dict_value(dict, "node", true);
    if (node_input == NULL) {
        goto out;
    }
    source_input = get_table_dict_value(dict, "source", true);
    if (source_input == NULL) {
        goto out;
    }
    dest_input = get_table_dict_value(dict, "dest", true);
    if (dest_input == NULL) {
        goto out;
    }
    time_input = get_table_dict_value(dict, "time", true);
    if (time_input == NULL) {
        goto out;
    }
    metadata_input = get_table_dict_value(dict, "metadata", false);
    if (metadata_input == NULL) {
        goto out;
    }
    metadata_offset_input = get_table_dict_value(dict, "metadata_offset", false);
    if (metadata_offset_input == NULL) {
        goto out;
    }
    metadata_schema_input = get_table_dict_value(dict, "metadata_schema", false);
    if (metadata_schema_input == NULL) {
        goto out;
    }

    /* Build the arrays */
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
    if ((metadata_input == Py_None) != (metadata_offset_input == Py_None)) {
        PyErr_SetString(PyExc_TypeError,
                "metadata and metadata_offset must be specified together");
        goto out;
    }
    if (metadata_input != Py_None) {
        metadata_array = table_read_column_array(metadata_input, NPY_INT8,
                &metadata_length, false);
        if (metadata_array == NULL) {
            goto out;
        }
        metadata_data = PyArray_DATA(metadata_array);
        metadata_offset_array = table_read_offset_array(metadata_offset_input, &num_rows,
                metadata_length, true);
        if (metadata_offset_array == NULL) {
            goto out;
        }
        metadata_offset_data = PyArray_DATA(metadata_offset_array);
    }
    if (metadata_schema_input != Py_None) {
        metadata_schema = parse_metadata_schema_arg(
            metadata_schema_input, &metadata_schema_length);
        if (metadata_schema == NULL) {
            goto out;
        }
        err = tsk_migration_table_set_metadata_schema(
            table, metadata_schema, metadata_schema_length);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }

    if (clear_table) {
        err = tsk_migration_table_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_migration_table_append_columns(table, num_rows,
        PyArray_DATA(left_array), PyArray_DATA(right_array), PyArray_DATA(node_array),
        PyArray_DATA(source_array), PyArray_DATA(dest_array), PyArray_DATA(time_array),
        metadata_data, metadata_offset_data);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(left_array);
    Py_XDECREF(right_array);
    Py_XDECREF(node_array);
    Py_XDECREF(source_array);
    Py_XDECREF(dest_array);
    Py_XDECREF(time_array);
    Py_XDECREF(metadata_array);
    Py_XDECREF(metadata_offset_array);
    return ret;
}

static int
parse_site_table_dict(tsk_site_table_t *table, PyObject *dict, bool clear_table)
{
    int err;
    int ret = -1;
    size_t num_rows = 0;
    size_t ancestral_state_length, metadata_length;
    PyObject *position_input = NULL;
    PyArrayObject *position_array = NULL;
    PyObject *ancestral_state_input = NULL;
    PyArrayObject *ancestral_state_array = NULL;
    PyObject *ancestral_state_offset_input = NULL;
    PyArrayObject *ancestral_state_offset_array = NULL;
    PyObject *metadata_input = NULL;
    PyArrayObject *metadata_array = NULL;
    PyObject *metadata_offset_input = NULL;
    PyArrayObject *metadata_offset_array = NULL;
    char *metadata_data;
    uint32_t *metadata_offset_data;
    PyObject *metadata_schema_input = NULL;
    const char *metadata_schema = NULL;
    Py_ssize_t metadata_schema_length = 0;


    /* Get the input values */
    position_input = get_table_dict_value(dict, "position", true);
    if (position_input == NULL) {
        goto out;
    }
    ancestral_state_input = get_table_dict_value(dict, "ancestral_state", true);
    if (ancestral_state_input == NULL) {
        goto out;
    }
    ancestral_state_offset_input = get_table_dict_value(dict, "ancestral_state_offset", true);
    if (ancestral_state_offset_input == NULL) {
        goto out;
    }
    metadata_input = get_table_dict_value(dict, "metadata", false);
    if (metadata_input == NULL) {
        goto out;
    }
    metadata_offset_input = get_table_dict_value(dict, "metadata_offset", false);
    if (metadata_offset_input == NULL) {
        goto out;
    }
    metadata_schema_input = get_table_dict_value(dict, "metadata_schema", false);
    if (metadata_schema_input == NULL) {
        goto out;
    }


    /* Get the arrays */
    position_array = table_read_column_array(position_input, NPY_FLOAT64, &num_rows, false);
    if (position_array == NULL) {
        goto out;
    }
    ancestral_state_array = table_read_column_array(ancestral_state_input, NPY_INT8,
            &ancestral_state_length, false);
    if (ancestral_state_array == NULL) {
        goto out;
    }
    ancestral_state_offset_array = table_read_offset_array(ancestral_state_offset_input,
            &num_rows, ancestral_state_length, true);
    if (ancestral_state_offset_array == NULL) {
        goto out;
    }

    metadata_data = NULL;
    metadata_offset_data = NULL;
    if ((metadata_input == Py_None) != (metadata_offset_input == Py_None)) {
        PyErr_SetString(PyExc_TypeError,
                "metadata and metadata_offset must be specified together");
        goto out;
    }
    if (metadata_input != Py_None) {
        metadata_array = table_read_column_array(metadata_input, NPY_INT8,
                &metadata_length, false);
        if (metadata_array == NULL) {
            goto out;
        }
        metadata_data = PyArray_DATA(metadata_array);
        metadata_offset_array = table_read_offset_array(metadata_offset_input, &num_rows,
                metadata_length, false);
        if (metadata_offset_array == NULL) {
            goto out;
        }
        metadata_offset_data = PyArray_DATA(metadata_offset_array);
    }
    if (metadata_schema_input != Py_None) {
        metadata_schema = parse_metadata_schema_arg(
            metadata_schema_input, &metadata_schema_length);
        if (metadata_schema == NULL) {
            goto out;
        }
        err = tsk_site_table_set_metadata_schema(
            table, metadata_schema, metadata_schema_length);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }

    if (clear_table) {
        err = tsk_site_table_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_site_table_append_columns(table, num_rows,
        PyArray_DATA(position_array), PyArray_DATA(ancestral_state_array),
        PyArray_DATA(ancestral_state_offset_array), metadata_data, metadata_offset_data);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(position_array);
    Py_XDECREF(ancestral_state_array);
    Py_XDECREF(ancestral_state_offset_array);
    Py_XDECREF(metadata_array);
    Py_XDECREF(metadata_offset_array);
    return ret;
}

static int
parse_mutation_table_dict(tsk_mutation_table_t *table, PyObject *dict, bool clear_table)
{
    int err;
    int ret = -1;
    size_t num_rows = 0;
    size_t derived_state_length = 0;
    size_t metadata_length = 0;
    PyObject *site_input = NULL;
    PyArrayObject *site_array = NULL;
    PyObject *derived_state_input = NULL;
    PyArrayObject *derived_state_array = NULL;
    PyObject *derived_state_offset_input = NULL;
    PyArrayObject *derived_state_offset_array = NULL;
    PyObject *node_input = NULL;
    PyArrayObject *node_array = NULL;
    PyObject *time_input = NULL;
    PyArrayObject *time_array = NULL;
    double *time_data;
    PyObject *parent_input = NULL;
    PyArrayObject *parent_array = NULL;
    tsk_id_t *parent_data;
    PyObject *metadata_input = NULL;
    PyArrayObject *metadata_array = NULL;
    PyObject *metadata_offset_input = NULL;
    PyArrayObject *metadata_offset_array = NULL;
    char *metadata_data;
    uint32_t *metadata_offset_data;
    PyObject *metadata_schema_input = NULL;
    const char *metadata_schema = NULL;
    Py_ssize_t metadata_schema_length = 0;

    /* Get the input values */
    site_input = get_table_dict_value(dict, "site", true);
    if (site_input == NULL) {
        goto out;
    }
    node_input = get_table_dict_value(dict, "node", true);
    if (node_input == NULL) {
        goto out;
    }
    parent_input = get_table_dict_value(dict, "parent", false);
    if (parent_input == NULL) {
        goto out;
    }
    time_input = get_table_dict_value(dict, "time", false);
    if (time_input == NULL) {
        goto out;
    }
    derived_state_input = get_table_dict_value(dict, "derived_state", true);
    if (derived_state_input == NULL) {
        goto out;
    }
    derived_state_offset_input = get_table_dict_value(dict, "derived_state_offset", true);
    if (derived_state_offset_input == NULL) {
        goto out;
    }
    metadata_input = get_table_dict_value(dict, "metadata", false);
    if (metadata_input == NULL) {
        goto out;
    }
    metadata_offset_input = get_table_dict_value(dict, "metadata_offset", false);
    if (metadata_offset_input == NULL) {
        goto out;
    }
    metadata_schema_input = get_table_dict_value(dict, "metadata_schema", false);
    if (metadata_schema_input == NULL) {
        goto out;
    }

    /* Get the arrays */
    site_array = table_read_column_array(site_input, NPY_INT32, &num_rows, false);
    if (site_array == NULL) {
        goto out;
    }
    derived_state_array = table_read_column_array(derived_state_input, NPY_INT8,
            &derived_state_length, false);
    if (derived_state_array == NULL) {
        goto out;
    }
    derived_state_offset_array = table_read_offset_array(derived_state_offset_input,
            &num_rows, derived_state_length, true);
    if (derived_state_offset_array == NULL) {
        goto out;
    }
    node_array = table_read_column_array(node_input, NPY_INT32, &num_rows, true);
    if (node_array == NULL) {
        goto out;
    }

    time_data = NULL;
    if (time_input != Py_None) {
        time_array = table_read_column_array(time_input, NPY_FLOAT64, &num_rows, true);
        if (time_array == NULL) {
            goto out;
        }
        time_data = PyArray_DATA(time_array);
    }

    parent_data = NULL;
    if (parent_input != Py_None) {
        parent_array = table_read_column_array(parent_input, NPY_INT32, &num_rows, true);
        if (parent_array == NULL) {
            goto out;
        }
        parent_data = PyArray_DATA(parent_array);
    }

    metadata_data = NULL;
    metadata_offset_data = NULL;
    if ((metadata_input == Py_None) != (metadata_offset_input == Py_None)) {
        PyErr_SetString(PyExc_TypeError,
                "metadata and metadata_offset must be specified together");
        goto out;
    }
    if (metadata_input != Py_None) {
        metadata_array = table_read_column_array(metadata_input, NPY_INT8,
                &metadata_length, false);
        if (metadata_array == NULL) {
            goto out;
        }
        metadata_data = PyArray_DATA(metadata_array);
        metadata_offset_array = table_read_offset_array(metadata_offset_input, &num_rows,
                metadata_length, false);
        if (metadata_offset_array == NULL) {
            goto out;
        }
        metadata_offset_data = PyArray_DATA(metadata_offset_array);
    }
    if (metadata_schema_input != Py_None) {
        metadata_schema = parse_metadata_schema_arg(
            metadata_schema_input, &metadata_schema_length);
        if (metadata_schema == NULL) {
            goto out;
        }
        err = tsk_mutation_table_set_metadata_schema(
            table, metadata_schema, metadata_schema_length);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }

    if (clear_table) {
        err = tsk_mutation_table_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_mutation_table_append_columns(table, num_rows,
            PyArray_DATA(site_array), PyArray_DATA(node_array),
            parent_data, time_data, PyArray_DATA(derived_state_array),
            PyArray_DATA(derived_state_offset_array),
            metadata_data, metadata_offset_data);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(site_array);
    Py_XDECREF(derived_state_array);
    Py_XDECREF(derived_state_offset_array);
    Py_XDECREF(metadata_array);
    Py_XDECREF(metadata_offset_array);
    Py_XDECREF(node_array);
    Py_XDECREF(parent_array);
    Py_XDECREF(time_array);
    return ret;
}

static int
parse_population_table_dict(tsk_population_table_t *table, PyObject *dict, bool clear_table)
{
    int err;
    int ret = -1;
    size_t num_rows, metadata_length;
    PyObject *metadata_input = NULL;
    PyArrayObject *metadata_array = NULL;
    PyObject *metadata_offset_input = NULL;
    PyArrayObject *metadata_offset_array = NULL;
    PyObject *metadata_schema_input = NULL;
    const char *metadata_schema = NULL;
    Py_ssize_t metadata_schema_length = 0;

    /* Get the inputs */
    metadata_input = get_table_dict_value(dict, "metadata", true);
    if (metadata_input == NULL) {
        goto out;
    }
    metadata_offset_input = get_table_dict_value(dict, "metadata_offset", true);
    if (metadata_offset_input == NULL) {
        goto out;
    }
    metadata_schema_input = get_table_dict_value(dict, "metadata_schema", false);
    if (metadata_schema_input == NULL) {
        goto out;
    }

    /* Get the arrays */
    metadata_array = table_read_column_array(metadata_input, NPY_INT8,
            &metadata_length, false);
    if (metadata_array == NULL) {
        goto out;
    }
    metadata_offset_array = table_read_offset_array(metadata_offset_input, &num_rows,
            metadata_length, false);
    if (metadata_offset_array == NULL) {
        goto out;
    }
    if (metadata_schema_input != Py_None) {
        metadata_schema = parse_metadata_schema_arg(
            metadata_schema_input, &metadata_schema_length);
        if (metadata_schema == NULL) {
            goto out;
        }
        err = tsk_population_table_set_metadata_schema(
            table, metadata_schema, metadata_schema_length);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }

    if (clear_table) {
        err = tsk_population_table_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_population_table_append_columns(table, num_rows,
            PyArray_DATA(metadata_array), PyArray_DATA(metadata_offset_array));
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(metadata_array);
    Py_XDECREF(metadata_offset_array);
    return ret;
}

static int
parse_provenance_table_dict(tsk_provenance_table_t *table, PyObject *dict, bool clear_table)
{
    int err;
    int ret = -1;
    size_t num_rows, timestamp_length, record_length;
    PyObject *timestamp_input = NULL;
    PyArrayObject *timestamp_array = NULL;
    PyObject *timestamp_offset_input = NULL;
    PyArrayObject *timestamp_offset_array = NULL;
    PyObject *record_input = NULL;
    PyArrayObject *record_array = NULL;
    PyObject *record_offset_input = NULL;
    PyArrayObject *record_offset_array = NULL;

    /* Get the inputs */
    timestamp_input = get_table_dict_value(dict, "timestamp", true);
    if (timestamp_input == NULL) {
        goto out;
    }
    timestamp_offset_input = get_table_dict_value(dict, "timestamp_offset", true);
    if (timestamp_offset_input == NULL) {
        goto out;
    }
    record_input = get_table_dict_value(dict, "record", true);
    if (record_input == NULL) {
        goto out;
    }
    record_offset_input = get_table_dict_value(dict, "record_offset", true);
    if (record_offset_input == NULL) {
        goto out;
    }

    timestamp_array = table_read_column_array(timestamp_input, NPY_INT8,
            &timestamp_length, false);
    if (timestamp_array == NULL) {
        goto out;
    }
    timestamp_offset_array = table_read_offset_array(timestamp_offset_input, &num_rows,
            timestamp_length, false);
    if (timestamp_offset_array == NULL) {
        goto out;
    }
    record_array = table_read_column_array(record_input, NPY_INT8,
            &record_length, false);
    if (record_array == NULL) {
        goto out;
    }
    record_offset_array = table_read_offset_array(record_offset_input, &num_rows,
            record_length, true);
    if (record_offset_array == NULL) {
        goto out;
    }

    if (clear_table) {
        err = tsk_provenance_table_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_provenance_table_append_columns(table, num_rows,
            PyArray_DATA(timestamp_array), PyArray_DATA(timestamp_offset_array),
            PyArray_DATA(record_array), PyArray_DATA(record_offset_array));
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(timestamp_array);
    Py_XDECREF(timestamp_offset_array);
    Py_XDECREF(record_array);
    Py_XDECREF(record_offset_array);
    return ret;
}

static int
parse_table_collection_dict(tsk_table_collection_t *tables, PyObject *tables_dict)
{
    int ret = -1;
    PyObject *value = NULL;
    int err;
    char *metadata = NULL;
    const char *metadata_schema = NULL;
    Py_ssize_t metadata_length, metadata_schema_length;

    value = get_table_dict_value(tables_dict, "sequence_length", true);
    if (value == NULL) {
        goto out;
    }
    if (!PyNumber_Check(value)) {
        PyErr_Format(PyExc_TypeError, "'sequence_length' is not number");
        goto out;
    }
    tables->sequence_length = PyFloat_AsDouble(value);

    /* metadata_schema */
    value = get_table_dict_value(tables_dict, "metadata_schema", false);
    if (value == NULL) {
        goto out;
    }
    if (value != Py_None) {
        if (!PyUnicode_Check(value)) {
            PyErr_Format(PyExc_TypeError, "'metadata_schema' is not a string");
            goto out;
        }
        metadata_schema = parse_metadata_schema_arg(value, &metadata_schema_length);
        if (metadata_schema == NULL) {
            goto out;
        }
        err = tsk_table_collection_set_metadata_schema(
            tables, metadata_schema, metadata_schema_length);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }

    /* metadata */
    value = get_table_dict_value(tables_dict, "metadata", false);
    if (value == NULL) {
        goto out;
    }
    if (value != Py_None) {
        if (!PyBytes_Check(value)) {
            PyErr_Format(PyExc_TypeError, "'metadata' is not bytes");
            goto out;
        }
        err = PyBytes_AsStringAndSize(value, &metadata, &metadata_length);
        if (err != 0) {
            goto out;
        }
        err = tsk_table_collection_set_metadata(
            tables, metadata, metadata_length);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }

    /* individuals */
    value = get_table_dict_value(tables_dict, "individuals", true);
    if (value == NULL) {
        goto out;
    }
    if (!PyDict_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "not a dictionary");
        goto out;
    }
    if (parse_individual_table_dict(&tables->individuals, value, true) != 0) {
        goto out;
    }

    /* nodes */
    value = get_table_dict_value(tables_dict, "nodes", true);
    if (value == NULL) {
        goto out;
    }
    if (!PyDict_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "not a dictionary");
        goto out;
    }
    if (parse_node_table_dict(&tables->nodes, value, true) != 0) {
        goto out;
    }

    /* edges */
    value = get_table_dict_value(tables_dict, "edges", true);
    if (value == NULL) {
        goto out;
    }
    if (!PyDict_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "not a dictionary");
        goto out;
    }
    if (parse_edge_table_dict(&tables->edges, value, true) != 0) {
        goto out;
    }

    /* migrations */
    value = get_table_dict_value(tables_dict, "migrations", true);
    if (value == NULL) {
        goto out;
    }
    if (!PyDict_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "not a dictionary");
        goto out;
    }
    if (parse_migration_table_dict(&tables->migrations, value, true) != 0) {
        goto out;
    }

    /* sites */
    value = get_table_dict_value(tables_dict, "sites", true);
    if (value == NULL) {
        goto out;
    }
    if (!PyDict_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "not a dictionary");
        goto out;
    }
    if (parse_site_table_dict(&tables->sites, value, true) != 0) {
        goto out;
    }

    /* mutations */
    value = get_table_dict_value(tables_dict, "mutations", true);
    if (value == NULL) {
        goto out;
    }
    if (!PyDict_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "not a dictionary");
        goto out;
    }
    if (parse_mutation_table_dict(&tables->mutations, value, true) != 0) {
        goto out;
    }

    /* populations */
    value = get_table_dict_value(tables_dict, "populations", true);
    if (value == NULL) {
        goto out;
    }
    if (!PyDict_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "not a dictionary");
        goto out;
    }
    if (parse_population_table_dict(&tables->populations, value, true) != 0) {
        goto out;
    }

    /* provenances */
    value = get_table_dict_value(tables_dict, "provenances", true);
    if (value == NULL) {
        goto out;
    }
    if (!PyDict_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "not a dictionary");
        goto out;
    }
    if (parse_provenance_table_dict(&tables->provenances, value, true) != 0) {
        goto out;
    }

    ret = 0;
out:
    return ret;
}

static int
write_table_arrays(tsk_table_collection_t *tables, PyObject *dict)
{
    struct table_col {
        const char *name;
        void *data;
        npy_intp num_rows;
        int type;
    };
    struct table_desc {
        const char *name;
        struct table_col *cols;
        char *metadata_schema;
        tsk_size_t metadata_schema_length;
    };
    int ret = -1;
    PyObject *array = NULL;
    PyObject *table_dict = NULL;
    size_t j;
    struct table_col *col;

    struct table_col individual_cols[] = {
        {"flags",
            (void *) tables->individuals.flags, tables->individuals.num_rows, NPY_UINT32},
        {"location",
            (void *) tables->individuals.location, tables->individuals.location_length,
            NPY_FLOAT64},
        {"location_offset",
            (void *) tables->individuals.location_offset, tables->individuals.num_rows + 1,
            NPY_UINT32},
        {"metadata",
            (void *) tables->individuals.metadata, tables->individuals.metadata_length,
            NPY_INT8},
        {"metadata_offset",
            (void *) tables->individuals.metadata_offset, tables->individuals.num_rows + 1,
            NPY_UINT32},
        {NULL},
    };

    struct table_col node_cols[] = {
        {"time",
            (void *) tables->nodes.time, tables->nodes.num_rows, NPY_FLOAT64},
        {"flags",
            (void *) tables->nodes.flags, tables->nodes.num_rows, NPY_UINT32},
        {"population",
            (void *) tables->nodes.population, tables->nodes.num_rows, NPY_INT32},
        {"individual",
            (void *) tables->nodes.individual, tables->nodes.num_rows, NPY_INT32},
        {"metadata",
            (void *) tables->nodes.metadata, tables->nodes.metadata_length, NPY_INT8},
        {"metadata_offset",
            (void *) tables->nodes.metadata_offset, tables->nodes.num_rows + 1, NPY_UINT32},
        {NULL},
    };

    struct table_col edge_cols[] = {
        {"left", (void *) tables->edges.left, tables->edges.num_rows, NPY_FLOAT64},
        {"right", (void *) tables->edges.right, tables->edges.num_rows, NPY_FLOAT64},
        {"parent", (void *) tables->edges.parent, tables->edges.num_rows, NPY_INT32},
        {"child", (void *) tables->edges.child, tables->edges.num_rows, NPY_INT32},
        {"metadata",
            (void *) tables->edges.metadata, tables->edges.metadata_length, NPY_INT8},
        {"metadata_offset",
            (void *) tables->edges.metadata_offset, tables->edges.num_rows + 1, NPY_UINT32},
        {NULL},
    };

    struct table_col migration_cols[] = {
        {"left",
            (void *) tables->migrations.left, tables->migrations.num_rows,  NPY_FLOAT64},
        {"right",
            (void *) tables->migrations.right, tables->migrations.num_rows,  NPY_FLOAT64},
        {"node",
            (void *) tables->migrations.node, tables->migrations.num_rows,  NPY_INT32},
        {"source",
            (void *) tables->migrations.source, tables->migrations.num_rows,  NPY_INT32},
        {"dest",
            (void *) tables->migrations.dest, tables->migrations.num_rows,  NPY_INT32},
        {"time",
            (void *) tables->migrations.time, tables->migrations.num_rows,  NPY_FLOAT64},
        {"metadata",
            (void *) tables->migrations.metadata, tables->migrations.metadata_length, NPY_INT8},
        {"metadata_offset",
            (void *) tables->migrations.metadata_offset, tables->migrations.num_rows + 1, NPY_UINT32},

        {NULL},
    };

    struct table_col site_cols[] = {
        {"position",
            (void *) tables->sites.position, tables->sites.num_rows, NPY_FLOAT64},
        {"ancestral_state",
            (void *) tables->sites.ancestral_state, tables->sites.ancestral_state_length,
            NPY_INT8},
        {"ancestral_state_offset",
            (void *) tables->sites.ancestral_state_offset, tables->sites.num_rows + 1,
            NPY_UINT32},
        {"metadata",
            (void *) tables->sites.metadata, tables->sites.metadata_length, NPY_INT8},
        {"metadata_offset",
            (void *) tables->sites.metadata_offset, tables->sites.num_rows + 1, NPY_UINT32},
        {NULL},
    };

    struct table_col mutation_cols[] = {
        {"site",
            (void *) tables->mutations.site, tables->mutations.num_rows, NPY_INT32},
        {"node",
            (void *) tables->mutations.node, tables->mutations.num_rows, NPY_INT32},
        {"time",
            (void *) tables->mutations.time, tables->mutations.num_rows, NPY_FLOAT64},
        {"parent",
            (void *) tables->mutations.parent, tables->mutations.num_rows, NPY_INT32},
        {"derived_state",
            (void *) tables->mutations.derived_state,
            tables->mutations.derived_state_length, NPY_INT8},
        {"derived_state_offset",
            (void *) tables->mutations.derived_state_offset,
            tables->mutations.num_rows + 1, NPY_UINT32},
        {"metadata",
            (void *) tables->mutations.metadata,
            tables->mutations.metadata_length, NPY_INT8},
        {"metadata_offset",
            (void *) tables->mutations.metadata_offset,
            tables->mutations.num_rows + 1, NPY_UINT32},
        {NULL},
    };

    struct table_col population_cols[] = {
        {"metadata", (void *) tables->populations.metadata,
            tables->populations.metadata_length, NPY_INT8},
        {"metadata_offset", (void *) tables->populations.metadata_offset,
            tables->populations.num_rows+ 1, NPY_UINT32},
        {NULL},
    };

    struct table_col provenance_cols[] = {
        {"timestamp", (void *) tables->provenances.timestamp,
            tables->provenances.timestamp_length, NPY_INT8},
        {"timestamp_offset", (void *) tables->provenances.timestamp_offset,
            tables->provenances.num_rows+ 1, NPY_UINT32},
        {"record", (void *) tables->provenances.record,
            tables->provenances.record_length, NPY_INT8},
        {"record_offset", (void *) tables->provenances.record_offset,
            tables->provenances.num_rows + 1, NPY_UINT32},
        {NULL},
    };

    struct table_desc table_descs[] = {
        {"individuals", individual_cols,
            tables->individuals.metadata_schema, tables->individuals.metadata_schema_length},
        {"nodes", node_cols,
            tables->nodes.metadata_schema, tables->nodes.metadata_schema_length},
        {"edges", edge_cols,
            tables->edges.metadata_schema, tables->edges.metadata_schema_length},
        {"migrations", migration_cols,
            tables->migrations.metadata_schema, tables->migrations.metadata_schema_length},
        {"sites", site_cols,
            tables->sites.metadata_schema, tables->sites.metadata_schema_length},
        {"mutations", mutation_cols,
            tables->mutations.metadata_schema, tables->mutations.metadata_schema_length},
        {"populations", population_cols,
            tables->populations.metadata_schema, tables->populations.metadata_schema_length},
        {"provenances", provenance_cols, NULL, 0},
    };

    for (j = 0; j < sizeof(table_descs) / sizeof(*table_descs); j++) {
        table_dict = PyDict_New();
        if (table_dict == NULL) {
            goto out;
        }
        col = table_descs[j].cols;
        while (col->name != NULL) {
            array = PyArray_SimpleNewFromData(1, &col->num_rows, col->type, col->data);
            if (array == NULL) {
                goto out;
            }
            if (PyDict_SetItemString(table_dict, col->name, array) != 0) {
                goto out;
            }
            Py_DECREF(array);
            array = NULL;
            col++;
        }
        if (table_descs[j].metadata_schema_length > 0) {
            array = make_Py_Unicode_FromStringAndLength(table_descs[j].metadata_schema,
                table_descs[j].metadata_schema_length);
            if (array == NULL) {
                goto out;
            }
            if (PyDict_SetItemString(table_dict, "metadata_schema", array) != 0) {
                goto out;
            }
            Py_DECREF(array);
            array = NULL;
        }

        if (PyDict_SetItemString(dict, table_descs[j].name, table_dict) != 0) {
            goto out;
        }
        Py_DECREF(table_dict);
        table_dict = NULL;
    }
    ret = 0;
out:
    Py_XDECREF(array);
    Py_XDECREF(table_dict);
    return ret;
}

/* Returns a dictionary encoding of the specified table collection */
static PyObject*
dump_tables_dict(tsk_table_collection_t *tables)
{
    PyObject *ret = NULL;
    PyObject *dict = NULL;
    PyObject *val = NULL;
    int err;

    dict = PyDict_New();
    if (dict == NULL) {
        goto out;
    }

    /* Dict representation version */
    val = Py_BuildValue("ll", 1, 1);
    if (val == NULL) {
        goto out;
    }
    if (PyDict_SetItemString(dict, "encoding_version", val) != 0) {
        goto out;
    }
    Py_DECREF(val);
    val = NULL;

    val = Py_BuildValue("d", tables->sequence_length);
    if (val == NULL) {
        goto out;
    }
    if (PyDict_SetItemString(dict, "sequence_length", val) != 0) {
        goto out;
    }
    Py_DECREF(val);
    val = NULL;

    if (tables->metadata_schema_length > 0) {
        val = make_Py_Unicode_FromStringAndLength(
            tables->metadata_schema, tables->metadata_schema_length);
        if (val == NULL) {
            goto out;
        }
        if (PyDict_SetItemString(dict, "metadata_schema", val) != 0) {
            goto out;
        }
        Py_DECREF(val);
        val = NULL;
    }

    if (tables->metadata_length > 0) {
        val = PyBytes_FromStringAndSize(tables->metadata, tables->metadata_length);
        if (val == NULL) {
            goto out;
        }
        if (PyDict_SetItemString(dict, "metadata", val) != 0) {
            goto out;
        }
        Py_DECREF(val);
        val = NULL;
    }

    err = write_table_arrays(tables, dict);
    if (err != 0) {
        goto out;
    }
    ret = dict;
    dict = NULL;
out:
    Py_XDECREF(dict);
    Py_XDECREF(val);
    return ret;
}

/*===================================================================
 * LightweightTableCollection
 *===================================================================
 */

static int
LightweightTableCollection_check_state(LightweightTableCollection *self)
{
    int ret = 0;
    if (self->tables == NULL) {
        PyErr_SetString(PyExc_SystemError, "LightweightTableCollection not initialised");
        ret = -1;
    }
    return ret;
}

static void
LightweightTableCollection_dealloc(LightweightTableCollection* self)
{
    if (self->tables != NULL) {
        tsk_table_collection_free(self->tables);
        PyMem_Free(self->tables);
        self->tables = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
LightweightTableCollection_init(LightweightTableCollection *self,
        PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"sequence_length", NULL};
    double sequence_length = -1;

    self->tables = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|d", kwlist, &sequence_length)) {
        goto out;
    }
    self->tables = PyMem_Malloc(sizeof(*self->tables));
    if (self->tables == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tsk_table_collection_init(self->tables, 0);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    self->tables->sequence_length = sequence_length;
    ret = 0;
out:
    return ret;
}

static PyObject *
LightweightTableCollection_asdict(LightweightTableCollection *self)
{
    PyObject *ret = NULL;

    if (LightweightTableCollection_check_state(self) != 0) {
        goto out;
    }
    ret = dump_tables_dict(self->tables);
out:
    return ret;
}

static PyObject *
LightweightTableCollection_fromdict(LightweightTableCollection *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    PyObject *dict = NULL;

    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) {
        goto out;
    }
    err = parse_table_collection_dict(self->tables, dict);
    if (err != 0) {
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyMemberDef LightweightTableCollection_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef LightweightTableCollection_methods[] = {
    {"asdict", (PyCFunction) LightweightTableCollection_asdict,
        METH_NOARGS, "Returns the tables encoded as a dictionary."},
    {"fromdict", (PyCFunction) LightweightTableCollection_fromdict,
        METH_VARARGS, "Populates the internal tables using the specified dictionary."},
    {NULL}  /* Sentinel */
};

static PyTypeObject LightweightTableCollectionType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_msprime.LightweightTableCollection",             /* tp_name */
    sizeof(LightweightTableCollection),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)LightweightTableCollection_dealloc, /* tp_dealloc */
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
    "LightweightTableCollection objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    LightweightTableCollection_methods,             /* tp_methods */
    LightweightTableCollection_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)LightweightTableCollection_init,      /* tp_init */
};

/***********************************************
 * End code copied from _tskitmodule.
 ***********************************************/

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

static PyObject *
RandomGenerator_flat(RandomGenerator *self, PyObject *args)
{
    PyObject *ret = NULL;
    double a, b;

    if (RandomGenerator_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "dd", &a, &b)) {
        goto out;
    }
    ret = Py_BuildValue("d", gsl_ran_flat(self->rng, a, b));
out:
    return ret;
}

static PyObject *
RandomGenerator_poisson(RandomGenerator *self, PyObject *args)
{
    PyObject *ret = NULL;
    double mu;

    if (RandomGenerator_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "d", &mu)) {
        goto out;
    }
    ret = Py_BuildValue("I", gsl_ran_poisson(self->rng, mu));
out:
    return ret;
}

static PyObject *
RandomGenerator_uniform_int(RandomGenerator *self, PyObject *args)
{
    PyObject *ret = NULL;
    unsigned long n;

    if (RandomGenerator_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "k", &n)) {
        goto out;
    }
    ret = Py_BuildValue("k", gsl_rng_uniform_int(self->rng, n));
out:
    return ret;
}

static PyMemberDef RandomGenerator_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef RandomGenerator_methods[] = {
    {"get_seed", (PyCFunction) RandomGenerator_get_seed,
        METH_NOARGS, "Returns the random seed for this generator."},
    {"flat", (PyCFunction) RandomGenerator_flat,
        METH_VARARGS, "Interface for gsl_ran_flat"},
    {"poisson", (PyCFunction) RandomGenerator_poisson,
        METH_VARARGS, "Interface for gsl_ran_poisson"},
    {"uniform_int", (PyCFunction) RandomGenerator_uniform_int,
        METH_VARARGS, "Interface for gsl_rng_uniform_int"},
    {NULL}  /* Sentinel */
};

static PyTypeObject RandomGeneratorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_msprime.RandomGenerator",
    .tp_basicsize = sizeof(RandomGenerator),
    .tp_dealloc = (destructor)RandomGenerator_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = "RandomGenerator objects",
    .tp_methods = RandomGenerator_methods,
    .tp_members = RandomGenerator_members,
    .tp_init = (initproc)RandomGenerator_init,
    .tp_new = PyType_GenericNew,
};

/*===================================================================
 * Base mutation model
 *===================================================================
 */

static PyTypeObject BaseMutationModelType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_msprime.BaseMutationModel",
    .tp_basicsize = sizeof(BaseMutationModel),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = "BaseMutationModel objects",
    .tp_new = PyType_GenericNew,
};

/*===================================================================
 * Matrix mutation model
 *===================================================================
 */

static int
MatrixMutationModel_check_state(MatrixMutationModel *self)
{
    int ret = 0;
    if (self->mutation_model == NULL) {
        PyErr_SetString(PyExc_SystemError, "MatrixMutationModel not initialised");
        ret = -1;
    }
    return ret;
}

static void
MatrixMutationModel_dealloc(MatrixMutationModel* self)
{
    if (self->mutation_model != NULL) {
        mutation_model_free(self->mutation_model);
        PyMem_Free(self->mutation_model);
        self->mutation_model = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
MatrixMutationModel_init(MatrixMutationModel *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"alleles", "root_distribution", "transition_matrix", NULL};
    Py_ssize_t j, num_alleles;
    PyObject *py_alleles = NULL;
    PyArrayObject *root_distribution_array = NULL;
    PyObject *py_transition_matrix = NULL;
    PyArrayObject *transition_matrix_array = NULL;
    char **alleles = NULL;
    size_t *allele_length = NULL;
    PyObject *item;
    npy_intp *shape;
    Py_ssize_t len;

    self->mutation_model = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O&O", kwlist,
            &PyList_Type, &py_alleles,
            double_PyArray_converter, &root_distribution_array,
            &py_transition_matrix)) {
        goto out;
    }
    num_alleles = PyList_Size(py_alleles);
    shape = PyArray_DIMS(root_distribution_array);
    if (num_alleles != shape[0]) {
        PyErr_SetString(PyExc_ValueError,
            "root distribution must have num_alleles elements");
        goto out;
    }
    transition_matrix_array = (PyArrayObject *) PyArray_FROMANY(
            py_transition_matrix, NPY_FLOAT64, 2, 2, NPY_ARRAY_IN_ARRAY);
    if (transition_matrix_array == NULL) {
        goto out;
    }
    shape = PyArray_DIMS(transition_matrix_array);
    if (shape[0] != shape[1]) {
        PyErr_SetString(PyExc_ValueError, "Square matrix required");
        goto out;
    }
    if (shape[0] != num_alleles) {
        PyErr_SetString(PyExc_ValueError,
            "transition matrix must be a square matrix with num_alleles rows");
        goto out;
    }

    /* Note: it's important we zero out mutation_model here because
     * we can error before we can mutation_model_alloc, leaving the
     * object in an uninitialised state */
    self->mutation_model = PyMem_Calloc(1, sizeof(*self->mutation_model));
    alleles = PyMem_Malloc(num_alleles * sizeof(*alleles));
    allele_length = PyMem_Malloc(num_alleles * sizeof(*allele_length));
    if (self->mutation_model == NULL || alleles == NULL || allele_length == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    /* This is a shortcut as it depends on using the underlying memory
     * of the unicode object. Because we don't execute any Python code
     * before calling mutation_model_alloc and we take a copies of
     * the allele strings in there, this should be safe. */
    for (j = 0; j < num_alleles; j++) {
        item = PyList_GetItem(py_alleles, j);
        if (!PyUnicode_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "alleles must be unicode strings");
            goto out;
        }
        alleles[j] = (char *) (intptr_t) PyUnicode_AsUTF8AndSize(item, &len);
        if (alleles[j] == NULL) {
            goto out;
        }
        allele_length[j] = (size_t) len;
    }
    err = matrix_mutation_model_alloc(self->mutation_model,
            num_alleles, alleles, allele_length,
            PyArray_DATA(root_distribution_array),
            PyArray_DATA(transition_matrix_array));
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    PyMem_Free(alleles);
    PyMem_Free(allele_length);
    Py_XDECREF(root_distribution_array);
    Py_XDECREF(transition_matrix_array);
    return ret;
}

static PyObject *
MatrixMutationModel_get_alleles(MatrixMutationModel *self, void *closure)
{
    PyObject *ret = NULL;
    size_t j;
    PyObject *item;
    mutation_matrix_t *params;
    size_t size;
    PyObject *list = NULL;

    if (MatrixMutationModel_check_state(self) != 0) {
        goto out;
    }
    params = &self->mutation_model->params.mutation_matrix;
    size = params->num_alleles;
    list = PyList_New(size);
    if (list == NULL) {
        goto out;
    }
    for (j = 0; j < size; j++) {
        item = PyUnicode_FromStringAndSize(params->alleles[j], params->allele_length[j]);
        if (item == NULL) {
            goto out;
        }
        PyList_SET_ITEM(list, j, item);
    }
    ret = list;
    list = NULL;
out:
    Py_XDECREF(list);
    return ret;
}

static PyObject *
MatrixMutationModel_get_root_distribution(MatrixMutationModel *self, void *closure)
{
    PyObject *ret = NULL;
    PyArrayObject *array;
    mutation_matrix_t *params = NULL;
    size_t size;
    npy_intp dims;

    if (MatrixMutationModel_check_state(self) != 0) {
        goto out;
    }
    params = &self->mutation_model->params.mutation_matrix;
    size = params->num_alleles;
    dims = (npy_intp) size;
    array = (PyArrayObject *) PyArray_EMPTY(1, &dims, NPY_FLOAT64, 0);
    if (array == NULL) {
        goto out;
    }
    memcpy(PyArray_DATA(array), params->root_distribution,
            size * sizeof(double));
    ret = (PyObject *) array;
out:
    return ret;
}

static PyObject *
MatrixMutationModel_get_transition_matrix(MatrixMutationModel *self, void *closure)
{
    PyObject *ret = NULL;
    PyArrayObject *array;
    mutation_matrix_t *params;
    size_t size;
    npy_intp dims[2];

    if (MatrixMutationModel_check_state(self) != 0) {
        goto out;
    }
    params = &self->mutation_model->params.mutation_matrix;
    size = params->num_alleles;
    dims[0] = size;
    dims[1] = size;
    array = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 0);
    if (array == NULL) {
        goto out;
    }
    memcpy(PyArray_DATA(array), params->transition_matrix,
            size * size * sizeof(double));
    ret = (PyObject *) array;
out:
    return ret;
}

static PyGetSetDef MatrixMutationModel_getsetters[] = {
    {"alleles", (getter) MatrixMutationModel_get_alleles, NULL,
        "A copy of the alleles list"},
    {"root_distribution", (getter) MatrixMutationModel_get_root_distribution, NULL,
        "A copy of the root_distribution array"},
    {"transition_matrix", (getter) MatrixMutationModel_get_transition_matrix, NULL,
        "A copy of the transition_matrix array"},
    {NULL}  /* Sentinel */
};

static PyTypeObject MatrixMutationModelType = {
    .tp_name = "_msprime.MatrixMutationModel",
    .tp_basicsize = sizeof(MatrixMutationModel),
    .tp_dealloc = (destructor)MatrixMutationModel_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_doc = "MatrixMutationModel objects",
    .tp_getset = MatrixMutationModel_getsetters,
    .tp_init = (initproc)MatrixMutationModel_init,
    .tp_new = PyType_GenericNew,
};

/*===================================================================
 * Slim mutation model
 *===================================================================
 */

static int
SLiMMutationModel_check_state(SLiMMutationModel *self)
{
    int ret = 0;
    if (self->mutation_model == NULL) {
        PyErr_SetString(PyExc_SystemError, "SLiMMutationModel not initialised");
        ret = -1;
    }
    return ret;
}

static void
SLiMMutationModel_dealloc(SLiMMutationModel* self)
{
    if (self->mutation_model != NULL) {
        mutation_model_free(self->mutation_model);
        PyMem_Free(self->mutation_model);
        self->mutation_model = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
SLiMMutationModel_init(SLiMMutationModel *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"type", "next_id", "block_size", NULL};
    long type;
    long long next_id = 0;
    Py_ssize_t block_size = 0;

    self->mutation_model = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|Ln", kwlist,
            &type, &next_id, &block_size)) {
        goto out;
    }

    /* Note: it's important we zero out mutation_model here because
     * we can error before we can mutation_model_alloc, leaving the
     * object in an uninitialised state */
    self->mutation_model = PyMem_Calloc(1, sizeof(*self->mutation_model));
    if (self->mutation_model == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = slim_mutation_model_alloc(self->mutation_model,
        (int32_t) type, (int64_t) next_id, (size_t) block_size);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
SLiMMutationModel_get_type(SLiMMutationModel *self, void *closure)
{
    slim_mutator_t *params;
    PyObject *ret = NULL;

    if (SLiMMutationModel_check_state(self) != 0) {
        goto out;
    }
    params = &self->mutation_model->params.slim_mutator;
    ret = Py_BuildValue("l", (long) params->mutation_type_id);
out:
    return ret;
}

static PyObject *
SLiMMutationModel_get_next_id(SLiMMutationModel *self, void *closure)
{
    slim_mutator_t *params;
    PyObject *ret = NULL;

    if (SLiMMutationModel_check_state(self) != 0) {
        goto out;
    }
    params = &self->mutation_model->params.slim_mutator;
    ret = Py_BuildValue("L", params->next_mutation_id);
out:
    return ret;
}

static PyGetSetDef SLiMMutationModel_getsetters[] = {
    {"type", (getter) SLiMMutationModel_get_type, NULL,
        "Return the mutation type"},
    {"next_id", (getter) SLiMMutationModel_get_next_id, NULL,
        "Return the next mutation id"},
    {NULL}  /* Sentinel */
};

static PyTypeObject SLiMMutationModelType = {
    .tp_name = "_msprime.SLiMMutationModel",
    .tp_basicsize = sizeof(SLiMMutationModel),
    .tp_dealloc = (destructor)SLiMMutationModel_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_doc = "SLiMMutationModel objects",
    .tp_getset = SLiMMutationModel_getsetters,
    .tp_init = (initproc)SLiMMutationModel_init,
    .tp_new = PyType_GenericNew,
};

/*===================================================================
 * Infinite alleles mutation model
 *===================================================================
 */

static int
InfiniteAllelesMutationModel_check_state(InfiniteAllelesMutationModel *self)
{
    int ret = 0;
    if (self->mutation_model == NULL) {
        PyErr_SetString(PyExc_SystemError, "InfiniteAllelesMutationModel not initialised");
        ret = -1;
    }
    return ret;
}

static void
InfiniteAllelesMutationModel_dealloc(InfiniteAllelesMutationModel* self)
{
    if (self->mutation_model != NULL) {
        mutation_model_free(self->mutation_model);
        PyMem_Free(self->mutation_model);
        self->mutation_model = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
InfiniteAllelesMutationModel_init(InfiniteAllelesMutationModel *self,
        PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"start_allele", NULL};
    unsigned long long start_allele = 0;

    self->mutation_model = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|K", kwlist, &start_allele)) {
        goto out;
    }

    /* Note: it's important we zero out mutation_model here because
     * we can error before we can mutation_model_alloc, leaving the
     * object in an uninitialised state */
    self->mutation_model = PyMem_Calloc(1, sizeof(*self->mutation_model));
    if (self->mutation_model == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = infinite_alleles_mutation_model_alloc(self->mutation_model,
        (uint64_t) start_allele, 0);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
InfiniteAllelesMutationModel_get_start_allele(InfiniteAllelesMutationModel *self,
        void *closure)
{
    PyObject *ret = NULL;
    infinite_alleles_t *params;

    if (InfiniteAllelesMutationModel_check_state(self) != 0) {
        goto out;
    }
    params = &self->mutation_model->params.infinite_alleles;
    ret = Py_BuildValue("K", (unsigned long long) params->start_allele);
out:
    return ret;
}

static PyObject *
InfiniteAllelesMutationModel_get_next_allele(InfiniteAllelesMutationModel *self,
        void *closure)
{
    PyObject *ret = NULL;
    infinite_alleles_t *params;

    if (InfiniteAllelesMutationModel_check_state(self) != 0) {
        goto out;
    }
    params = &self->mutation_model->params.infinite_alleles;
    ret = Py_BuildValue("K", (unsigned long long) params->next_allele);
out:
    return ret;
}

static PyGetSetDef InfiniteAllelesMutationModel_getsetters[] = {
    {"start_allele", (getter) InfiniteAllelesMutationModel_get_start_allele, NULL,
        "Returns the initial allele"},
    {"next_allele", (getter) InfiniteAllelesMutationModel_get_next_allele, NULL,
        "Return the next allele"},
    {NULL}  /* Sentinel */
};

static PyTypeObject InfiniteAllelesMutationModelType = {
    .tp_name = "_msprime.InfiniteAllelesMutationModel",
    .tp_basicsize = sizeof(InfiniteAllelesMutationModel),
    .tp_dealloc = (destructor)InfiniteAllelesMutationModel_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_doc = "InfiniteAllelesMutationModel objects",
    .tp_getset = InfiniteAllelesMutationModel_getsetters,
    .tp_init = (initproc)InfiniteAllelesMutationModel_init,
    .tp_new = PyType_GenericNew,
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
Simulator_parse_population_configuration(Simulator *self, PyObject *py_pop_config)
{
    int ret = -1;
    Py_ssize_t j;
    double initial_size, growth_rate;
    int err;
    PyObject *item, *value;

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
            handle_input_error("population configuration", err);
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}

static int
Simulator_parse_migration_matrix(Simulator *self, PyObject *py_migration_matrix)
{
    int ret = -1;
    int err;
    npy_intp *shape;
    PyArrayObject *migration_matrix_array = NULL;
    size_t num_populations = msp_get_num_populations(self->sim);
    const char *err_msg =
        "migration matrix must be a N x N square matrix encoded "
        "as a list-of-lists or numpy array, where N is the number of populations "
        "defined in the population_configurations. The diagonal "
        "elements of this matrix must be zero. For example, a "
        "valid matrix for a 3 population system is "
        "[[0, 1, 1], [1, 0, 1], [1, 1, 0]]";

    migration_matrix_array = (PyArrayObject *) PyArray_FROMANY(
            py_migration_matrix, NPY_FLOAT64, 2, 2, NPY_ARRAY_IN_ARRAY);
    if (migration_matrix_array == NULL) {
        goto out;
    }
    shape = PyArray_DIMS(migration_matrix_array);
    if (shape[0] != shape[1]) {
        PyErr_SetString(PyExc_ValueError, err_msg);
        goto out;
    }
    if (shape[0] != (npy_intp) num_populations) {
        PyErr_SetString(PyExc_ValueError, err_msg);
        goto out;
    }
    err = msp_set_migration_matrix(self->sim,
            num_populations * num_populations,
            PyArray_DATA(migration_matrix_array));
    if (err != 0) {
        handle_input_error("migration matrix", err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(migration_matrix_array);
    return ret;
}

static int
Simulator_parse_sweep_genic_selection_model(Simulator *self, PyObject *py_model)
{
    int ret = -1;
    int err;
    double position, start_frequency, end_frequency, alpha, dt;
    PyObject *value;

    value = get_dict_number(py_model, "position");
    if (value == NULL) {
        goto out;
    }
    position = PyFloat_AsDouble(value);

    value = get_dict_number(py_model, "start_frequency");
    if (value == NULL) {
        goto out;
    }
    start_frequency = PyFloat_AsDouble(value);

    value = get_dict_number(py_model, "end_frequency");
    if (value == NULL) {
        goto out;
    }
    end_frequency = PyFloat_AsDouble(value);

    value = get_dict_number(py_model, "alpha");
    if (value == NULL) {
        goto out;
    }
    alpha = PyFloat_AsDouble(value);

    value = get_dict_number(py_model, "dt");
    if (value == NULL) {
        goto out;
    }
    dt = PyFloat_AsDouble(value);

    err = msp_set_simulation_model_sweep_genic_selection(self->sim,
            position, start_frequency, end_frequency, alpha, dt);
    if (err != 0) {
        handle_input_error("sweep genic selection", err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static int
Simulator_parse_simulation_model(Simulator *self, PyObject *py_model)
{
    int ret = -1;
    int err = 0;
    PyObject *py_name = NULL;
    PyObject *hudson_s = NULL;
    PyObject *smc_s = NULL;
    PyObject *smc_prime_s = NULL;
    PyObject *dtwf_s = NULL;
    PyObject *wf_ped_s = NULL;
    PyObject *dirac_s = NULL;
    PyObject *beta_s = NULL;
    PyObject *sweep_genic_selection_s = NULL;
    PyObject *value;
    int is_hudson, is_dtwf, is_smc, is_smc_prime, is_dirac, is_beta, is_sweep_genic_selection;
    int is_wf_ped;
    double psi, c, alpha, truncation_point;

    hudson_s = Py_BuildValue("s", "hudson");
    if (hudson_s == NULL) {
        goto out;
    }
    dtwf_s = Py_BuildValue("s", "dtwf");
    if (dtwf_s == NULL) {
        goto out;
    }
    wf_ped_s = Py_BuildValue("s", "wf_ped");
    if (wf_ped_s == NULL) {
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
    sweep_genic_selection_s = Py_BuildValue("s", "sweep_genic_selection");
    if (sweep_genic_selection_s == NULL) {
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
        err = msp_set_simulation_model_hudson(self->sim);
    }

    is_dtwf = PyObject_RichCompareBool(py_name, dtwf_s, Py_EQ);
    if (is_dtwf == -1) {
        goto out;
    }
    if (is_dtwf) {
        err = msp_set_simulation_model_dtwf(self->sim);
    }
    is_wf_ped = PyObject_RichCompareBool(py_name, wf_ped_s, Py_EQ);
    if (is_wf_ped == -1) {
        goto out;
    }
    if (is_wf_ped) {
        err = msp_set_simulation_model_wf_ped(self->sim);
    }

    is_smc = PyObject_RichCompareBool(py_name, smc_s, Py_EQ);
    if (is_smc == -1) {
        goto out;
    }
    if (is_smc) {
        err = msp_set_simulation_model_smc(self->sim);
    }

    is_smc_prime = PyObject_RichCompareBool(py_name, smc_prime_s, Py_EQ);
    if (is_smc_prime == -1) {
        goto out;
    }
    if (is_smc_prime) {
        err = msp_set_simulation_model_smc_prime(self->sim);
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
        value = get_dict_number(py_model, "c");
        if (value == NULL) {
            goto out;
        }
        c = PyFloat_AsDouble(value);
        if (psi <= 0 || psi >= 1.0) {
            PyErr_SetString(PyExc_ValueError, "Must have 0 < psi < 1");
            goto out;
        }
        if (c < 0) {
            PyErr_SetString(PyExc_ValueError, "c >= 0");
            goto out;
        }
        err = msp_set_simulation_model_dirac(self->sim, psi, c);
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
        err = msp_set_simulation_model_beta(self->sim,
                alpha, truncation_point);
    }

    is_sweep_genic_selection = PyObject_RichCompareBool(py_name,
            sweep_genic_selection_s, Py_EQ);
    if (is_sweep_genic_selection == -1) {
        goto out;
    }
    if (is_sweep_genic_selection) {
        ret = Simulator_parse_sweep_genic_selection_model(self, py_model);
        if (ret != 0) {
            goto out;
        }
    }

    if (! (is_hudson || is_dtwf || is_smc || is_smc_prime || is_dirac
                || is_beta || is_sweep_genic_selection || is_wf_ped)) {
        PyErr_SetString(PyExc_ValueError, "Unknown simulation model");
        goto out;
    }
    if (err != 0) {
        handle_input_error("simulation model", err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(hudson_s);
    Py_XDECREF(dtwf_s);
    Py_XDECREF(wf_ped_s);
    Py_XDECREF(smc_s);
    Py_XDECREF(smc_prime_s);
    Py_XDECREF(beta_s);
    Py_XDECREF(dirac_s);
    Py_XDECREF(sweep_genic_selection_s);
    return ret;
}

static int
Simulator_parse_recombination_map(Simulator *self, PyObject *py_recomb_map)
{
    int ret = -1;
    int err = 0;
    PyArrayObject *position_array = NULL;
    PyArrayObject *rate_array = NULL;
    size_t size;

    err = parse_rate_map(py_recomb_map, &size, &position_array, &rate_array);
    if (err != 0) {
        goto out;
    }
    err = msp_set_recombination_map(self->sim,
            (size_t) size,
            PyArray_DATA(position_array),
            PyArray_DATA(rate_array));
    if (err != 0) {
        handle_input_error("recombination map", err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(position_array);
    Py_XDECREF(rate_array);
    return ret;
}

static int
Simulator_parse_demographic_events(Simulator *self, PyObject *py_events)
{
    int ret = -1;
    Py_ssize_t j;
    double time, initial_size, growth_rate, migration_rate, proportion,
           strength;
    int err, population_id, source, destination;
    int is_population_parameter_change, is_migration_rate_change, is_mass_migration,
        is_simple_bottleneck, is_instantaneous_bottleneck, is_census_event;
    PyObject *item, *value, *type;
    PyObject *population_parameter_change_s = NULL;
    PyObject *migration_rate_change_s = NULL;
    PyObject *mass_migration_s = NULL;
    PyObject *simple_bottleneck_s = NULL;
    PyObject *instantaneous_bottleneck_s = NULL;
    PyObject *census_event_s = NULL;
    PyObject *initial_size_s = NULL;
    PyObject *growth_rate_s = NULL;

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
    census_event_s = Py_BuildValue("s", "census_event");
    if (census_event_s == NULL) {
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
        is_census_event = PyObject_RichCompareBool(
                type, census_event_s, Py_EQ);
        if (is_census_event == -1) {
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
            value = get_dict_number(item, "population");
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
            value = get_dict_number(item, "source");
            if (value == NULL) {
                goto out;
            }
            source = (int) PyLong_AsLong(value);
            value = get_dict_number(item, "dest");
            if (value == NULL) {
                goto out;
            }
            destination = (int) PyLong_AsLong(value);
            err = msp_add_migration_rate_change(self->sim, time, source, destination,
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
            value = get_dict_number(item, "dest");
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
            value = get_dict_number(item, "population");
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
            value = get_dict_number(item, "population");
            if (value == NULL) {
                goto out;
            }
            population_id = (int) PyLong_AsLong(value);
            err = msp_add_instantaneous_bottleneck(self->sim, time, population_id,
                    strength);
        } else if (is_census_event) {
            err = msp_add_census_event(self->sim, time);
        } else {
            PyErr_Format(PyExc_ValueError, "Unknown demographic event type");
            goto out;
        }
        if (err != 0) {
            PyErr_Format(MsprimeInputError,
                    "Input error in demographic_events[%d]: %s", j,
                    msp_strerror(err));
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
    Py_XDECREF(census_event_s);
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
    Py_XDECREF(self->tables);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
Simulator_init(Simulator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int sim_ret;
    static char *kwlist[] = {
        "tables", "random_generator", "recombination_map",
        "population_configuration", "migration_matrix",
        "demographic_events", "model", "avl_node_block_size", "segment_block_size",
        "node_mapping_block_size", "store_migrations", "start_time",
        "store_full_arg", "num_labels", "gene_conversion_rate",
        "gene_conversion_track_length", "discrete_genome",
        "ploidy", NULL};
    PyObject *migration_matrix = NULL;
    PyObject *population_configuration = NULL;
    PyObject *demographic_events = NULL;
    PyObject *py_model = NULL;
    LightweightTableCollection *tables = NULL;
    RandomGenerator *random_generator = NULL;
    PyObject *recombination_map = NULL;
    /* parameter defaults */
    Py_ssize_t avl_node_block_size = 10;
    Py_ssize_t segment_block_size = 10;
    Py_ssize_t node_mapping_block_size = 10;
    Py_ssize_t num_labels = 1;
    Py_ssize_t num_populations = 1;
    int store_migrations = false;
    int store_full_arg = false;
    int discrete_genome = true;
    double start_time = -1;
    double gene_conversion_rate = 0;
    double gene_conversion_track_length = 1.0;
    int ploidy = 2;

    self->sim = NULL;
    self->random_generator = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds,
            "O!O!|O!O!OO!O!nnnidinddii", kwlist,
            &LightweightTableCollectionType, &tables,
            &RandomGeneratorType, &random_generator,
            /* optional */
            &PyDict_Type, &recombination_map,
            &PyList_Type, &population_configuration,
            &migration_matrix,
            &PyList_Type, &demographic_events,
            &PyDict_Type, &py_model,
            &avl_node_block_size, &segment_block_size,
            &node_mapping_block_size, &store_migrations, &start_time,
            &store_full_arg, &num_labels,
            &gene_conversion_rate, &gene_conversion_track_length,
            &discrete_genome, &ploidy)) {
        goto out;
    }
    self->random_generator = random_generator;
    self->tables = tables;
    Py_INCREF(self->random_generator);
    Py_INCREF(self->tables);

    if (RandomGenerator_check_state(self->random_generator) != 0) {
        goto out;
    }
    if (population_configuration != NULL) {
        num_populations = PyList_Size(population_configuration);
        if (num_populations == 0) {
            PyErr_SetString(PyExc_ValueError, "Empty population configuration");
            goto out;
        }
    }
    self->sim = PyMem_Malloc(sizeof(msp_t));
    if (self->sim == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    sim_ret = msp_alloc(self->sim, tables->tables, self->random_generator->rng);
    if (sim_ret != 0) {
        handle_input_error("simulator alloc", sim_ret);
        goto out;
    }
    if (py_model != NULL) {
        if (Simulator_parse_simulation_model(self, py_model) != 0) {
            goto out;
        }
    }
    if (start_time >= 0) {
        sim_ret = msp_set_start_time(self->sim, start_time);
        if (sim_ret != 0) {
            handle_input_error("start time", sim_ret);
            goto out;
        }
    }
    sim_ret = msp_set_store_migrations(self->sim, (bool) store_migrations);
    if (sim_ret != 0) {
        handle_input_error("store migrations", sim_ret);
        goto out;
    }
    sim_ret = msp_set_avl_node_block_size(self->sim,
            (size_t) avl_node_block_size);
    if (sim_ret != 0) {
        handle_input_error("avl_node_block_size", sim_ret);
        goto out;
    }
    sim_ret = msp_set_segment_block_size(self->sim,
            (size_t) segment_block_size);
    if (sim_ret != 0) {
        handle_input_error("segment_block_size", sim_ret);
        goto out;
    }
    sim_ret = msp_set_node_mapping_block_size(self->sim,
            (size_t) node_mapping_block_size);
    if (sim_ret != 0) {
        handle_input_error("node_mapping_block_size", sim_ret);
        goto out;
    }
    if (gene_conversion_rate != 0) {
        sim_ret = msp_set_gene_conversion_rate(self->sim, gene_conversion_rate);
        if (sim_ret != 0) {
            handle_input_error("set_gene_conversion_rate", sim_ret);
            goto out;
        }
        sim_ret = msp_set_gene_conversion_track_length(self->sim,
                gene_conversion_track_length);
        if (sim_ret != 0) {
            handle_input_error("set_gene_conversion_track_length", sim_ret);
            goto out;
        }
    }
    if (recombination_map != NULL) {
        if (Simulator_parse_recombination_map(self, recombination_map) != 0) {
            goto out;
        }
    }
    msp_set_discrete_genome(self->sim, discrete_genome);

    sim_ret = msp_set_ploidy(self->sim, ploidy);
    if (sim_ret != 0) {
        handle_input_error("set_ploidy", sim_ret);
        goto out;
    }

    sim_ret = msp_set_num_labels(self->sim, (size_t) num_labels);
    if (sim_ret != 0) {
        handle_input_error("set_num_labels", sim_ret);
        goto out;
    }
    if (population_configuration != NULL) {
        if (Simulator_parse_population_configuration(self, population_configuration) != 0) {
            goto out;
        }
        if (migration_matrix == NULL) {
            PyErr_SetString(PyExc_ValueError,
                "A migration matrix must be provided when a non-default "
                "population configuration is used.");
            goto out;
        }
        if (Simulator_parse_migration_matrix(self, migration_matrix) != 0) {
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
    msp_set_store_full_arg(self->sim, store_full_arg);

    sim_ret = msp_initialise(self->sim);
    if (sim_ret != 0) {
        handle_input_error("initialise", sim_ret);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
Simulator_get_model(Simulator *self, void *closure)
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
        value = Py_BuildValue("d", model->params.dirac_coalescent.c);
        if (value == NULL) {
            goto out;
        }
        if (PyDict_SetItemString(d, "c", value) != 0) {
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
    } else if (model->type == MSP_MODEL_SWEEP) {
        value = Py_BuildValue("d", model->params.sweep.position);
        if (value == NULL) {
            goto out;
        }
        if (PyDict_SetItemString(d, "locus", value) != 0) {
            goto out;
        }
        Py_DECREF(value);
        value = NULL;
        /* TODO fill in the parameters for the different types of trajectories. */
    }
    ret = d;
    d = NULL;
out:
    Py_XDECREF(d);
    Py_XDECREF(value);
    return ret;
}

static int
Simulator_set_model(Simulator *self, PyObject *args, void *closure)
{
    int ret = -1;
    PyObject *py_model = args;

    if (py_model == NULL) {
        /* deleting the model attribute isn't supported */
        goto out;
    }

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (Simulator_parse_simulation_model(self, py_model) != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
Simulator_get_discrete_genome(Simulator *self, void *closure)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("i", self->sim->discrete_genome);
out:
    return ret;
}

static PyObject *
Simulator_get_ploidy(Simulator *self, void *closure)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("i", self->sim->ploidy);
out:
    return ret;
}

static PyObject *
Simulator_get_store_migrations(Simulator *self, void *closure)
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
Simulator_get_num_populations(Simulator *self, void *closure)
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
Simulator_get_num_labels(Simulator *self, void *closure)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_num_labels(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_sequence_length(Simulator *self, void *closure)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sim->sequence_length);
out:
    return ret;
}

static PyObject *
Simulator_get_segment_block_size(Simulator  *self, void *closure)
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
Simulator_get_avl_node_block_size(Simulator  *self, void *closure)
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
Simulator_get_node_mapping_block_size(Simulator  *self, void *closure)
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
Simulator_get_time(Simulator  *self, void *closure)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", msp_get_time(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_ancestors(Simulator *self, void *closure)
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
Simulator_get_num_common_ancestor_events(Simulator  *self, void *closure)
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
Simulator_get_num_rejected_common_ancestor_events(Simulator  *self, void *closure)
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
Simulator_get_num_recombination_events(Simulator  *self, void *closure)
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
Simulator_get_num_gene_conversion_events(Simulator  *self, void *closure)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n",
        (Py_ssize_t) msp_get_num_gene_conversion_events(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_migration_events(Simulator  *self, void *closure)
{
    PyObject *ret = NULL;
    PyObject *arr = NULL;
    int err;
    npy_intp size[2];

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    size[0] = msp_get_num_populations(self->sim);
    size[1] = size[0];
    arr = PyArray_SimpleNew(2, size, NPY_UINTP);
    if (arr == NULL) {
        goto out;
    }
    err = msp_get_num_migration_events(self->sim, PyArray_DATA((PyArrayObject *) arr));
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = arr;
    arr = NULL;
out:
    Py_XDECREF(arr);
    return ret;
}

static PyObject *
Simulator_get_num_multiple_recombination_events(Simulator  *self, void *closure)
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
Simulator_get_num_avl_node_blocks(Simulator  *self, void *closure)
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
Simulator_get_num_node_mapping_blocks(Simulator  *self, void *closure)
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
Simulator_get_num_segment_blocks(Simulator  *self, void *closure)
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
Simulator_get_num_fenwick_rebuilds(Simulator  *self, void *closure)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sim->num_fenwick_rebuilds);
out:
    return ret;
}

static PyObject *
Simulator_get_num_breakpoints(Simulator  *self, void *closure)
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
Simulator_get_num_edges(Simulator  *self, void *closure)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_num_edges(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_nodes(Simulator  *self, void *closure)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) msp_get_num_nodes(self->sim));
out:
    return ret;
}

static PyObject *
Simulator_get_num_migrations(Simulator  *self, void *closure)
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
        t = Py_BuildValue("(d,d,I,I)", u->left, u->right, u->value,
                u->population);
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
Simulator_get_ancestors(Simulator *self, void *closure)
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
Simulator_get_breakpoints(Simulator *self, void *closure)
{
    PyObject *ret = NULL;
    PyObject *arr = NULL;
    npy_intp num_breakpoints;
    int err;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    num_breakpoints = msp_get_num_breakpoints(self->sim);
    arr = PyArray_SimpleNew(1, &num_breakpoints, NPY_UINTP);
    if (arr == NULL) {
        goto out;
    }
    err = msp_get_breakpoints(self->sim, PyArray_DATA((PyArrayObject *) arr));
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = arr;
    arr = NULL;
out:
    Py_XDECREF(arr);
    return ret;
}

static PyObject *
Simulator_get_recombination_map(Simulator *self, void *closure)
{
    PyObject *ret = NULL;
    rate_map_t rate_map;
    PyArrayObject *position = NULL;
    PyArrayObject *rate = NULL;
    npy_intp dims;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    rate_map = self->sim->recomb_map;
    dims = rate_map.size + 1;
    position = (PyArrayObject *) PyArray_SimpleNew(1, &dims, NPY_FLOAT64);
    dims = rate_map.size;
    rate = (PyArrayObject *) PyArray_SimpleNew(1, &dims, NPY_FLOAT64);
    if (position == NULL || rate == NULL) {
        goto out;
    }
    memcpy(PyArray_DATA(position), rate_map.position,
            (rate_map.size + 1) * (sizeof(*rate_map.position)));
    memcpy(PyArray_DATA(rate), rate_map.rate,
            (rate_map.size) * (sizeof(*rate_map.rate)));
    ret = Py_BuildValue("{s:O,s:O}",
        "position", position,
        "rate", rate);
out:
    Py_XDECREF(position);
    Py_XDECREF(rate);
    return ret;
}

static PyObject *
Simulator_get_migration_matrix(Simulator *self, void *closure)
{
    PyObject *ret = NULL;
    PyObject *arr = NULL;
    npy_intp N[2];
    int err;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    N[0] = msp_get_num_populations(self->sim);
    N[1] = N[0];
    arr = PyArray_SimpleNew(2, N, NPY_FLOAT64);
    if (arr == NULL) {
        goto out;
    }
    err = msp_get_migration_matrix(self->sim, PyArray_DATA((PyArrayObject *)arr));
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = arr;
    arr = NULL;
out:
    Py_XDECREF(arr);
    return ret;
}

static PyObject *
Simulator_get_population_configuration(Simulator *self, void *closure)
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
Simulator_get_random_generator(Simulator *self, void *closure)
{
    return Py_BuildValue("O", self->random_generator);
}

static PyObject *
Simulator_get_tables(Simulator *self, void *closure)
{
    return Py_BuildValue("O", self->tables);
}

static PyObject *
Simulator_run(Simulator *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    static char *kwlist[] = {"end_time", "max_events", NULL};
    int status;
    unsigned long max_events = UINT32_MAX;
    double end_time = DBL_MAX;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dk", kwlist,
                &end_time, &max_events)) {
        goto out;
    }
    if (end_time < 0) {
        PyErr_SetString(PyExc_ValueError, "end_time must be > 0");
        goto out;
    }
    if (max_events == 0) {
        PyErr_SetString(PyExc_ValueError, "max_events must be > 0");
        goto out;
    }

    Py_BEGIN_ALLOW_THREADS
    status = msp_run(self->sim, end_time, max_events);
    Py_END_ALLOW_THREADS
    if (status < 0) {
        handle_library_error(status);
        goto out;
    }
    ret = Py_BuildValue("i", status);
out:
    return ret;
}

static PyObject *
Simulator_finalise_tables(Simulator *self)
{
    PyObject *ret = NULL;
    int status;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    /* finalise the tables so that any uncoalesced segments are recorded */
    status = msp_finalise_tables(self->sim);
    if (status != 0) {
        handle_library_error(status);
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

static PyObject *
Simulator_compute_population_size(Simulator *self, PyObject *args)
{

    PyObject *ret = NULL;
    int sim_ret, population_id;
    double time, size;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "id", &population_id, &time)) {
        goto out;
    }
    sim_ret = msp_compute_population_size(self->sim, population_id, time, &size);
    if (sim_ret != 0) {
        handle_library_error(sim_ret);
        goto out;
    }
    ret = Py_BuildValue("d", size);
out:
    return ret;
}

static PyObject *
Simulator_fenwick_drift(Simulator *self, PyObject *args)
{

    PyObject *ret = NULL;
    int label;
    double drift = 0;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "i", &label)) {
        goto out;
    }
    if (label < 0 || label >= (int) msp_get_num_labels(self->sim)) {
        PyErr_SetString(PyExc_ValueError, "bad label ID");
        goto out;
    }
    /* TODO need a better API for this, as we should also think about the
     * drift in the GC map. */
    if (self->sim->recomb_mass_index != NULL) {
        drift = fenwick_get_numerical_drift(&self->sim->recomb_mass_index[label]);
    }
    ret = Py_BuildValue("d", drift);
out:
    return ret;
}

static FILE *
make_file(PyObject *fileobj, const char *mode)
{
    FILE *ret = NULL;
    FILE *file = NULL;
    int fileobj_fd, new_fd;

    fileobj_fd = PyObject_AsFileDescriptor(fileobj);
    if (fileobj_fd == -1) {
        goto out;
    }
    new_fd = dup(fileobj_fd);
    if (new_fd == -1) {
        PyErr_SetFromErrno(PyExc_OSError);
        goto out;
    }
    file = fdopen(new_fd, mode);
    if (file == NULL) {
        (void) close(new_fd);
        PyErr_SetFromErrno(PyExc_OSError);
        goto out;
    }
    ret = file;
out:
    return ret;
}

static PyObject *
Simulator_print_state(Simulator *self, PyObject *args)
{
    PyObject *ret = NULL;
    PyObject *fileobj;
    FILE *file = NULL;;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "O", &fileobj)) {
        goto out;
    }
    file = make_file(fileobj, "w");
    if (file == NULL) {
        goto out;
    }
    msp_print_state(self->sim, file);
    ret = Py_BuildValue("");
out:
    if (file != NULL) {
        (void) fclose(file);
    }
    return ret;
}

static PyObject *
Simulator_verify(Simulator *self, PyObject *args)
{
    PyObject *ret = NULL;
    int verify_breakpoints = 0;
    tsk_flags_t options = 0;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "|i", &verify_breakpoints)) {
        goto out;
    }
    if (verify_breakpoints) {
        options |= MSP_VERIFY_BREAKPOINTS;
    }
    msp_verify(self->sim, options);
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyMethodDef Simulator_methods[] = {
    {"run", (PyCFunction) Simulator_run, METH_VARARGS|METH_KEYWORDS,
            "Simulates until at most the specified time. Returns True "
            "if sample has coalesced and False otherwise." },
    {"reset", (PyCFunction) Simulator_reset, METH_NOARGS,
            "Resets the simulation so it's ready for another replicate."},
    {"finalise_tables", (PyCFunction) Simulator_finalise_tables, METH_NOARGS,
            "Finalises the tables so they're ready for export."},
    {"debug_demography", (PyCFunction) Simulator_debug_demography, METH_NOARGS,
            "Runs the state of the simulator forward for one demographic event."},
    {"compute_population_size",
            (PyCFunction) Simulator_compute_population_size, METH_VARARGS,
            "Computes the size of a population at a given time. Debug method."},
    {"fenwick_drift",
            (PyCFunction) Simulator_fenwick_drift, METH_VARARGS,
            "Return the numerical drift in the specified label's recombination tree. "
            "Debug method."},
    {"print_state",
            (PyCFunction) Simulator_print_state, METH_VARARGS,
            "Prints out the state of the low-level simulator. Debug method."},
    {"verify",
            (PyCFunction) Simulator_verify, METH_VARARGS,
            "Runs low-level integrity checks on the simulator's internal state."
            "This is a *debugging method only* and can result in assertions"
            "failing."},
    {NULL}  /* Sentinel */
};

static PyMemberDef Simulator_members[] = {
    {NULL}  /* Sentinel */
};


static PyGetSetDef Simulator_getsetters[] = {
    {"ancestors", (getter) Simulator_get_ancestors, NULL,
            "The ancestors" },
    {"avl_node_block_size",
            (getter) Simulator_get_avl_node_block_size, NULL,
            "The avl_node block size" },
    {"breakpoints",
            (getter) Simulator_get_breakpoints, NULL,
            "The recombination breakpoints in physical coordinates" },
    {"recombination_map",
            (getter) Simulator_get_recombination_map, NULL,
            "The recombination map" },
    {"model",
            (getter) Simulator_get_model, (setter) Simulator_set_model, NULL,
            "The simulation model." },
    {"migration_matrix",
            (getter) Simulator_get_migration_matrix, NULL,
            "The migration matrix." },
    {"node_mapping_block_size",
            (getter) Simulator_get_node_mapping_block_size, NULL,
            "The node_mapping block size" },
    {"num_ancestors",
            (getter) Simulator_get_num_ancestors, NULL,
            "The number of ancestors" },
    {"num_avl_node_blocks",
            (getter) Simulator_get_num_avl_node_blocks, NULL,
            "The number of avl_node memory blocks"},
    {"num_breakpoints",
            (getter) Simulator_get_num_breakpoints, NULL,
            "The number of recombination breakpoints" },
    {"num_common_ancestor_events",
            (getter) Simulator_get_num_common_ancestor_events, NULL,
            "The number of common_ancestor_events" },
    {"num_edges",
            (getter) Simulator_get_num_edges, NULL,
            "The number of coalescence records" },
    {"num_gene_conversion_events",
            (getter) Simulator_get_num_gene_conversion_events, NULL,
            "The number of gene_conversion_events" },
    {"num_labels",
            (getter) Simulator_get_num_labels, NULL,
            "The number of labels." },
    {"num_migration_events",
            (getter) Simulator_get_num_migration_events, NULL,
            "The number of migration events" },
    {"num_migrations",
            (getter) Simulator_get_num_migrations, NULL,
            "The number of migration records" },
    {"num_multiple_recombination_events",
            (getter) Simulator_get_num_multiple_recombination_events, NULL,
            "The number of recombination_events that occur at an "
            "existing breakpoint" },
    {"num_node_mapping_blocks",
            (getter) Simulator_get_num_node_mapping_blocks, NULL,
            "The number of node_mapping memory blocks"},
    {"num_nodes",
            (getter) Simulator_get_num_nodes, NULL,
            "The number of coalescence records" },
    {"num_populations",
            (getter) Simulator_get_num_populations, NULL,
            "The number of populations." },
    {"num_recombination_events",
            (getter) Simulator_get_num_recombination_events, NULL,
            "The number of recombination_events" },
    {"num_rejected_common_ancestor_events",
            (getter) Simulator_get_num_rejected_common_ancestor_events, NULL,
            "The number of rejected common_ancestor_events" },
    {"num_segment_blocks",
            (getter) Simulator_get_num_segment_blocks, NULL,
            "The number of segment memory blocks"},
    {"num_fenwick_rebuilds",
            (getter) Simulator_get_num_fenwick_rebuilds, NULL,
            "The number of times fenwick_rebuild was called."},
    {"population_configuration",
            (getter) Simulator_get_population_configuration, NULL,
            "The population configurations"},
    {"random_generator",
            (getter) Simulator_get_random_generator, NULL,
            "The random generator"},
    {"segment_block_size",
            (getter) Simulator_get_segment_block_size, NULL,
            "The segment block size." },
    {"sequence_length",
            (getter) Simulator_get_sequence_length, NULL,
            "The sequence length for this simulator."},
    {"store_migrations",
            (getter) Simulator_get_store_migrations, NULL,
            "True if the simulator should store migration records." },
    {"discrete_genome",
            (getter) Simulator_get_discrete_genome, NULL,
            "True if the simulator has a discrete genome." },
    {"ploidy",
            (getter) Simulator_get_ploidy, NULL,
            "Returns the simulation ploidy." },
    {"tables",
            (getter) Simulator_get_tables, NULL,
            "The tables"},
    {"time", (getter) Simulator_get_time, NULL,
            "The current simulation time" },
    {NULL}  /* Sentinel */
};

static PyTypeObject SimulatorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_msprime.Simulator",
    .tp_doc = "Simulator objects",
    .tp_basicsize = sizeof(Simulator),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor)Simulator_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_methods = Simulator_methods,
    .tp_members = Simulator_members,
    .tp_getset = Simulator_getsetters,
    .tp_init = (initproc)Simulator_init,
    .tp_new = PyType_GenericNew,
};

/*===================================================================
 * Module level functions
 *===================================================================
 */

mutation_model_t *
parse_mutation_model(PyObject *py_model)
{
    mutation_model_t *model = NULL;
    MatrixMutationModel *matrix_mutation_model = NULL;
    SLiMMutationModel *slim_mutation_model = NULL;
    InfiniteAllelesMutationModel *infinite_alleles_model = NULL;

    if (PyObject_TypeCheck(py_model, &MatrixMutationModelType)) {
        matrix_mutation_model = (MatrixMutationModel *) py_model;
        if (MatrixMutationModel_check_state(matrix_mutation_model) != 0) {
            goto out;
        }
        model = matrix_mutation_model->mutation_model;
    } else if (PyObject_TypeCheck(py_model, &SLiMMutationModelType)) {
        slim_mutation_model = (SLiMMutationModel *) py_model;
        if (SLiMMutationModel_check_state(slim_mutation_model) != 0) {
            goto out;
        }
        model = slim_mutation_model->mutation_model;
    } else if (PyObject_TypeCheck(py_model, &InfiniteAllelesMutationModelType)) {
        infinite_alleles_model = (InfiniteAllelesMutationModel *) py_model;
        if (InfiniteAllelesMutationModel_check_state(
                    infinite_alleles_model) != 0) {
            goto out;
        }
        model = infinite_alleles_model->mutation_model;
    } else {
        PyErr_SetString(PyExc_TypeError,
            "model must be an instance of MatrixMutationModel, "
            "SLiMMutationModel or InfiniteAllelesMutationModel.");
        goto out;
    }
out:
    return model;
}

static PyObject *
msprime_sim_mutations(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int flags = 0;
    int keep = 0;
    double start_time = -DBL_MAX;
    double end_time = DBL_MAX;
    LightweightTableCollection *tables = NULL;
    RandomGenerator *random_generator = NULL;
    PyObject *rate_map = NULL;
    PyObject *py_model = NULL;
    PyArrayObject *position_array = NULL;
    PyArrayObject *rate_array = NULL;
    size_t size;
    mutation_model_t *model = NULL;
    int discrete_sites = false;
    static char *kwlist[] = {
        "tables", "random_generator", "rate_map", "model",
        "discrete_sites", "keep", "start_time", "end_time", NULL};
    mutgen_t mutgen;
    int err;

    memset(&mutgen, 0, sizeof(mutgen));
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!O!O|iidd", kwlist,
            &LightweightTableCollectionType, &tables,
            &RandomGeneratorType, &random_generator,
            &PyDict_Type, &rate_map,
            &py_model, &discrete_sites, &keep, &start_time, &end_time)) {
        goto out;
    }
    if (LightweightTableCollection_check_state(tables) != 0
            || RandomGenerator_check_state(random_generator) != 0) {
        goto out;
    }
    model = parse_mutation_model(py_model);
    if (model == NULL) {
        goto out;
    }
    err = mutgen_alloc(&mutgen,
            random_generator->rng,
            tables->tables,
            model, 0);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    if (parse_rate_map(rate_map, &size, &position_array, &rate_array) != 0) {
        goto out;
    }
    err = mutgen_set_rate_map(&mutgen,
            size,
            PyArray_DATA(position_array),
            PyArray_DATA(rate_array));
    if (err != 0) {
        handle_input_error("mutation rate map", err);
        goto out;
    }
    err = mutgen_set_time_interval(&mutgen, start_time, end_time);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    if (discrete_sites) {
        flags |= MSP_DISCRETE_SITES;
    }
    if (keep) {
        flags |= MSP_KEEP_SITES;
    }
    err = mutgen_generate(&mutgen, flags);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    mutgen_free(&mutgen);
    Py_XDECREF(position_array);
    Py_XDECREF(rate_array);
    return ret;
}

static PyObject *
msprime_log_likelihood_arg(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    LightweightTableCollection *tables = NULL;
    double recombination_rate, Ne, ret_likelihood;
    static char *kwlist[] = {"tables", "Ne", "recombination_rate", NULL};
    tsk_treeseq_t ts;

    memset(&ts, 0, sizeof(ts));
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!dd", kwlist,
            &LightweightTableCollectionType, &tables, &Ne, &recombination_rate)) {
        goto out;
    }

    if (recombination_rate < 0) {
        PyErr_SetString(PyExc_ValueError, "recombination_rate must be >= 0");
        goto out;
    }

    /* Note: this will be inefficient here if we're building indexes for large
     * tables. */
    err = tsk_treeseq_init(&ts, tables->tables, TSK_BUILD_INDEXES);
    if (err != 0) {
        handle_tskit_library_error(err);
        goto out;
    }

    err = msp_log_likelihood_arg(&ts, recombination_rate, Ne, &ret_likelihood);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("d", ret_likelihood);
out:
    tsk_treeseq_free(&ts);
    return ret;
}

static PyObject *
msprime_get_gsl_version(PyObject *self)
{
    return Py_BuildValue("ii", GSL_MAJOR_VERSION, GSL_MINOR_VERSION);
}

static PyObject *
msprime_restore_gsl_error_handler(PyObject *self)
{
    gsl_set_error_handler(old_gsl_error_handler);
    return Py_BuildValue("");
}

static PyObject *
msprime_unset_gsl_error_handler(PyObject *self)
{
    /* turn off GSL error handler so we don't abort on errors. Can be restored
     * by calling restore_gsl_error_handler() */
    old_gsl_error_handler = gsl_set_error_handler_off();
    return Py_BuildValue("");
}

static PyMethodDef msprime_methods[] = {
    {"sim_mutations", (PyCFunction) msprime_sim_mutations,
            METH_VARARGS|METH_KEYWORDS,
            "Simulate mutations on a set of tables." },
    {"log_likelihood_arg", (PyCFunction) msprime_log_likelihood_arg,
            METH_VARARGS|METH_KEYWORDS,
            "Computes the log-likelihood of an ARG." },
    {"get_gsl_version", (PyCFunction) msprime_get_gsl_version, METH_NOARGS,
            "Returns the version of GSL we are linking against." },
    {"restore_gsl_error_handler", (PyCFunction) msprime_restore_gsl_error_handler,
            METH_NOARGS, "Restores the GSL error handler to its value before module import." },
    {"unset_gsl_error_handler", (PyCFunction) msprime_unset_gsl_error_handler,
            METH_NOARGS, "Unsets the GSL error handler (and stores the current value)." },
    {NULL}        /* Sentinel */
};

static struct PyModuleDef msprimemodule = {
    PyModuleDef_HEAD_INIT,
    "_msprime",
    "Low level interface for msprime",
    -1,
    msprime_methods,
    NULL, NULL, NULL, NULL
};

PyObject *
PyInit__msprime(void)
{
    PyObject *module = PyModule_Create(&msprimemodule);
    if (module == NULL) {
        return NULL;
    }
    import_array();

    /* LightweightTableCollection type */
    /* Note: we use the pre C-99 way of initialising this for the sake
     * of keeping in sync with the tskit code. This should all be
     * abstracted out into another module soon.
     */
    LightweightTableCollectionType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&LightweightTableCollectionType) < 0) {
        return NULL;
    }
    Py_INCREF(&LightweightTableCollectionType);
    PyModule_AddObject(module, "LightweightTableCollection",
            (PyObject *) &LightweightTableCollectionType);

    /* RandomGenerator type */
    if (PyType_Ready(&RandomGeneratorType) < 0) {
        return NULL;
    }
    Py_INCREF(&RandomGeneratorType);
    PyModule_AddObject(module, "RandomGenerator", (PyObject *) &RandomGeneratorType);

    /* Simulator type */
    if (PyType_Ready(&SimulatorType) < 0) {
        return NULL;
    }
    Py_INCREF(&SimulatorType);
    PyModule_AddObject(module, "Simulator", (PyObject *) &SimulatorType);

    /* BaseMutationModel type */
    if (PyType_Ready(&BaseMutationModelType) < 0) {
        return NULL;
    }
    Py_INCREF(&BaseMutationModelType);
    PyModule_AddObject(module, "BaseMutationModel",
            (PyObject *) &BaseMutationModelType);

    /* MatrixMutationModel type */
    MatrixMutationModelType.tp_base = &BaseMutationModelType;
    if (PyType_Ready(&MatrixMutationModelType) < 0) {
        return NULL;
    }
    Py_INCREF(&MatrixMutationModelType);
    PyModule_AddObject(module, "MatrixMutationModel",
            (PyObject *) &MatrixMutationModelType);

    /* SLiMMutationModel type */
    SLiMMutationModelType.tp_base = &BaseMutationModelType;
    if (PyType_Ready(&SLiMMutationModelType) < 0) {
        return NULL;
    }
    Py_INCREF(&SLiMMutationModelType);
    PyModule_AddObject(module, "SLiMMutationModel",
            (PyObject *) &SLiMMutationModelType);

    /* InfiniteAllelesMutationModel type */
    InfiniteAllelesMutationModelType.tp_base = &BaseMutationModelType;
    if (PyType_Ready(&InfiniteAllelesMutationModelType) < 0) {
        return NULL;
    }
    Py_INCREF(&InfiniteAllelesMutationModelType);
    PyModule_AddObject(module, "InfiniteAllelesMutationModel",
            (PyObject *) &InfiniteAllelesMutationModelType);

    /* Errors and constants */
    MsprimeInputError = PyErr_NewException("_msprime.InputError", NULL, NULL);
    Py_INCREF(MsprimeInputError);
    PyModule_AddObject(module, "InputError", MsprimeInputError);
    MsprimeLibraryError = PyErr_NewException("_msprime.LibraryError", NULL, NULL);
    Py_INCREF(MsprimeLibraryError);
    PyModule_AddObject(module, "LibraryError", MsprimeLibraryError);

    PyModule_AddIntConstant(module, "NODE_IS_CA_EVENT", MSP_NODE_IS_CA_EVENT);
    PyModule_AddIntConstant(module, "NODE_IS_RE_EVENT", MSP_NODE_IS_RE_EVENT);
    PyModule_AddIntConstant(module, "NODE_IS_MIG_EVENT", MSP_NODE_IS_MIG_EVENT);
    PyModule_AddIntConstant(module, "NODE_IS_CEN_EVENT", MSP_NODE_IS_CEN_EVENT);

    PyModule_AddIntConstant(module, "EXIT_COALESCENCE", MSP_EXIT_COALESCENCE);
    PyModule_AddIntConstant(module, "EXIT_MAX_EVENTS", MSP_EXIT_MAX_EVENTS);
    PyModule_AddIntConstant(module, "EXIT_MAX_TIME", MSP_EXIT_MAX_TIME);

    /* The function unset_gsl_error_handler should be called at import time,
     * ensuring we capture the value of the handler. However, just in case
     * someone calls restore_gsl_error_handler before this is called, we
     * set it to null. */
    old_gsl_error_handler = NULL;

    return module;
}
