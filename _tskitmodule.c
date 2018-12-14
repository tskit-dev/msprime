/*
** Copyright (C) 2014-2018 University of Oxford
**
** This file is part of tskit.
**
** tskit is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** tskit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with tskit.  If not, see <http://www.gnu.org/licenses/>.
*/

#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <float.h>

/* Allow C code to use kastore's dynamic API. This MUST be compiled with
 * -DKASTORE_DYNAMIC_API */
#include <kastore.h>
kas_funcptr *kas_dynamic_api;

#include "tskit.h"

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

#define MODULE_DOC \
"Low level interface for tskit"

#define SET_COLS 0
#define APPEND_COLS 1

static PyObject *TskitLibraryError;
static PyObject *TskitFileFormatError;
static PyObject *TskitVersionTooOldError;
static PyObject *TskitVersionTooNewError;


/* A lightweight wrapper for a table collection. This serves only as a wrapper
 * around a pointer and a way to data in-and-out of the low level structures
 * via the canonical dictionary encoding.
 */
typedef struct {
    PyObject_HEAD
    tsk_tbl_collection_t *tables;
} LightweightTableCollection;

/* The XTable classes each have 'lock' attribute, which is used to
 * raise an error if a Python thread attempts to access a table
 * while another Python thread is operating on it. Because tables
 * allocate memory dynamically, we cannot gaurantee safety otherwise.
 * The locks are set before the GIL is released and unset afterwards.
 * Because C code executed here represents atomic Python operations
 * (while the GIL is held), this should be safe */

typedef struct _TableCollection {
    PyObject_HEAD
    tsk_tbl_collection_t *tables;
} TableCollection;

 /* The table pointer in each of the Table classes either points to locally
  * allocated memory or to the table stored in a tbl_collection_t. If we're
  * using the memory in a tbl_collection_t, we keep a reference to the
  * TableCollection object to ensure that the memory isn't free'd while a
  * reference to the table itself is live. */
typedef struct {
    PyObject_HEAD
    bool locked;
    tsk_individual_tbl_t *table;
    TableCollection *tables;
} IndividualTable;

typedef struct {
    PyObject_HEAD
    bool locked;
    tsk_node_tbl_t *table;
    TableCollection *tables;
} NodeTable;

typedef struct {
    PyObject_HEAD
    bool locked;
    tsk_edge_tbl_t *table;
    TableCollection *tables;
} EdgeTable;

typedef struct {
    PyObject_HEAD
    bool locked;
    tsk_site_tbl_t *table;
    TableCollection *tables;
} SiteTable;

typedef struct {
    PyObject_HEAD
    bool locked;
    tsk_mutation_tbl_t *table;
    TableCollection *tables;
} MutationTable;

typedef struct {
    PyObject_HEAD
    bool locked;
    tsk_migration_tbl_t *table;
    TableCollection *tables;
} MigrationTable;

typedef struct {
    PyObject_HEAD
    bool locked;
    tsk_population_tbl_t *table;
    TableCollection *tables;
} PopulationTable;

typedef struct {
    PyObject_HEAD
    bool locked;
    tsk_provenance_tbl_t *table;
    TableCollection *tables;
} ProvenanceTable;

typedef struct {
    PyObject_HEAD
    tsk_treeseq_t *tree_sequence;
} TreeSequence;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    tsk_tree_t *sparse_tree;
} SparseTree;

typedef struct {
    PyObject_HEAD
    SparseTree *sparse_tree;
    int first;
} SparseTreeIterator;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    tsk_diff_iter_t *tree_diff_iterator;
} TreeDiffIterator;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    tsk_vcf_converter_t *tsk_vcf_converter;
} VcfConverter;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    tsk_hapgen_t *haplotype_generator;
} HaplotypeGenerator;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    tsk_vargen_t *variant_generator;
} VariantGenerator;

typedef struct {
    PyObject_HEAD
    TreeSequence *tree_sequence;
    tsk_ld_calc_t *ld_calc;
} LdCalculator;

static void
handle_library_error(int err)
{
    if (tsk_is_kas_error(err)) {
        PyErr_SetString(TskitFileFormatError, tsk_strerror(err));
    } else {
        switch (err) {
            case TSK_ERR_FILE_VERSION_TOO_NEW:
                PyErr_SetString(TskitVersionTooNewError, tsk_strerror(err));
                break;
            case TSK_ERR_FILE_VERSION_TOO_OLD:
                PyErr_SetString(TskitVersionTooOldError, tsk_strerror(err));
                break;
            case TSK_ERR_FILE_FORMAT:
                PyErr_SetString(TskitFileFormatError, tsk_strerror(err));
                break;
            case TSK_ERR_OUT_OF_BOUNDS:
                PyErr_SetString(PyExc_IndexError, tsk_strerror(err));
                break;
            default:
                PyErr_SetString(TskitLibraryError, tsk_strerror(err));
        }
    }
}

static int
parse_sample_ids(PyObject *py_samples, tsk_treeseq_t *ts, size_t *num_samples,
        tsk_id_t **samples)
{
    int ret = -1;
    PyObject *item;
    Py_ssize_t j, num_samples_local;
    tsk_id_t *samples_local = NULL;

    num_samples_local = PyList_Size(py_samples);
    if (num_samples_local < 2) {
        PyErr_SetString(PyExc_ValueError, "Must provide at least 2 samples");
        goto out;
    }
    samples_local = PyMem_Malloc(num_samples_local * sizeof(tsk_id_t));
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
        samples_local[j] = (tsk_id_t) PyLong_AsLong(item);
        if (samples_local[j] < 0
                || samples_local[j] >= (tsk_id_t) tsk_treeseq_get_num_nodes(ts)) {
            PyErr_SetString(PyExc_ValueError, "node ID out of bounds");
            goto out;
        }
        if (! tsk_treeseq_is_sample(ts, samples_local[j])) {
            PyErr_SetString(PyExc_ValueError, "Specified node is not a sample");
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


static PyObject *
convert_node_id_list(tsk_id_t *children, size_t num_children)
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
make_metadata(const char *metadata, Py_ssize_t length)
{
    const char *m = metadata == NULL? "": metadata;
    return PyBytes_FromStringAndSize(m, length);
}

static PyObject *
make_mutation(tsk_mutation_t *mutation)
{
    PyObject *ret = NULL;
    PyObject* metadata = NULL;

    metadata = make_metadata(mutation->metadata, (Py_ssize_t) mutation->metadata_length);
    if (metadata == NULL) {
        goto out;
    }
    ret = Py_BuildValue("iis#iO", mutation->site, mutation->node, mutation->derived_state,
            (Py_ssize_t) mutation->derived_state_length, mutation->parent,
            metadata);
out:
    Py_XDECREF(metadata);
    return ret;
}

static PyObject *
make_mutation_id_list(tsk_mutation_t *mutations, size_t length)
{
    PyObject *ret = NULL;
    PyObject *t;
    PyObject *item;
    size_t j;

    t = PyTuple_New(length);
    if (t == NULL) {
        goto out;
    }
    for (j = 0; j < length; j++) {
        item = Py_BuildValue("i", mutations[j].id);
        if (item == NULL) {
            Py_DECREF(t);
            goto out;
        }
        PyTuple_SET_ITEM(t, j, item);
    }
    ret = t;
out:
    return ret;
}

static PyObject *
make_population(tsk_population_t *population)
{
    PyObject *ret = NULL;
    PyObject *metadata = make_metadata(population->metadata,
            (Py_ssize_t) population->metadata_length);

    ret = Py_BuildValue("(O)", metadata);
    return ret;
}

static PyObject *
make_provenance(tsk_provenance_t *provenance)
{
    PyObject *ret = NULL;

    ret = Py_BuildValue("s#s#",
            provenance->timestamp, (Py_ssize_t) provenance->timestamp_length,
            provenance->record, (Py_ssize_t) provenance->record_length);
    return ret;
}

static PyObject *
make_individual_row(tsk_individual_t *r)
{
    PyObject *ret = NULL;
    PyObject *metadata = make_metadata(r->metadata, (Py_ssize_t) r->metadata_length);
    PyArrayObject *location = NULL;
    npy_intp dims;

    dims = (npy_intp) r->location_length;
    location = (PyArrayObject *) PyArray_SimpleNew(1, &dims, NPY_FLOAT64);
    if (metadata == NULL || location == NULL) {
        goto out;
    }
    memcpy(PyArray_DATA(location), r->location, r->location_length * sizeof(double));
    ret = Py_BuildValue("IOO", (unsigned int) r->flags, location, metadata);
out:
    Py_XDECREF(location);
    Py_XDECREF(metadata);
    return ret;
}

static PyObject *
make_individual_object(tsk_individual_t *r)
{
    PyObject *ret = NULL;
    PyObject *metadata = make_metadata(r->metadata, (Py_ssize_t) r->metadata_length);
    PyArrayObject *location = NULL;
    PyArrayObject *nodes = NULL;
    npy_intp dims;

    dims = (npy_intp) r->location_length;
    location = (PyArrayObject *) PyArray_SimpleNew(1, &dims, NPY_FLOAT64);
    dims = (npy_intp) r->nodes_length;
    nodes = (PyArrayObject *) PyArray_SimpleNew(1, &dims, NPY_INT32);
    if (metadata == NULL || location == NULL || nodes == NULL) {
        goto out;
    }
    memcpy(PyArray_DATA(location), r->location, r->location_length * sizeof(double));
    memcpy(PyArray_DATA(nodes), r->nodes, r->nodes_length * sizeof(tsk_id_t));
    ret = Py_BuildValue("IOOO", (unsigned int) r->flags, location, metadata, nodes);
out:
    Py_XDECREF(location);
    Py_XDECREF(metadata);
    Py_XDECREF(nodes);
    return ret;
}

static PyObject *
make_node(tsk_node_t *r)
{
    PyObject *ret = NULL;
    PyObject* metadata = make_metadata(r->metadata, (Py_ssize_t) r->metadata_length);
    if (metadata == NULL) {
        goto out;
    }
    ret = Py_BuildValue("IdiiO",
        (unsigned int) r->flags, r->time, (int) r->population, (int) r->individual, metadata);
out:
    Py_XDECREF(metadata);
    return ret;
}

static PyObject *
make_edge(tsk_edge_t *edge)
{
    return Py_BuildValue("ddii",
            edge->left, edge->right, (int) edge->parent, (int) edge->child);
}

static PyObject *
make_migration(tsk_migration_t *r)
{
    int source = r->source == TSK_NULL_POPULATION ? -1: r->source;
    int dest = r->dest == TSK_NULL_POPULATION ? -1: r->dest;
    PyObject *ret = NULL;

    ret = Py_BuildValue("ddiiid",
            r->left, r->right, (int) r->node, source, dest, r->time);
    return ret;
}

static PyObject *
make_site_row(tsk_site_t *site)
{
    PyObject *ret = NULL;
    PyObject* metadata = NULL;

    metadata = make_metadata(site->metadata, (Py_ssize_t) site->metadata_length);
    if (metadata == NULL) {
        goto out;
    }
    ret = Py_BuildValue("ds#O", site->position, site->ancestral_state,
            (Py_ssize_t) site->ancestral_state_length, metadata);
out:
    Py_XDECREF(metadata);
    return ret;
}

static PyObject *
make_site_object(tsk_site_t *site)
{
    PyObject *ret = NULL;
    PyObject *mutations = NULL;
    PyObject* metadata = NULL;

    metadata = make_metadata(site->metadata, (Py_ssize_t) site->metadata_length);
    if (metadata == NULL) {
        goto out;
    }
    mutations = make_mutation_id_list(site->mutations, site->mutations_length);
    if (mutations == NULL) {
        goto out;
    }
    /* TODO should reorder this tuple, as it's not very logical. */
    ret = Py_BuildValue("ds#OnO", site->position, site->ancestral_state,
            (Py_ssize_t) site->ancestral_state_length, mutations,
            (Py_ssize_t) site->id, metadata);
out:
    Py_XDECREF(mutations);
    Py_XDECREF(metadata);
    return ret;
}

static PyObject *
make_alleles(tsk_variant_t *variant)
{
    PyObject *ret = NULL;
    PyObject *item, *t;
    size_t j;

    t = PyTuple_New(variant->num_alleles);
    if (t == NULL) {
        goto out;
    }
    for (j = 0; j < variant->num_alleles; j++) {
        item = Py_BuildValue("s#", variant->alleles[j], variant->allele_lengths[j]);
        if (item == NULL) {
            Py_DECREF(t);
            goto out;
        }
        PyTuple_SET_ITEM(t, j, item);
    }
    ret = t;
out:
    return ret;
}

static PyObject *
make_variant(tsk_variant_t *variant, size_t num_samples)
{
    PyObject *ret = NULL;
    npy_intp dims = num_samples;
    PyObject *alleles = make_alleles(variant);
    PyArrayObject *genotypes = (PyArrayObject *) PyArray_SimpleNew(1, &dims, NPY_UINT8);

    /* TODO update this to account for 16 bit variants when we provide the
     * high-level interface. */
    if (genotypes == NULL || alleles == NULL) {
        goto out;
    }
    memcpy(PyArray_DATA(genotypes), variant->genotypes.u8, num_samples * sizeof(uint8_t));
    ret = Py_BuildValue("iOO", variant->site->id, genotypes, alleles);
out:
    Py_XDECREF(genotypes);
    Py_XDECREF(alleles);
    return ret;
}

static PyObject *
convert_sites(tsk_site_t *sites, size_t num_sites)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_site = NULL;
    size_t j;

    l = PyList_New(num_sites);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < num_sites; j++) {
        py_site = make_site_object(&sites[j]);
        if (py_site == NULL) {
            Py_DECREF(l);
            goto out;
        }
        PyList_SET_ITEM(l, j, py_site);
    }
    ret = l;
out:
    return ret;
}

/*===================================================================
 * General table code.
 *===================================================================
 */

/*
 * Retrieves the PyObject* corresponding the specified key in the
 * specified dictionary. If required is true, raise a TypeError if the
 * value is None.
 *
 * NB This returns a *borrowed reference*, so don't DECREF it!
 */
static PyObject *
get_table_dict_value(PyObject *dict, const char *key_str, bool required)
{
    PyObject *ret = NULL;

    ret = PyDict_GetItemString(dict, key_str);
    if (ret == NULL) {
        PyErr_Format(PyExc_ValueError, "'%s' not specified", key_str);
    }
    if (required && ret == Py_None) {
        PyErr_Format(PyExc_TypeError, "'%s' is required", key_str);
        ret = NULL;
    }
    return ret;
}

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
    if (shape[0] != *num_rows + 1) {
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

static int
parse_individual_table_dict(tsk_individual_tbl_t *table, PyObject *dict, bool clear_table)
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

    if (clear_table) {
        err = tsk_individual_tbl_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_individual_tbl_append_columns(table, num_rows,
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
parse_node_table_dict(tsk_node_tbl_t *table, PyObject *dict, bool clear_table)
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
    if (clear_table) {
        err = tsk_node_tbl_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_node_tbl_append_columns(table, num_rows,
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
parse_edge_table_dict(tsk_edge_tbl_t *table, PyObject *dict, bool clear_table)
{
    int ret = -1;
    int err;
    size_t num_rows = 0;
    PyObject *left_input = NULL;
    PyArrayObject *left_array = NULL;
    PyObject *right_input = NULL;
    PyArrayObject *right_array = NULL;
    PyObject *parent_input = NULL;
    PyArrayObject *parent_array = NULL;
    PyObject *child_input = NULL;
    PyArrayObject *child_array = NULL;

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

    if (clear_table) {
        err = tsk_edge_tbl_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_edge_tbl_append_columns(table, num_rows,
            PyArray_DATA(left_array), PyArray_DATA(right_array),
            PyArray_DATA(parent_array), PyArray_DATA(child_array));
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
    return ret;
}

static int
parse_migration_table_dict(tsk_migration_tbl_t *table, PyObject *dict, bool clear_table)
{
    int err;
    int ret = -1;
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

    if (clear_table) {
        err = tsk_migration_tbl_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_migration_tbl_append_columns(table, num_rows,
        PyArray_DATA(left_array), PyArray_DATA(right_array), PyArray_DATA(node_array),
        PyArray_DATA(source_array), PyArray_DATA(dest_array), PyArray_DATA(time_array));
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
    return ret;
}

static int
parse_site_table_dict(tsk_site_tbl_t *table, PyObject *dict, bool clear_table)
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

    if (clear_table) {
        err = tsk_site_tbl_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_site_tbl_append_columns(table, num_rows,
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
parse_mutation_table_dict(tsk_mutation_tbl_t *table, PyObject *dict, bool clear_table)
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
    PyObject *parent_input = NULL;
    PyArrayObject *parent_array = NULL;
    tsk_id_t *parent_data;
    PyObject *metadata_input = NULL;
    PyArrayObject *metadata_array = NULL;
    PyObject *metadata_offset_input = NULL;
    PyArrayObject *metadata_offset_array = NULL;
    char *metadata_data;
    uint32_t *metadata_offset_data;

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

    if (clear_table) {
        err = tsk_mutation_tbl_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_mutation_tbl_append_columns(table, num_rows,
            PyArray_DATA(site_array), PyArray_DATA(node_array),
            parent_data, PyArray_DATA(derived_state_array),
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
    return ret;
}

static int
parse_population_table_dict(tsk_population_tbl_t *table, PyObject *dict, bool clear_table)
{
    int err;
    int ret = -1;
    size_t num_rows, metadata_length;
    PyObject *metadata_input = NULL;
    PyArrayObject *metadata_array = NULL;
    PyObject *metadata_offset_input = NULL;
    PyArrayObject *metadata_offset_array = NULL;

    /* Get the inputs */
    metadata_input = get_table_dict_value(dict, "metadata", true);
    if (metadata_input == NULL) {
        goto out;
    }
    metadata_offset_input = get_table_dict_value(dict, "metadata_offset", true);
    if (metadata_offset_input == NULL) {
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

    if (clear_table) {
        err = tsk_population_tbl_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_population_tbl_append_columns(table, num_rows,
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
parse_provenance_table_dict(tsk_provenance_tbl_t *table, PyObject *dict, bool clear_table)
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
        err = tsk_provenance_tbl_clear(table);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    err = tsk_provenance_tbl_append_columns(table, num_rows,
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
parse_table_collection_dict(tsk_tbl_collection_t *tables, PyObject *tables_dict)
{
    int ret = -1;
    PyObject *value = NULL;

    value = get_table_dict_value(tables_dict, "sequence_length", true);
    if (value == NULL) {
        goto out;
    }
    if (!PyNumber_Check(value)) {
        PyErr_Format(PyExc_TypeError, "'sequence_length' is not number");
        goto out;
    }
    tables->sequence_length = PyFloat_AsDouble(value);

    /* individuals */
    value = get_table_dict_value(tables_dict, "individuals", true);
    if (value == NULL) {
        goto out;
    }
    if (!PyDict_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "not a dictionary");
        goto out;
    }
    if (parse_individual_table_dict(tables->individuals, value, true) != 0) {
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
    if (parse_node_table_dict(tables->nodes, value, true) != 0) {
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
    if (parse_edge_table_dict(tables->edges, value, true) != 0) {
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
    if (parse_migration_table_dict(tables->migrations, value, true) != 0) {
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
    if (parse_site_table_dict(tables->sites, value, true) != 0) {
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
    if (parse_mutation_table_dict(tables->mutations, value, true) != 0) {
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
    if (parse_population_table_dict(tables->populations, value, true) != 0) {
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
    if (parse_provenance_table_dict(tables->provenances, value, true) != 0) {
        goto out;
    }

    ret = 0;
out:
    return ret;
}

static int
write_table_arrays(tsk_tbl_collection_t *tables, PyObject *dict)
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
    };
    int ret = -1;
    PyObject *array = NULL;
    PyObject *table_dict = NULL;
    size_t j;
    struct table_col *col;

    struct table_col individual_cols[] = {
        {"flags",
            (void *) tables->individuals->flags, tables->individuals->num_rows, NPY_UINT32},
        {"location",
            (void *) tables->individuals->location, tables->individuals->location_length,
            NPY_FLOAT64},
        {"location_offset",
            (void *) tables->individuals->location_offset, tables->individuals->num_rows + 1,
            NPY_UINT32},
        {"metadata",
            (void *) tables->individuals->metadata, tables->individuals->metadata_length,
            NPY_INT8},
        {"metadata_offset",
            (void *) tables->individuals->metadata_offset, tables->individuals->num_rows + 1,
            NPY_UINT32},
        {NULL},
    };

    struct table_col node_cols[] = {
        {"time",
            (void *) tables->nodes->time, tables->nodes->num_rows, NPY_FLOAT64},
        {"flags",
            (void *) tables->nodes->flags, tables->nodes->num_rows, NPY_UINT32},
        {"population",
            (void *) tables->nodes->population, tables->nodes->num_rows, NPY_INT32},
        {"individual",
            (void *) tables->nodes->individual, tables->nodes->num_rows, NPY_INT32},
        {"metadata",
            (void *) tables->nodes->metadata, tables->nodes->metadata_length, NPY_INT8},
        {"metadata_offset",
            (void *) tables->nodes->metadata_offset, tables->nodes->num_rows + 1, NPY_UINT32},
        {NULL},
    };

    struct table_col edge_cols[] = {
        {"left", (void *) tables->edges->left, tables->edges->num_rows, NPY_FLOAT64},
        {"right", (void *) tables->edges->right, tables->edges->num_rows, NPY_FLOAT64},
        {"parent", (void *) tables->edges->parent, tables->edges->num_rows, NPY_INT32},
        {"child", (void *) tables->edges->child, tables->edges->num_rows, NPY_INT32},
        {NULL},
    };

    struct table_col migration_cols[] = {
        {"left",
            (void *) tables->migrations->left, tables->migrations->num_rows,  NPY_FLOAT64},
        {"right",
            (void *) tables->migrations->right, tables->migrations->num_rows,  NPY_FLOAT64},
        {"node",
            (void *) tables->migrations->node, tables->migrations->num_rows,  NPY_INT32},
        {"source",
            (void *) tables->migrations->source, tables->migrations->num_rows,  NPY_INT32},
        {"dest",
            (void *) tables->migrations->dest, tables->migrations->num_rows,  NPY_INT32},
        {"time",
            (void *) tables->migrations->time, tables->migrations->num_rows,  NPY_FLOAT64},
        {NULL},
    };

    struct table_col site_cols[] = {
        {"position",
            (void *) tables->sites->position, tables->sites->num_rows, NPY_FLOAT64},
        {"ancestral_state",
            (void *) tables->sites->ancestral_state, tables->sites->ancestral_state_length,
            NPY_INT8},
        {"ancestral_state_offset",
            (void *) tables->sites->ancestral_state_offset, tables->sites->num_rows + 1,
            NPY_UINT32},
        {"metadata",
            (void *) tables->sites->metadata, tables->sites->metadata_length, NPY_INT8},
        {"metadata_offset",
            (void *) tables->sites->metadata_offset, tables->sites->num_rows + 1, NPY_UINT32},
        {NULL},
    };

    struct table_col mutation_cols[] = {
        {"site",
            (void *) tables->mutations->site, tables->mutations->num_rows, NPY_INT32},
        {"node",
            (void *) tables->mutations->node, tables->mutations->num_rows, NPY_INT32},
        {"parent",
            (void *) tables->mutations->parent, tables->mutations->num_rows, NPY_INT32},
        {"derived_state",
            (void *) tables->mutations->derived_state,
            tables->mutations->derived_state_length, NPY_INT8},
        {"derived_state_offset",
            (void *) tables->mutations->derived_state_offset,
            tables->mutations->num_rows + 1, NPY_UINT32},
        {"metadata",
            (void *) tables->mutations->metadata,
            tables->mutations->metadata_length, NPY_INT8},
        {"metadata_offset",
            (void *) tables->mutations->metadata_offset,
            tables->mutations->num_rows + 1, NPY_UINT32},
        {NULL},
    };

    struct table_col population_cols[] = {
        {"metadata", (void *) tables->populations->metadata,
            tables->populations->metadata_length, NPY_INT8},
        {"metadata_offset", (void *) tables->populations->metadata_offset,
            tables->populations->num_rows+ 1, NPY_UINT32},
        {NULL},
    };

    struct table_col provenance_cols[] = {
        {"timestamp", (void *) tables->provenances->timestamp,
            tables->provenances->timestamp_length, NPY_INT8},
        {"timestamp_offset", (void *) tables->provenances->timestamp_offset,
            tables->provenances->num_rows+ 1, NPY_UINT32},
        {"record", (void *) tables->provenances->record,
            tables->provenances->record_length, NPY_INT8},
        {"record_offset", (void *) tables->provenances->record_offset,
            tables->provenances->num_rows + 1, NPY_UINT32},
        {NULL},
    };

    struct table_desc table_descs[] = {
        {"individuals", individual_cols},
        {"nodes", node_cols},
        {"edges", edge_cols},
        {"migrations", migration_cols},
        {"sites", site_cols},
        {"mutations", mutation_cols},
        {"populations", population_cols},
        {"provenances", provenance_cols},
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
dump_tables_dict(tsk_tbl_collection_t *tables)
{
    PyObject *ret = NULL;
    PyObject *dict = NULL;
    PyObject *val = NULL;
    int err;

    dict = PyDict_New();
    if (dict == NULL) {
        goto out;
    }
    val = Py_BuildValue("d", tables->sequence_length);
    if (val == NULL) {
        goto out;
    }
    if (PyDict_SetItemString(dict, "sequence_length", val) != 0) {
        goto out;
    }
    Py_DECREF(val);
    val = NULL;

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
        tsk_tbl_collection_free(self->tables);
        PyMem_Free(self->tables);
        self->tables = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
LightweightTableCollection_init(LightweightTableCollection *self, PyObject *args, PyObject *kwds)
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
    err = tsk_tbl_collection_alloc(self->tables, TSK_ALLOC_TABLES);
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

/*===================================================================
 * IndividualTable
 *===================================================================
 */

static int
IndividualTable_check_state(IndividualTable *self)
{
    int ret = -1;
    if (self->table == NULL) {
        PyErr_SetString(PyExc_SystemError, "IndividualTable not initialised");
        goto out;
    }
    if (self->locked) {
        PyErr_SetString(PyExc_RuntimeError, "IndividualTable in use by other thread.");
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static void
IndividualTable_dealloc(IndividualTable* self)
{
    if (self->tables != NULL) {
        Py_DECREF(self->tables);
    } else if (self->table != NULL) {
        tsk_individual_tbl_free(self->table);
        PyMem_Free(self->table);
        self->table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
IndividualTable_init(IndividualTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"max_rows_increment", NULL};
    Py_ssize_t max_rows_increment = 0;
    Py_ssize_t max_position_length_increment = 0;
    Py_ssize_t max_metadata_length_increment = 0;

    self->table = NULL;
    self->locked = false;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist,
                &max_rows_increment)) {
        goto out;
    }
    if (max_rows_increment < 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    self->table = PyMem_Malloc(sizeof(tsk_individual_tbl_t));
    if (self->table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tsk_individual_tbl_alloc(self->table,
            (size_t) max_rows_increment,
            (size_t) max_position_length_increment,
            (size_t) max_metadata_length_increment);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
IndividualTable_add_row(IndividualTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    unsigned int flags = 0;
    PyObject *py_metadata = Py_None;
    PyObject *py_location = Py_None;
    PyArrayObject *location_array = NULL;
    double *location_data = NULL;
    tsk_tbl_size_t location_length = 0;
    char *metadata = "";
    Py_ssize_t metadata_length = 0;
    npy_intp *shape;
    static char *kwlist[] = {"flags", "location", "metadata", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iOO", kwlist,
                &flags, &py_location, &py_metadata)) {
        goto out;
    }
    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    if (py_metadata != Py_None) {
        if (PyBytes_AsStringAndSize(py_metadata, &metadata, &metadata_length) < 0) {
            goto out;
        }
    }
    if (py_location != Py_None) {
        /* This ensures that only 1D arrays are accepted. */
        location_array = (PyArrayObject *) PyArray_FromAny(py_location,
                PyArray_DescrFromType(NPY_FLOAT64), 1, 1,
                NPY_ARRAY_IN_ARRAY, NULL);
        if (location_array == NULL) {
            goto out;
        }
        shape = PyArray_DIMS(location_array);
        location_length = (tsk_tbl_size_t) shape[0];
        location_data = PyArray_DATA(location_array);
    }
    err = tsk_individual_tbl_add_row(self->table, (uint32_t) flags,
            location_data, location_length, metadata, metadata_length);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("i", err);
out:
    Py_XDECREF(location_array);
    return ret;
}


/* Forward declaration */
static PyTypeObject IndividualTableType;

static PyObject *
IndividualTable_equals(IndividualTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    IndividualTable *other = NULL;

    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "O!", &IndividualTableType, &other)) {
        goto out;
    }
    ret = Py_BuildValue("i", tsk_individual_tbl_equals(self->table, other->table));
out:
    return ret;
}

static PyObject *
IndividualTable_get_row(IndividualTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    int err;
    Py_ssize_t row_id;
    tsk_individual_t individual;

    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &row_id)) {
        goto out;
    }
    err = tsk_individual_tbl_get_row(self->table, (size_t) row_id, &individual);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_individual_row(&individual);
out:
    return ret;
}

static PyObject *
IndividualTable_parse_dict_arg(IndividualTable *self, PyObject *args, bool clear_table)
{
    int err;
    PyObject *ret = NULL;
    PyObject *dict = NULL;

    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) {
        goto out;
    }
    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    err = parse_individual_table_dict(self->table, dict, clear_table);
    if (err != 0) {
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
IndividualTable_append_columns(IndividualTable *self, PyObject *args)
{
    return IndividualTable_parse_dict_arg(self, args, false);
}

static PyObject *
IndividualTable_set_columns(IndividualTable *self, PyObject *args)
{
    return IndividualTable_parse_dict_arg(self, args, true);
}

static PyObject *
IndividualTable_clear(IndividualTable *self)
{
    PyObject *ret = NULL;
    int err;

    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    err = tsk_individual_tbl_clear(self->table);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
IndividualTable_truncate(IndividualTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t num_rows;
    int err;

    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &num_rows)) {
        goto out;
    }
    if (num_rows < 0 || num_rows > (Py_ssize_t) self->table->num_rows) {
        PyErr_SetString(PyExc_ValueError, "num_rows out of bounds");
        goto out;
    }
    err = tsk_individual_tbl_truncate(self->table, (size_t) num_rows);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
IndividualTable_get_max_rows_increment(IndividualTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows_increment);
out:
    return ret;
}

static PyObject *
IndividualTable_get_num_rows(IndividualTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->num_rows);
out:
    return ret;
}

static PyObject *
IndividualTable_get_max_rows(IndividualTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows);
out:
    return ret;
}

static PyObject *
IndividualTable_get_flags(IndividualTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->num_rows, self->table->flags,
            NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyObject *
IndividualTable_get_location(IndividualTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->location_length,
            self->table->location, NPY_FLOAT64, sizeof(double));
out:
    return ret;
}

static PyObject *
IndividualTable_get_location_offset(IndividualTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->num_rows + 1,
            self->table->location_offset, NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyObject *
IndividualTable_get_metadata(IndividualTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->metadata_length,
            self->table->metadata, NPY_INT8, sizeof(char));
out:
    return ret;
}

static PyObject *
IndividualTable_get_metadata_offset(IndividualTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (IndividualTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->num_rows + 1,
            self->table->metadata_offset, NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyGetSetDef IndividualTable_getsetters[] = {
    {"max_rows_increment",
        (getter) IndividualTable_get_max_rows_increment, NULL, "The size increment"},
    {"num_rows", (getter) IndividualTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"max_rows", (getter) IndividualTable_get_max_rows, NULL,
        "The current maximum number of rows in the table."},
    {"flags", (getter) IndividualTable_get_flags, NULL, "The flags array"},
    {"location", (getter) IndividualTable_get_location, NULL, "The location array"},
    {"location_offset", (getter) IndividualTable_get_location_offset, NULL,
        "The location offset array"},
    {"metadata", (getter) IndividualTable_get_metadata, NULL, "The metadata array"},
    {"metadata_offset", (getter) IndividualTable_get_metadata_offset, NULL,
        "The metadata offset array"},
    {NULL}  /* Sentinel */
};

static PyMethodDef IndividualTable_methods[] = {
    {"add_row", (PyCFunction) IndividualTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    {"get_row", (PyCFunction) IndividualTable_get_row, METH_VARARGS,
        "Returns the kth row in this table."},
    {"equals", (PyCFunction) IndividualTable_equals, METH_VARARGS,
        "Returns true if the specified individual table is equal."},
    {"append_columns", (PyCFunction) IndividualTable_append_columns,
        METH_VARARGS|METH_KEYWORDS,
        "Appends the data in the specified arrays into the columns."},
    {"set_columns", (PyCFunction) IndividualTable_set_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the specified arrays into the columns."},
    {"clear", (PyCFunction) IndividualTable_clear, METH_NOARGS,
        "Clears this table."},
    {"truncate", (PyCFunction) IndividualTable_truncate, METH_VARARGS,
        "Truncates this table to the specified number of rows."},
    {NULL}  /* Sentinel */
};

static PyTypeObject IndividualTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_tskit.IndividualTable",             /* tp_name */
    sizeof(IndividualTable),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)IndividualTable_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "IndividualTable objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    IndividualTable_methods,             /* tp_methods */
    0,                             /* tp_members */
    IndividualTable_getsetters,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)IndividualTable_init,      /* tp_init */
};


/*===================================================================
 * NodeTable
 *===================================================================
 */

static int
NodeTable_check_state(NodeTable *self)
{
    int ret = -1;
    if (self->table == NULL) {
        PyErr_SetString(PyExc_SystemError, "NodeTable not initialised");
        goto out;
    }
    if (self->locked) {
        PyErr_SetString(PyExc_RuntimeError, "NodeTable in use by other thread.");
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static void
NodeTable_dealloc(NodeTable* self)
{
    if (self->tables != NULL) {
        Py_DECREF(self->tables);
    } else if (self->table != NULL) {
        tsk_node_tbl_free(self->table);
        PyMem_Free(self->table);
        self->table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
NodeTable_init(NodeTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"max_rows_increment", NULL};
    Py_ssize_t max_rows_increment = 0;
    Py_ssize_t max_metadata_length_increment = 0;

    self->table = NULL;
    self->locked = false;
    self->tables = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist,
                &max_rows_increment)) {
        goto out;
    }
    if (max_rows_increment < 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    self->table = PyMem_Malloc(sizeof(tsk_node_tbl_t));
    if (self->table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tsk_node_tbl_alloc(self->table, (size_t) max_rows_increment,
            (size_t) max_metadata_length_increment);
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
    int individual = -1;
    PyObject *py_metadata = Py_None;
    char *metadata = "";
    Py_ssize_t metadata_length = 0;
    static char *kwlist[] = {"flags", "time", "population", "individual", "metadata", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|idiiO", kwlist,
                &flags, &time, &population, &individual, &py_metadata)) {
        goto out;
    }
    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    if (py_metadata != Py_None) {
        if (PyBytes_AsStringAndSize(py_metadata, &metadata, &metadata_length) < 0) {
            goto out;
        }
    }
    err = tsk_node_tbl_add_row(self->table, (uint32_t) flags, time,
            (tsk_id_t) population, individual, metadata, metadata_length);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("i", err);
out:
    return ret;
}

/* Forward declaration */
static PyTypeObject NodeTableType;

static PyObject *
NodeTable_equals(NodeTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    NodeTable *other = NULL;

    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "O!", &NodeTableType, &other)) {
        goto out;
    }
    ret = Py_BuildValue("i", tsk_node_tbl_equals(self->table, other->table));
out:
    return ret;
}

static PyObject *
NodeTable_get_row(NodeTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    int err;
    Py_ssize_t row_id;
    tsk_node_t node;

    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &row_id)) {
        goto out;
    }
    err = tsk_node_tbl_get_row(self->table, (size_t) row_id, &node);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_node(&node);
out:
    return ret;
}

static PyObject *
NodeTable_parse_dict_arg(NodeTable *self, PyObject *args, bool clear_table)
{
    int err;
    PyObject *ret = NULL;
    PyObject *dict = NULL;

    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) {
        goto out;
    }
    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    err = parse_node_table_dict(self->table, dict, clear_table);
    if (err != 0) {
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
NodeTable_append_columns(NodeTable *self, PyObject *args)
{
    return NodeTable_parse_dict_arg(self, args, false);
}

static PyObject *
NodeTable_set_columns(NodeTable *self, PyObject *args)
{
    return NodeTable_parse_dict_arg(self, args, true);
}

static PyObject *
NodeTable_clear(NodeTable *self)
{
    PyObject *ret = NULL;
    int err;

    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    err = tsk_node_tbl_clear(self->table);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
NodeTable_truncate(NodeTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t num_rows;
    int err;

    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &num_rows)) {
        goto out;
    }
    if (num_rows < 0 || num_rows > (Py_ssize_t) self->table->num_rows) {
        PyErr_SetString(PyExc_ValueError, "num_rows out of bounds");
        goto out;
    }
    err = tsk_node_tbl_truncate(self->table, (size_t) num_rows);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
NodeTable_get_max_rows_increment(NodeTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows_increment);
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
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->num_rows);
out:
    return ret;
}

static PyObject *
NodeTable_get_max_rows(NodeTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows);
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
    ret = table_get_column_array(self->table->num_rows, self->table->time,
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
    ret = table_get_column_array(self->table->num_rows, self->table->flags,
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
    ret = table_get_column_array(self->table->num_rows, self->table->population,
            NPY_INT32, sizeof(int32_t));
out:
    return ret;
}

static PyObject *
NodeTable_get_individual(NodeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->num_rows, self->table->individual,
            NPY_INT32, sizeof(int32_t));
out:
    return ret;
}

static PyObject *
NodeTable_get_metadata(NodeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->metadata_length,
            self->table->metadata, NPY_INT8, sizeof(char));
out:
    return ret;
}

static PyObject *
NodeTable_get_metadata_offset(NodeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (NodeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->num_rows + 1,
            self->table->metadata_offset, NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyGetSetDef NodeTable_getsetters[] = {
    {"max_rows_increment",
        (getter) NodeTable_get_max_rows_increment, NULL, "The size increment"},
    {"num_rows", (getter) NodeTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"max_rows", (getter) NodeTable_get_max_rows, NULL,
        "The current maximum number of rows in the table."},
    {"time", (getter) NodeTable_get_time, NULL, "The time array"},
    {"flags", (getter) NodeTable_get_flags, NULL, "The flags array"},
    {"population", (getter) NodeTable_get_population, NULL, "The population array"},
    {"individual", (getter) NodeTable_get_individual, NULL, "The individual array"},
    {"metadata", (getter) NodeTable_get_metadata, NULL, "The metadata array"},
    {"metadata_offset", (getter) NodeTable_get_metadata_offset, NULL,
        "The metadata offset array"},
    {NULL}  /* Sentinel */
};

static PyMethodDef NodeTable_methods[] = {
    {"add_row", (PyCFunction) NodeTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    {"equals", (PyCFunction) NodeTable_equals, METH_VARARGS,
        "Returns True if the specified NodeTable is equal to this one."},
    {"get_row", (PyCFunction) NodeTable_get_row, METH_VARARGS,
        "Returns the kth row in this table."},
    {"append_columns", (PyCFunction) NodeTable_append_columns, METH_VARARGS,
        "Appends the data in the specified arrays into the columns."},
    {"set_columns", (PyCFunction) NodeTable_set_columns, METH_VARARGS,
        "Copies the data in the specified arrays into the columns."},
    {"clear", (PyCFunction) NodeTable_clear, METH_NOARGS,
        "Clears this table."},
    {"truncate", (PyCFunction) NodeTable_truncate, METH_VARARGS,
        "Truncates this table to the specified number of rows."},
    {NULL}  /* Sentinel */
};

static PyTypeObject NodeTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_tskit.NodeTable",             /* tp_name */
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
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
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
 * EdgeTable
 *===================================================================
 */

static int
EdgeTable_check_state(EdgeTable *self)
{
    int ret = -1;
    if (self->table == NULL) {
        PyErr_SetString(PyExc_SystemError, "EdgeTable not initialised");
        goto out;
    }
    if (self->locked) {
        PyErr_SetString(PyExc_RuntimeError, "EdgeTable in use by other thread.");
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static void
EdgeTable_dealloc(EdgeTable* self)
{
    if (self->tables != NULL) {
        Py_DECREF(self->tables);
    } else if (self->table != NULL) {
        tsk_edge_tbl_free(self->table);
        PyMem_Free(self->table);
        self->table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
EdgeTable_init(EdgeTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"max_rows_increment", NULL};
    Py_ssize_t max_rows_increment = 0;

    self->table = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist, &max_rows_increment)) {
        goto out;
    }
    if (max_rows_increment < 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    self->table = PyMem_Malloc(sizeof(tsk_edge_tbl_t));
    if (self->table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tsk_edge_tbl_alloc(self->table, (size_t) max_rows_increment);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
EdgeTable_add_row(EdgeTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    double left = 0.0;
    double right = 1.0;
    int parent;
    int child;
    static char *kwlist[] = {"left", "right", "parent", "child", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddii", kwlist,
                &left, &right, &parent, &child)) {
        goto out;
    }
    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    err = tsk_edge_tbl_add_row(self->table, left, right, parent, child);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("i", err);
out:
    return ret;
}

/* Forward declaration */
static PyTypeObject EdgeTableType;

static PyObject *
EdgeTable_equals(EdgeTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    EdgeTable *other = NULL;

    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "O!", &EdgeTableType, &other)) {
        goto out;
    }
    ret = Py_BuildValue("i", tsk_edge_tbl_equals(self->table, other->table));
out:
    return ret;
}

static PyObject *
EdgeTable_get_row(EdgeTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t row_id;
    int err;
    tsk_edge_t edge;

    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &row_id)) {
        goto out;
    }
    err = tsk_edge_tbl_get_row(self->table, (size_t) row_id, &edge);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_edge(&edge);
out:
    return ret;
}

static PyObject *
EdgeTable_parse_dict_arg(EdgeTable *self, PyObject *args, bool clear_table)
{
    int err;
    PyObject *ret = NULL;
    PyObject *dict = NULL;

    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) {
        goto out;
    }
    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    err = parse_edge_table_dict(self->table, dict, clear_table);
    if (err != 0) {
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
EdgeTable_append_columns(EdgeTable *self, PyObject *args)
{
    return EdgeTable_parse_dict_arg(self, args, false);
}

static PyObject *
EdgeTable_set_columns(EdgeTable *self, PyObject *args)
{
    return EdgeTable_parse_dict_arg(self, args, true);
}

static PyObject *
EdgeTable_clear(EdgeTable *self)
{
    PyObject *ret = NULL;
    int err;

    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    err = tsk_edge_tbl_clear(self->table);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
EdgeTable_truncate(EdgeTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t num_rows;
    int err;

    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &num_rows)) {
        goto out;
    }
    if (num_rows < 0 || num_rows > (Py_ssize_t) self->table->num_rows) {
        PyErr_SetString(PyExc_ValueError, "num_rows out of bounds");
        goto out;
    }
    err = tsk_edge_tbl_truncate(self->table, (size_t) num_rows);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
EdgeTable_get_max_rows_increment(EdgeTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows_increment);
out:
    return ret;
}

static PyObject *
EdgeTable_get_num_rows(EdgeTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->num_rows);
out:
    return ret;
}

static PyObject *
EdgeTable_get_max_rows(EdgeTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows);
out:
    return ret;
}

static PyObject *
EdgeTable_get_left(EdgeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows, self->table->left, NPY_FLOAT64,
            sizeof(double));
out:
    return ret;
}

static PyObject *
EdgeTable_get_right(EdgeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows, self->table->right, NPY_FLOAT64,
            sizeof(double));
out:
    return ret;
}

static PyObject *
EdgeTable_get_parent(EdgeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows, self->table->parent, NPY_INT32,
            sizeof(int32_t));
out:
    return ret;
}

static PyObject *
EdgeTable_get_child(EdgeTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (EdgeTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows, self->table->child, NPY_INT32,
            sizeof(int32_t));
out:
    return ret;
}

static PyGetSetDef EdgeTable_getsetters[] = {
    {"max_rows_increment",
        (getter) EdgeTable_get_max_rows_increment, NULL,
        "The size increment"},
    {"num_rows", (getter) EdgeTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"max_rows", (getter) EdgeTable_get_max_rows, NULL,
        "The current maximum number of rows in the table."},
    {"left", (getter) EdgeTable_get_left, NULL, "The left array"},
    {"right", (getter) EdgeTable_get_right, NULL, "The right array"},
    {"parent", (getter) EdgeTable_get_parent, NULL, "The parent array"},
    {"child", (getter) EdgeTable_get_child, NULL, "The child array"},
    {NULL}  /* Sentinel */
};

static PyMethodDef EdgeTable_methods[] = {
    {"add_row", (PyCFunction) EdgeTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    {"equals", (PyCFunction) EdgeTable_equals, METH_VARARGS,
        "Returns True if the specified EdgeTable is equal to this one."},
    {"get_row", (PyCFunction) EdgeTable_get_row, METH_VARARGS,
        "Returns the kth row in this table."},
    {"set_columns", (PyCFunction) EdgeTable_set_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the specified arrays into the columns."},
    {"append_columns", (PyCFunction) EdgeTable_append_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the specified arrays into the columns."},
    {"clear", (PyCFunction) EdgeTable_clear, METH_NOARGS,
        "Clears this table."},
    {"truncate", (PyCFunction) EdgeTable_truncate, METH_VARARGS,
        "Truncates this table to the specified number of rows."},
    {NULL}  /* Sentinel */
};

static PyTypeObject EdgeTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_tskit.EdgeTable",             /* tp_name */
    sizeof(EdgeTable),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)EdgeTable_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "EdgeTable objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    EdgeTable_methods,             /* tp_methods */
    0,                             /* tp_members */
    EdgeTable_getsetters,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)EdgeTable_init,      /* tp_init */
};

/*===================================================================
 * MigrationTable
 *===================================================================
 */

static int
MigrationTable_check_state(MigrationTable *self)
{
    int ret = -1;
    if (self->table == NULL) {
        PyErr_SetString(PyExc_SystemError, "MigrationTable not initialised");
        goto out;
    }
    if (self->locked) {
        PyErr_SetString(PyExc_RuntimeError, "MigrationTable in use by other thread.");
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static void
MigrationTable_dealloc(MigrationTable* self)
{
    if (self->tables != NULL) {
        Py_DECREF(self->tables);
    } else if (self->table != NULL) {
        tsk_migration_tbl_free(self->table);
        PyMem_Free(self->table);
        self->table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
MigrationTable_init(MigrationTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"max_rows_increment", NULL};
    Py_ssize_t max_rows_increment = 0;

    self->table = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist,
                &max_rows_increment)) {
        goto out;
    }
    if (max_rows_increment < 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    self->table = PyMem_Malloc(sizeof(tsk_migration_tbl_t));
    if (self->table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tsk_migration_tbl_alloc(self->table, (size_t) max_rows_increment);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
MigrationTable_add_row(MigrationTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    double left, right, time;
    int node, source, dest;
    static char *kwlist[] = {"left", "right", "node", "source", "dest", "time", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddiiid", kwlist,
                &left, &right, &node, &source, &dest, &time)) {
        goto out;
    }
    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    err = tsk_migration_tbl_add_row(self->table, left, right, node,
            source, dest, time);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("i", err);
out:
    return ret;
}

/* Forward declaration */
static PyTypeObject MigrationTableType;

static PyObject *
MigrationTable_equals(MigrationTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    MigrationTable *other = NULL;

    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "O!", &MigrationTableType, &other)) {
        goto out;
    }
    ret = Py_BuildValue("i", tsk_migration_tbl_equals(self->table, other->table));
out:
    return ret;
}

static PyObject *
MigrationTable_get_row(MigrationTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t row_id;
    int err;
    tsk_migration_t migration;

    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &row_id)) {
        goto out;
    }
    err = tsk_migration_tbl_get_row(self->table, (size_t) row_id, &migration);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_migration(&migration);
out:
    return ret;
}

static PyObject *
MigrationTable_parse_dict_arg(MigrationTable *self, PyObject *args, bool clear_table)
{
    int err;
    PyObject *ret = NULL;
    PyObject *dict = NULL;

    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) {
        goto out;
    }
    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    err = parse_migration_table_dict(self->table, dict, clear_table);
    if (err != 0) {
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
MigrationTable_append_columns(MigrationTable *self, PyObject *args)
{
    return MigrationTable_parse_dict_arg(self, args, false);
}

static PyObject *
MigrationTable_set_columns(MigrationTable *self, PyObject *args)
{
    return MigrationTable_parse_dict_arg(self, args, true);
}

static PyObject *
MigrationTable_clear(MigrationTable *self)
{
    PyObject *ret = NULL;
    int err;

    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    err = tsk_migration_tbl_clear(self->table);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
MigrationTable_truncate(MigrationTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t num_rows;
    int err;

    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &num_rows)) {
        goto out;
    }
    if (num_rows < 0 || num_rows > (Py_ssize_t) self->table->num_rows) {
        PyErr_SetString(PyExc_ValueError, "num_rows out of bounds");
        goto out;
    }
    err = tsk_migration_tbl_truncate(self->table, (size_t) num_rows);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
MigrationTable_get_max_rows_increment(MigrationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows_increment);
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
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->num_rows);
out:
    return ret;
}

static PyObject *
MigrationTable_get_max_rows(MigrationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (MigrationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows);
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
    ret = table_get_column_array(self->table->num_rows, self->table->left,
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
    ret = table_get_column_array(self->table->num_rows, self->table->right,
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
    ret = table_get_column_array(self->table->num_rows, self->table->time,
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
    ret = table_get_column_array(self->table->num_rows, self->table->node,
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
    ret = table_get_column_array(self->table->num_rows, self->table->source,
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
    ret = table_get_column_array(self->table->num_rows, self->table->dest,
            NPY_INT32, sizeof(uint32_t));
out:
    return ret;
}

static PyGetSetDef MigrationTable_getsetters[] = {
    {"max_rows_increment",
        (getter) MigrationTable_get_max_rows_increment, NULL, "The size increment"},
    {"num_rows", (getter) MigrationTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"max_rows", (getter) MigrationTable_get_max_rows, NULL,
        "The current maximum number of rows in the table."},
    {"left", (getter) MigrationTable_get_left, NULL, "The left array"},
    {"right", (getter) MigrationTable_get_right, NULL, "The right array"},
    {"node", (getter) MigrationTable_get_node, NULL, "The node array"},
    {"source", (getter) MigrationTable_get_source, NULL, "The source array"},
    {"dest", (getter) MigrationTable_get_dest, NULL, "The dest array"},
    {"time", (getter) MigrationTable_get_time, NULL, "The time array"},
    {NULL}  /* Sentinel */
};

static PyMethodDef MigrationTable_methods[] = {
    {"add_row", (PyCFunction) MigrationTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    {"equals", (PyCFunction) MigrationTable_equals, METH_VARARGS,
        "Returns True if the specified MigrationTable is equal to this one."},
    {"get_row", (PyCFunction) MigrationTable_get_row, METH_VARARGS,
        "Returns the kth row in this table."},
    {"set_columns", (PyCFunction) MigrationTable_set_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the specified arrays into the columns."},
    {"append_columns", (PyCFunction) MigrationTable_append_columns, METH_VARARGS|METH_KEYWORDS,
        "Appends the data in the specified arrays into the columns."},
    {"clear", (PyCFunction) MigrationTable_clear, METH_NOARGS,
        "Clears this table."},
    {"truncate", (PyCFunction) MigrationTable_truncate, METH_VARARGS,
        "Truncates this table to the specified number of rows."},
    {NULL}  /* Sentinel */
};

static PyTypeObject MigrationTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_tskit.MigrationTable",             /* tp_name */
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
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
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
 * SiteTable
 *===================================================================
 */

static int
SiteTable_check_state(SiteTable *self)
{
    int ret = -1;
    if (self->table == NULL) {
        PyErr_SetString(PyExc_SystemError, "SiteTable not initialised");
        goto out;
    }
    if (self->locked) {
        PyErr_SetString(PyExc_RuntimeError, "SiteTable in use by other thread.");
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static void
SiteTable_dealloc(SiteTable* self)
{
    if (self->tables != NULL) {
        Py_DECREF(self->tables);
    } else if (self->table != NULL) {
        tsk_site_tbl_free(self->table);
        PyMem_Free(self->table);
        self->table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
SiteTable_init(SiteTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"max_rows_increment", NULL};
    Py_ssize_t max_rows_increment = 0;

    self->table = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist, &max_rows_increment)) {
        goto out;
    }
    if (max_rows_increment < 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    self->table = PyMem_Malloc(sizeof(tsk_site_tbl_t));
    if (self->table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tsk_site_tbl_alloc(self->table, (size_t) max_rows_increment, 0, 0);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
SiteTable_add_row(SiteTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    double position;
    char *ancestral_state = NULL;
    Py_ssize_t ancestral_state_length = 0;
    PyObject *py_metadata = Py_None;
    char *metadata = NULL;
    Py_ssize_t metadata_length = 0;
    static char *kwlist[] = {"position", "ancestral_state", "metadata", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ds#|O", kwlist,
                &position, &ancestral_state, &ancestral_state_length, &py_metadata)) {
        goto out;
    }
    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    if (py_metadata != Py_None) {
        if (PyBytes_AsStringAndSize(py_metadata, &metadata, &metadata_length) < 0) {
            goto out;
        }
    }
    err = tsk_site_tbl_add_row(self->table, position, ancestral_state,
            ancestral_state_length, metadata, metadata_length);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("i", err);
out:
    return ret;
}

/* Forward declaration */
static PyTypeObject SiteTableType;

static PyObject *
SiteTable_equals(SiteTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    SiteTable *other = NULL;

    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "O!", &SiteTableType, &other)) {
        goto out;
    }
    ret = Py_BuildValue("i", tsk_site_tbl_equals(self->table, other->table));
out:
    return ret;
}

static PyObject *
SiteTable_get_row(SiteTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t row_id;
    int err;
    tsk_site_t site;

    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &row_id)) {
        goto out;
    }
    err = tsk_site_tbl_get_row(self->table, (size_t) row_id, &site);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_site_row(&site);
out:
    return ret;
}

static PyObject *
SiteTable_parse_dict_arg(SiteTable *self, PyObject *args, bool clear_table)
{
    int err;
    PyObject *ret = NULL;
    PyObject *dict = NULL;

    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) {
        goto out;
    }
    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    err = parse_site_table_dict(self->table, dict, clear_table);
    if (err != 0) {
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
SiteTable_append_columns(SiteTable *self, PyObject *args)
{
    return SiteTable_parse_dict_arg(self, args, false);
}

static PyObject *
SiteTable_set_columns(SiteTable *self, PyObject *args)
{
    return SiteTable_parse_dict_arg(self, args, true);
}

static PyObject *
SiteTable_clear(SiteTable *self)
{
    PyObject *ret = NULL;
    int err;

    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    err = tsk_site_tbl_clear(self->table);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
SiteTable_truncate(SiteTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t num_rows;
    int err;

    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &num_rows)) {
        goto out;
    }
    if (num_rows < 0 || num_rows > (Py_ssize_t) self->table->num_rows) {
        PyErr_SetString(PyExc_ValueError, "num_rows out of bounds");
        goto out;
    }
    err = tsk_site_tbl_truncate(self->table, (size_t) num_rows);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
SiteTable_get_max_rows_increment(SiteTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows_increment);
out:
    return ret;
}

static PyObject *
SiteTable_get_num_rows(SiteTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->num_rows);
out:
    return ret;
}

static PyObject *
SiteTable_get_max_rows(SiteTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows);
out:
    return ret;
}

static PyObject *
SiteTable_get_position(SiteTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows,
            self->table->position, NPY_FLOAT64, sizeof(double));
out:
    return ret;
}

static PyObject *
SiteTable_get_ancestral_state(SiteTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->ancestral_state_length,
            self->table->ancestral_state, NPY_INT8, sizeof(char));
out:
    return ret;
}

static PyObject *
SiteTable_get_ancestral_state_offset(SiteTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows + 1,
            self->table->ancestral_state_offset, NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyObject *
SiteTable_get_metadata(SiteTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->metadata_length,
            self->table->metadata, NPY_INT8, sizeof(char));
out:
    return ret;
}

static PyObject *
SiteTable_get_metadata_offset(SiteTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (SiteTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows + 1,
            self->table->metadata_offset, NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyGetSetDef SiteTable_getsetters[] = {
    {"max_rows_increment",
        (getter) SiteTable_get_max_rows_increment, NULL,
        "The size increment"},
    {"num_rows",
        (getter) SiteTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"max_rows",
        (getter) SiteTable_get_max_rows, NULL,
        "The current maximum number of rows in the table."},
    {"position", (getter) SiteTable_get_position, NULL,
        "The position array."},
    {"ancestral_state", (getter) SiteTable_get_ancestral_state, NULL,
        "The ancestral state array."},
    {"ancestral_state_offset", (getter) SiteTable_get_ancestral_state_offset, NULL,
        "The ancestral state offset array."},
    {"metadata", (getter) SiteTable_get_metadata, NULL,
        "The metadata array."},
    {"metadata_offset", (getter) SiteTable_get_metadata_offset, NULL,
        "The metadata offset array."},
    {NULL}  /* Sentinel */
};

static PyMethodDef SiteTable_methods[] = {
    {"add_row", (PyCFunction) SiteTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    {"equals", (PyCFunction) SiteTable_equals, METH_VARARGS,
        "Returns True if the specified SiteTable is equal to this one."},
    {"get_row", (PyCFunction) SiteTable_get_row, METH_VARARGS,
        "Returns the kth row in this table."},
    {"set_columns", (PyCFunction) SiteTable_set_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the specified arrays into the columns."},
    {"append_columns", (PyCFunction) SiteTable_append_columns, METH_VARARGS|METH_KEYWORDS,
        "Appends the data in the specified arrays into the columns."},
    {"clear", (PyCFunction) SiteTable_clear, METH_NOARGS,
        "Clears this table."},
    {"truncate", (PyCFunction) SiteTable_truncate, METH_VARARGS,
        "Truncates this table to the specified number of rows."},
    {NULL}  /* Sentinel */
};

static PyTypeObject SiteTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_tskit.SiteTable",             /* tp_name */
    sizeof(SiteTable),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)SiteTable_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "SiteTable objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    SiteTable_methods,             /* tp_methods */
    0,                             /* tp_members */
    SiteTable_getsetters,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)SiteTable_init,      /* tp_init */
};


/*===================================================================
 * MutationTable
 *===================================================================
 */

static int
MutationTable_check_state(MutationTable *self)
{
    int ret = -1;
    if (self->table == NULL) {
        PyErr_SetString(PyExc_SystemError, "MutationTable not initialised");
        goto out;
    }
    if (self->locked) {
        PyErr_SetString(PyExc_RuntimeError, "MutationTable in use by other thread.");
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static void
MutationTable_dealloc(MutationTable* self)
{
    if (self->tables != NULL) {
        Py_DECREF(self->tables);
    } else if (self->table != NULL) {
        tsk_mutation_tbl_free(self->table);
        PyMem_Free(self->table);
        self->table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
MutationTable_init(MutationTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"max_rows_increment", NULL};
    Py_ssize_t max_rows_increment = 0;

    self->table = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist, &max_rows_increment)) {
        goto out;
    }
    if (max_rows_increment < 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    self->table = PyMem_Malloc(sizeof(tsk_mutation_tbl_t));
    if (self->table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tsk_mutation_tbl_alloc(self->table, (size_t) max_rows_increment, 0, 0);
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
    int site;
    int node;
    int parent = TSK_NULL_MUTATION;
    char *derived_state;
    Py_ssize_t derived_state_length;
    PyObject *py_metadata = Py_None;
    char *metadata = NULL;
    Py_ssize_t metadata_length = 0;
    static char *kwlist[] = {"site", "node", "derived_state", "parent", "metadata", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "iis#|iO", kwlist,
                &site, &node, &derived_state, &derived_state_length, &parent,
                &py_metadata)) {
        goto out;
    }
    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    if (py_metadata != Py_None) {
        if (PyBytes_AsStringAndSize(py_metadata, &metadata, &metadata_length) < 0) {
            goto out;
        }
    }
    err = tsk_mutation_tbl_add_row(self->table, (tsk_id_t) site,
            (tsk_id_t) node, (tsk_id_t) parent,
            derived_state, derived_state_length,
            metadata, metadata_length);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("i", err);
out:
    return ret;
}

/* Forward declaration */
static PyTypeObject MutationTableType;

static PyObject *
MutationTable_equals(MutationTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    MutationTable *other = NULL;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "O!", &MutationTableType, &other)) {
        goto out;
    }
    ret = Py_BuildValue("i", tsk_mutation_tbl_equals(self->table, other->table));
out:
    return ret;
}

static PyObject *
MutationTable_get_row(MutationTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t row_id;
    int err;
    tsk_mutation_t mutation;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &row_id)) {
        goto out;
    }
    err = tsk_mutation_tbl_get_row(self->table, (size_t) row_id, &mutation);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_mutation(&mutation);
out:
    return ret;
}

static PyObject *
MutationTable_parse_dict_arg(MutationTable *self, PyObject *args, bool clear_table)
{
    int err;
    PyObject *ret = NULL;
    PyObject *dict = NULL;

    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) {
        goto out;
    }
    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    err = parse_mutation_table_dict(self->table, dict, clear_table);
    if (err != 0) {
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
MutationTable_append_columns(MutationTable *self, PyObject *args)
{
    return MutationTable_parse_dict_arg(self, args, false);
}

static PyObject *
MutationTable_set_columns(MutationTable *self, PyObject *args)
{
    return MutationTable_parse_dict_arg(self, args, true);
}

static PyObject *
MutationTable_clear(MutationTable *self)
{
    PyObject *ret = NULL;
    int err;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    err = tsk_mutation_tbl_clear(self->table);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
MutationTable_truncate(MutationTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t num_rows;
    int err;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &num_rows)) {
        goto out;
    }
    if (num_rows < 0 || num_rows > (Py_ssize_t) self->table->num_rows) {
        PyErr_SetString(PyExc_ValueError, "num_rows out of bounds");
        goto out;
    }
    err = tsk_mutation_tbl_truncate(self->table, (size_t) num_rows);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
MutationTable_get_max_rows_increment(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows_increment);
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
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->num_rows);
out:
    return ret;
}

static PyObject *
MutationTable_get_max_rows(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows);
out:
    return ret;
}

static PyObject *
MutationTable_get_site(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows, self->table->site, NPY_INT32,
            sizeof(int32_t));
out:
    return ret;
}

static PyObject *
MutationTable_get_node(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows, self->table->node, NPY_INT32,
            sizeof(int32_t));
out:
    return ret;
}

static PyObject *
MutationTable_get_parent(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows, self->table->parent, NPY_INT32,
            sizeof(int32_t));
out:
    return ret;
}

static PyObject *
MutationTable_get_derived_state(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->derived_state_length, self->table->derived_state,
            NPY_INT8, sizeof(char));
out:
    return ret;
}

static PyObject *
MutationTable_get_derived_state_offset(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows + 1, self->table->derived_state_offset,
            NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyObject *
MutationTable_get_metadata(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->metadata_length, self->table->metadata,
            NPY_INT8, sizeof(char));
out:
    return ret;
}

static PyObject *
MutationTable_get_metadata_offset(MutationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (MutationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(
            self->table->num_rows + 1, self->table->metadata_offset,
            NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyGetSetDef MutationTable_getsetters[] = {
    {"max_rows_increment",
        (getter) MutationTable_get_max_rows_increment, NULL,
        "The size increment"},
    {"num_rows",
        (getter) MutationTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"max_rows",
        (getter) MutationTable_get_max_rows, NULL,
        "The curret maximum number of rows in the table."},
    {"site", (getter) MutationTable_get_site, NULL, "The site array"},
    {"node", (getter) MutationTable_get_node, NULL, "The node array"},
    {"parent", (getter) MutationTable_get_parent, NULL, "The parent array"},
    {"derived_state", (getter) MutationTable_get_derived_state, NULL,
        "The derived_state array"},
    {"derived_state_offset", (getter) MutationTable_get_derived_state_offset, NULL,
        "The derived_state_offset array"},
    {"metadata", (getter) MutationTable_get_metadata, NULL,
        "The metadata array"},
    {"metadata_offset", (getter) MutationTable_get_metadata_offset, NULL,
        "The metadata_offset array"},
    {NULL}  /* Sentinel */
};

static PyMethodDef MutationTable_methods[] = {
    {"add_row", (PyCFunction) MutationTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    {"equals", (PyCFunction) MutationTable_equals, METH_VARARGS,
        "Returns True if the specified MutationTable is equal to this one."},
    {"get_row", (PyCFunction) MutationTable_get_row, METH_VARARGS,
        "Returns the kth row in this table."},
    {"set_columns", (PyCFunction) MutationTable_set_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the specified arrays into the columns."},
    {"append_columns", (PyCFunction) MutationTable_append_columns, METH_VARARGS|METH_KEYWORDS,
        "Appends the data in the specified  arrays into the columns."},
    {"clear", (PyCFunction) MutationTable_clear, METH_NOARGS,
        "Clears this table."},
    {"truncate", (PyCFunction) MutationTable_truncate, METH_VARARGS,
        "Truncates this table to the specified number of rows."},
    {NULL}  /* Sentinel */
};

static PyTypeObject MutationTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_tskit.MutationTable",             /* tp_name */
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
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
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
 * PopulationTable
 *===================================================================
 */

static int
PopulationTable_check_state(PopulationTable *self)
{
    int ret = -1;
    if (self->table == NULL) {
        PyErr_SetString(PyExc_SystemError, "PopulationTable not initialised");
        goto out;
    }
    if (self->locked) {
        PyErr_SetString(PyExc_RuntimeError, "PopulationTable in use by other thread.");
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static void
PopulationTable_dealloc(PopulationTable* self)
{
    if (self->tables != NULL) {
        Py_DECREF(self->tables);
    } else if (self->table != NULL) {
        tsk_population_tbl_free(self->table);
        PyMem_Free(self->table);
        self->table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
PopulationTable_init(PopulationTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"max_rows_increment", NULL};
    Py_ssize_t max_rows_increment = 0;

    self->table = NULL;
    self->locked = false;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist,
                &max_rows_increment)) {
        goto out;
    }
    if (max_rows_increment < 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    self->table = PyMem_Malloc(sizeof(tsk_population_tbl_t));
    if (self->table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    /* Take the default size increments for metadata and record */
    err = tsk_population_tbl_alloc(self->table,
            (size_t) max_rows_increment, 0);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
PopulationTable_add_row(PopulationTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    PyObject *py_metadata = Py_None;
    char *metadata = NULL;
    Py_ssize_t metadata_length = 0;
    static char *kwlist[] = {"metadata", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &py_metadata)) {
        goto out;
    }
    if (PopulationTable_check_state(self) != 0) {
        goto out;
    }

    if (py_metadata != Py_None) {
        if (PyBytes_AsStringAndSize(py_metadata, &metadata, &metadata_length) < 0) {
            goto out;
        }
    }
    err = tsk_population_tbl_add_row(self->table, metadata, metadata_length);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("i", err);
out:
    return ret;
}

/* Forward declaration */
static PyTypeObject PopulationTableType;

static PyObject *
PopulationTable_equals(PopulationTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    PopulationTable *other = NULL;

    if (PopulationTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "O!", &PopulationTableType, &other)) {
        goto out;
    }
    ret = Py_BuildValue("i", tsk_population_tbl_equals(self->table, other->table));
out:
    return ret;
}

static PyObject *
PopulationTable_get_row(PopulationTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t row_id;
    int err;
    tsk_population_t population;

    if (PopulationTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &row_id)) {
        goto out;
    }
    err = tsk_population_tbl_get_row(self->table, (size_t) row_id, &population);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_population(&population);
out:
    return ret;
}

static PyObject *
PopulationTable_parse_dict_arg(PopulationTable *self, PyObject *args, bool clear_table)
{
    int err;
    PyObject *ret = NULL;
    PyObject *dict = NULL;

    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) {
        goto out;
    }
    if (PopulationTable_check_state(self) != 0) {
        goto out;
    }
    err = parse_population_table_dict(self->table, dict, clear_table);
    if (err != 0) {
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
PopulationTable_append_columns(PopulationTable *self, PyObject *args)
{
    return PopulationTable_parse_dict_arg(self, args, false);
}

static PyObject *
PopulationTable_set_columns(PopulationTable *self, PyObject *args)
{
    return PopulationTable_parse_dict_arg(self, args, true);
}

static PyObject *
PopulationTable_clear(PopulationTable *self)
{
    PyObject *ret = NULL;
    int err;

    if (PopulationTable_check_state(self) != 0) {
        goto out;
    }
    err = tsk_population_tbl_clear(self->table);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
PopulationTable_truncate(PopulationTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t num_rows;
    int err;

    if (PopulationTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &num_rows)) {
        goto out;
    }
    if (num_rows < 0 || num_rows > (Py_ssize_t) self->table->num_rows) {
        PyErr_SetString(PyExc_ValueError, "num_rows out of bounds");
        goto out;
    }
    err = tsk_population_tbl_truncate(self->table, (size_t) num_rows);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
PopulationTable_get_max_rows_increment(PopulationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (PopulationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows_increment);
out:
    return ret;
}

static PyObject *
PopulationTable_get_num_rows(PopulationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (PopulationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->num_rows);
out:
    return ret;
}

static PyObject *
PopulationTable_get_max_rows(PopulationTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (PopulationTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows);
out:
    return ret;
}

static PyObject *
PopulationTable_get_metadata(PopulationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (PopulationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->metadata_length,
            self->table->metadata, NPY_INT8, sizeof(char));
out:
    return ret;
}

static PyObject *
PopulationTable_get_metadata_offset(PopulationTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (PopulationTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->num_rows + 1,
            self->table->metadata_offset, NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyGetSetDef PopulationTable_getsetters[] = {
    {"max_rows_increment",
        (getter) PopulationTable_get_max_rows_increment, NULL, "The size increment"},
    {"num_rows", (getter) PopulationTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"max_rows", (getter) PopulationTable_get_max_rows, NULL,
        "The current maximum number of rows in the table."},
    {"metadata", (getter) PopulationTable_get_metadata, NULL, "The metadata array"},
    {"metadata_offset", (getter) PopulationTable_get_metadata_offset, NULL,
        "The metadata offset array"},
    {NULL}  /* Sentinel */
};

static PyMethodDef PopulationTable_methods[] = {
    {"add_row", (PyCFunction) PopulationTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    {"equals", (PyCFunction) PopulationTable_equals, METH_VARARGS,
        "Returns True if the specified PopulationTable is equal to this one."},
    {"get_row", (PyCFunction) PopulationTable_get_row, METH_VARARGS,
        "Returns the kth row in this table."},
    {"append_columns", (PyCFunction) PopulationTable_append_columns,
        METH_VARARGS|METH_KEYWORDS,
        "Appends the data in the specified arrays into the columns."},
    {"set_columns", (PyCFunction) PopulationTable_set_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the specified arrays into the columns."},
    {"clear", (PyCFunction) PopulationTable_clear, METH_NOARGS,
        "Clears this table."},
    {"truncate", (PyCFunction) PopulationTable_truncate, METH_VARARGS,
        "Truncates this table to the specified number of rows."},
    {NULL}  /* Sentinel */
};

static PyTypeObject PopulationTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_tskit.PopulationTable",             /* tp_name */
    sizeof(PopulationTable),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)PopulationTable_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "PopulationTable objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    PopulationTable_methods,             /* tp_methods */
    0,                             /* tp_members */
    PopulationTable_getsetters,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PopulationTable_init,      /* tp_init */
};


/*===================================================================
 * ProvenanceTable
 *===================================================================
 */

static int
ProvenanceTable_check_state(ProvenanceTable *self)
{
    int ret = -1;
    if (self->table == NULL) {
        PyErr_SetString(PyExc_SystemError, "ProvenanceTable not initialised");
        goto out;
    }
    if (self->locked) {
        PyErr_SetString(PyExc_RuntimeError, "ProvenanceTable in use by other thread.");
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static void
ProvenanceTable_dealloc(ProvenanceTable* self)
{
    if (self->tables != NULL) {
        Py_DECREF(self->tables);
    } else if (self->table != NULL) {
        tsk_provenance_tbl_free(self->table);
        PyMem_Free(self->table);
        self->table = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
ProvenanceTable_init(ProvenanceTable *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"max_rows_increment", NULL};
    Py_ssize_t max_rows_increment = 0;

    self->table = NULL;
    self->locked = false;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist,
                &max_rows_increment)) {
        goto out;
    }
    if (max_rows_increment < 0) {
        PyErr_SetString(PyExc_ValueError, "max_rows_increment must be positive");
        goto out;
    }
    self->table = PyMem_Malloc(sizeof(tsk_provenance_tbl_t));
    if (self->table == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    /* Take the default size increments for timestamp and record */
    err = tsk_provenance_tbl_alloc(self->table,
            (size_t) max_rows_increment, 0, 0);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    return ret;
}
static PyObject *
ProvenanceTable_add_row(ProvenanceTable *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    int err;
    char *timestamp = "";
    Py_ssize_t timestamp_length = 0;
    char *record = "";
    Py_ssize_t record_length = 0;
    static char *kwlist[] = {"timestamp", "record", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s#s#", kwlist,
                &timestamp, &timestamp_length, &record, &record_length)){
        goto out;
    }
    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    err = tsk_provenance_tbl_add_row(self->table,
            timestamp, timestamp_length, record, record_length);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("i", err);
out:
    return ret;
}

/* Forward declaration */
static PyTypeObject ProvenanceTableType;

static PyObject *
ProvenanceTable_equals(ProvenanceTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    ProvenanceTable *other = NULL;

    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "O!", &ProvenanceTableType, &other)) {
        goto out;
    }
    ret = Py_BuildValue("i", tsk_provenance_tbl_equals(self->table, other->table));
out:
    return ret;
}

static PyObject *
ProvenanceTable_get_row(ProvenanceTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t row_id;
    int err;
    tsk_provenance_t provenance;

    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &row_id)) {
        goto out;
    }
    err = tsk_provenance_tbl_get_row(self->table, (size_t) row_id, &provenance);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_provenance(&provenance);
out:
    return ret;
}

static PyObject *
ProvenanceTable_parse_dict_arg(ProvenanceTable *self, PyObject *args, bool clear_table)
{
    int err;
    PyObject *ret = NULL;
    PyObject *dict = NULL;

    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) {
        goto out;
    }
    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    err = parse_provenance_table_dict(self->table, dict, clear_table);
    if (err != 0) {
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
ProvenanceTable_append_columns(ProvenanceTable *self, PyObject *args)
{
    return ProvenanceTable_parse_dict_arg(self, args, false);
}

static PyObject *
ProvenanceTable_set_columns(ProvenanceTable *self, PyObject *args)
{
    return ProvenanceTable_parse_dict_arg(self, args, true);
}

static PyObject *
ProvenanceTable_clear(ProvenanceTable *self)
{
    PyObject *ret = NULL;
    int err;

    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    err = tsk_provenance_tbl_clear(self->table);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
ProvenanceTable_truncate(ProvenanceTable *self, PyObject *args)
{
    PyObject *ret = NULL;
    Py_ssize_t num_rows;
    int err;

    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &num_rows)) {
        goto out;
    }
    if (num_rows < 0 || num_rows > (Py_ssize_t) self->table->num_rows) {
        PyErr_SetString(PyExc_ValueError, "num_rows out of bounds");
        goto out;
    }
    err = tsk_provenance_tbl_truncate(self->table, (size_t) num_rows);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
ProvenanceTable_get_max_rows_increment(ProvenanceTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows_increment);
out:
    return ret;
}

static PyObject *
ProvenanceTable_get_num_rows(ProvenanceTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->num_rows);
out:
    return ret;
}

static PyObject *
ProvenanceTable_get_max_rows(ProvenanceTable *self, void *closure)
{
    PyObject *ret = NULL;
    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->table->max_rows);
out:
    return ret;
}

static PyObject *
ProvenanceTable_get_timestamp(ProvenanceTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->timestamp_length,
            self->table->timestamp, NPY_INT8, sizeof(char));
out:
    return ret;
}

static PyObject *
ProvenanceTable_get_timestamp_offset(ProvenanceTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->num_rows + 1,
            self->table->timestamp_offset, NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyObject *
ProvenanceTable_get_record(ProvenanceTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->record_length,
            self->table->record, NPY_INT8, sizeof(char));
out:
    return ret;
}

static PyObject *
ProvenanceTable_get_record_offset(ProvenanceTable *self, void *closure)
{
    PyObject *ret = NULL;

    if (ProvenanceTable_check_state(self) != 0) {
        goto out;
    }
    ret = table_get_column_array(self->table->num_rows + 1,
            self->table->record_offset, NPY_UINT32, sizeof(uint32_t));
out:
    return ret;
}

static PyGetSetDef ProvenanceTable_getsetters[] = {
    {"max_rows_increment",
        (getter) ProvenanceTable_get_max_rows_increment, NULL, "The size increment"},
    {"num_rows", (getter) ProvenanceTable_get_num_rows, NULL,
        "The number of rows in the table."},
    {"max_rows", (getter) ProvenanceTable_get_max_rows, NULL,
        "The current maximum number of rows in the table."},
    {"timestamp", (getter) ProvenanceTable_get_timestamp, NULL, "The timestamp array"},
    {"timestamp_offset", (getter) ProvenanceTable_get_timestamp_offset, NULL,
        "The timestamp offset array"},
    {"record", (getter) ProvenanceTable_get_record, NULL, "The record array"},
    {"record_offset", (getter) ProvenanceTable_get_record_offset, NULL,
        "The record offset array"},
    {NULL}  /* Sentinel */
};

static PyMethodDef ProvenanceTable_methods[] = {
    {"add_row", (PyCFunction) ProvenanceTable_add_row, METH_VARARGS|METH_KEYWORDS,
        "Adds a new row to this table."},
    {"equals", (PyCFunction) ProvenanceTable_equals, METH_VARARGS,
        "Returns True if the specified ProvenanceTable is equal to this one."},
    {"get_row", (PyCFunction) ProvenanceTable_get_row, METH_VARARGS,
        "Returns the kth row in this table."},
    {"append_columns", (PyCFunction) ProvenanceTable_append_columns,
        METH_VARARGS|METH_KEYWORDS,
        "Appends the data in the specified arrays into the columns."},
    {"set_columns", (PyCFunction) ProvenanceTable_set_columns, METH_VARARGS|METH_KEYWORDS,
        "Copies the data in the specified arrays into the columns."},
    {"clear", (PyCFunction) ProvenanceTable_clear, METH_NOARGS,
        "Clears this table."},
    {"truncate", (PyCFunction) ProvenanceTable_truncate, METH_VARARGS,
        "Truncates this table to the specified number of rows."},
    {NULL}  /* Sentinel */
};

static PyTypeObject ProvenanceTableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_tskit.ProvenanceTable",             /* tp_name */
    sizeof(ProvenanceTable),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)ProvenanceTable_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "ProvenanceTable objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    ProvenanceTable_methods,             /* tp_methods */
    0,                             /* tp_members */
    ProvenanceTable_getsetters,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)ProvenanceTable_init,      /* tp_init */
};

/*===================================================================
 * TableCollection
 *===================================================================
 */

static void
TableCollection_dealloc(TableCollection* self)
{
    if (self->tables != NULL) {
        tsk_tbl_collection_free(self->tables);
        PyMem_Free(self->tables);
        self->tables = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
TableCollection_init(TableCollection *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"sequence_length", NULL};
    double sequence_length = -1;

    self->tables = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|d", kwlist, &sequence_length)) {
        goto out;
    }

    self->tables = PyMem_Malloc(sizeof(tsk_tbl_collection_t));
    if (self->tables == NULL) {
        PyErr_NoMemory();
    }
    err = tsk_tbl_collection_alloc(self->tables, TSK_ALLOC_TABLES);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    self->tables->sequence_length = sequence_length;
    ret = 0;
out:
    return ret;
}

/* The getters for each of the tables returns a new reference which we
 * set up here. These references use a pointer to the table stored in
 * the table collection, so to guard against this memory getting freed
 * we the Python Table classes keep a reference to the TableCollection
 * and INCREF it. We don't keep permanent references to the Table classes
 * in the TableCollection as this gives a circular references which would
 * require implementing support for cyclic garbage collection.
 */

static PyObject *
TableCollection_get_individuals(TableCollection *self, void *closure)
{
    IndividualTable *individuals = NULL;

    individuals = PyObject_New(IndividualTable, &IndividualTableType);
    if (individuals == NULL) {
        goto out;
    }
    individuals->table = self->tables->individuals;
    individuals->locked = false;
    individuals->tables = self;
    Py_INCREF(self);
out:
    return (PyObject *) individuals;
}

static PyObject *
TableCollection_get_nodes(TableCollection *self, void *closure)
{
    NodeTable *nodes = NULL;

    nodes = PyObject_New(NodeTable, &NodeTableType);
    if (nodes == NULL) {
        goto out;
    }
    nodes->table = self->tables->nodes;
    nodes->locked = false;
    nodes->tables = self;
    Py_INCREF(self);
out:
    return (PyObject *) nodes;
}

static PyObject *
TableCollection_get_edges(TableCollection *self, void *closure)
{
    EdgeTable *edges = NULL;

    edges = PyObject_New(EdgeTable, &EdgeTableType);
    if (edges == NULL) {
        goto out;
    }
    edges->table = self->tables->edges;
    edges->locked = false;
    edges->tables = self;
    Py_INCREF(self);
out:
    return (PyObject *) edges;
}

static PyObject *
TableCollection_get_migrations(TableCollection *self, void *closure)
{
    MigrationTable *migrations = NULL;

    migrations = PyObject_New(MigrationTable, &MigrationTableType);
    if (migrations == NULL) {
        goto out;
    }
    migrations->table = self->tables->migrations;
    migrations->locked = false;
    migrations->tables = self;
    Py_INCREF(self);
out:
    return (PyObject *) migrations;
}

static PyObject *
TableCollection_get_sites(TableCollection *self, void *closure)
{
    SiteTable *sites = NULL;

    sites = PyObject_New(SiteTable, &SiteTableType);
    if (sites == NULL) {
        goto out;
    }
    sites->table = self->tables->sites;
    sites->locked = false;
    sites->tables = self;
    Py_INCREF(self);
out:
    return (PyObject *) sites;
}

static PyObject *
TableCollection_get_mutations(TableCollection *self, void *closure)
{
    MutationTable *mutations = NULL;

    mutations = PyObject_New(MutationTable, &MutationTableType);
    if (mutations == NULL) {
        goto out;
    }
    mutations->table = self->tables->mutations;
    mutations->locked = false;
    mutations->tables = self;
    Py_INCREF(self);
out:
    return (PyObject *) mutations;
}

static PyObject *
TableCollection_get_populations(TableCollection *self, void *closure)
{
    PopulationTable *populations = NULL;

    populations = PyObject_New(PopulationTable, &PopulationTableType);
    if (populations == NULL) {
        goto out;
    }
    populations->table = self->tables->populations;
    populations->locked = false;
    populations->tables = self;
    Py_INCREF(self);
out:
    return (PyObject *) populations;
}

static PyObject *
TableCollection_get_provenances(TableCollection *self, void *closure)
{
    ProvenanceTable *provenances = NULL;

    provenances = PyObject_New(ProvenanceTable, &ProvenanceTableType);
    if (provenances == NULL) {
        goto out;
    }
    provenances->table = self->tables->provenances;
    provenances->locked = false;
    provenances->tables = self;
    Py_INCREF(self);
out:
    return (PyObject *) provenances;
}

static PyObject *
TableCollection_get_sequence_length(TableCollection *self, void *closure)
{
    return Py_BuildValue("f", self->tables->sequence_length);
}

static PyObject *
TableCollection_get_file_uuid(TableCollection *self, void *closure)
{
    return Py_BuildValue("s", self->tables->file_uuid);
}

static PyObject *
TableCollection_simplify(TableCollection *self, PyObject *args, PyObject *kwds)
{
    int err;
    PyObject *ret = NULL;
    PyObject *samples = NULL;
    PyArrayObject *samples_array = NULL;
    PyArrayObject *node_map_array = NULL;
    npy_intp *shape, dims;
    size_t num_samples;
    int flags = 0;
    int filter_sites = true;
    int filter_individuals = false;
    int filter_populations = false;
    int reduce_to_site_topology = false;
    static char *kwlist[] = {
        "samples", "filter_sites", "filter_populations", "filter_individuals",
        "reduce_to_site_topology", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|iiii", kwlist,
            &samples, &filter_sites, &filter_populations, &filter_individuals,
            &reduce_to_site_topology)) {
        goto out;
    }
    samples_array = (PyArrayObject *) PyArray_FROMANY(samples, NPY_INT32, 1, 1,
            NPY_ARRAY_IN_ARRAY);
    if (samples_array == NULL) {
        goto out;
    }
    shape = PyArray_DIMS(samples_array);
    num_samples = shape[0];
    if (filter_sites) {
        flags |= TSK_FILTER_SITES;
    }
    if (filter_individuals) {
        flags |= TSK_FILTER_INDIVIDUALS;
    }
    if (filter_populations) {
        flags |= TSK_FILTER_POPULATIONS;
    }
    if (reduce_to_site_topology) {
        flags |= TSK_REDUCE_TO_SITE_TOPOLOGY;
    }

    /* Allocate a new array to hold the node map. */
    dims = self->tables->nodes->num_rows;
    node_map_array = (PyArrayObject *) PyArray_SimpleNew(1, &dims, NPY_INT32);
    if (node_map_array == NULL) {
        goto out;
    }
    err = tsk_tbl_collection_simplify(self->tables,
            PyArray_DATA(samples_array), num_samples, flags,
            PyArray_DATA(node_map_array));
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = (PyObject *) node_map_array;
    node_map_array = NULL;
out:
    Py_XDECREF(samples_array);
    Py_XDECREF(node_map_array);
    return ret;
}

static PyObject *
TableCollection_sort(TableCollection *self, PyObject *args, PyObject *kwds)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t edge_start = 0;
    int flags = 0;

    static char *kwlist[] = {"edge_start", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|n", kwlist, &edge_start)) {
        goto out;
    }
    err = tsk_tbl_collection_sort(self->tables, (size_t) edge_start, flags);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
TableCollection_compute_mutation_parents(TableCollection *self)
{
    int err;
    PyObject *ret = NULL;

    err = tsk_tbl_collection_compute_mutation_parents(self->tables, 0);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
TableCollection_deduplicate_sites(TableCollection *self)
{
    int err;
    PyObject *ret = NULL;

    err = tsk_tbl_collection_deduplicate_sites(self->tables, 0);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

/* Forward declaration */
static PyTypeObject TableCollectionType;

static PyObject *
TableCollection_equals(TableCollection *self, PyObject *args)
{
    PyObject *ret = NULL;
    TableCollection *other = NULL;

    if (!PyArg_ParseTuple(args, "O!", &TableCollectionType, &other)) {
        goto out;
    }
    ret = Py_BuildValue("i", tsk_tbl_collection_equals(self->tables, other->tables));
out:
    return ret;
}

static PyGetSetDef TableCollection_getsetters[] = {
    {"individuals", (getter) TableCollection_get_individuals, NULL, "The individual table."},
    {"nodes", (getter) TableCollection_get_nodes, NULL, "The node table."},
    {"edges", (getter) TableCollection_get_edges, NULL, "The edge table."},
    {"migrations", (getter) TableCollection_get_migrations, NULL, "The migration table."},
    {"sites", (getter) TableCollection_get_sites, NULL, "The site table."},
    {"mutations", (getter) TableCollection_get_mutations, NULL, "The mutation table."},
    {"populations", (getter) TableCollection_get_populations, NULL, "The population table."},
    {"provenances", (getter) TableCollection_get_provenances, NULL, "The provenance table."},
    {"sequence_length", (getter) TableCollection_get_sequence_length, NULL,
        "The sequence length."},
    {"file_uuid", (getter) TableCollection_get_file_uuid, NULL,
        "The UUID of the corresponding file."},
    {NULL}  /* Sentinel */
};

static PyMethodDef TableCollection_methods[] = {
    {"simplify", (PyCFunction) TableCollection_simplify, METH_VARARGS|METH_KEYWORDS,
            "Simplifies for a given sample subset." },
    {"sort", (PyCFunction) TableCollection_sort, METH_VARARGS|METH_KEYWORDS,
            "Sorts the tables to satisfy tree sequence requirements." },
    {"equals", (PyCFunction) TableCollection_equals, METH_VARARGS,
            "Returns True if the parameter table collection is equal to this one." },
    {"compute_mutation_parents", (PyCFunction) TableCollection_compute_mutation_parents,
        METH_NOARGS, "Computes the mutation parents for a the tables." },
    {"deduplicate_sites", (PyCFunction) TableCollection_deduplicate_sites,
        METH_NOARGS, "Removes sites with duplicate positions." },
    {NULL}  /* Sentinel */
};

static PyTypeObject TableCollectionType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_tskit.TableCollection",             /* tp_name */
    sizeof(TableCollection),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)TableCollection_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "TableCollection objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    TableCollection_methods,             /* tp_methods */
    0,                             /* tp_members */
    TableCollection_getsetters,           /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)TableCollection_init,      /* tp_init */
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
        tsk_treeseq_free(self->tree_sequence);
        PyMem_Free(self->tree_sequence);
        self->tree_sequence = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
TreeSequence_alloc(TreeSequence *self)
{
    int ret = -1;

    if (self->tree_sequence != NULL) {
        tsk_treeseq_free(self->tree_sequence);
        PyMem_Free(self->tree_sequence);
    }
    self->tree_sequence = PyMem_Malloc(sizeof(tsk_treeseq_t));
    if (self->tree_sequence == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->tree_sequence, 0, sizeof(*self->tree_sequence));
    ret = 0;
out:
    return ret;
}

static int
TreeSequence_init(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    self->tree_sequence = NULL;
    return 0;
}

static PyObject *
TreeSequence_dump(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    int err;
    char *path;
    PyObject *ret = NULL;
    int flags = 0;
    static char *kwlist[] = {"path", NULL};

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &path)) {
        goto out;
    }
    err = tsk_treeseq_dump(self->tree_sequence, path, flags);
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
    TableCollection *tables = NULL;
    static char *kwlist[] = {"tables", NULL};
    /* TODO add an interface to turn this on and off. */
    int flags = TSK_BUILD_INDEXES;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
            &TableCollectionType, &tables)) {
        goto out;
    }
    err = TreeSequence_alloc(self);
    if (err != 0) {
        goto out;
    }
    err = tsk_treeseq_load_tables(self->tree_sequence, tables->tables, flags);
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
    TableCollection *tables = NULL;
    static char *kwlist[] = {"tables", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
            &TableCollectionType, &tables)) {
        goto out;
    }
    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    err = tsk_treeseq_dump_tables(self->tree_sequence, tables->tables, 0);
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

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &path)) {
        goto out;
    }
    err = TreeSequence_alloc(self);
    if (err != 0) {
        goto out;
    }
    err = tsk_treeseq_load(self->tree_sequence, path, flags);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("");
out:
    return ret;
}

static PyObject *
TreeSequence_get_node(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    tsk_node_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tsk_treeseq_get_num_nodes(self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tsk_treeseq_get_node(self->tree_sequence, (size_t) record_index, &record);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_node(&record);
out:
    return ret;
}

static PyObject *
TreeSequence_get_edge(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    tsk_edge_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tsk_treeseq_get_num_edges(self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tsk_treeseq_get_edge(self->tree_sequence, (size_t) record_index, &record);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_edge(&record);
out:
    return ret;
}

static PyObject *
TreeSequence_get_migration(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    tsk_migration_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tsk_treeseq_get_num_migrations(
        self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tsk_treeseq_get_migration(self->tree_sequence,
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
TreeSequence_get_site(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    tsk_site_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tsk_treeseq_get_num_sites(self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tsk_treeseq_get_site(self->tree_sequence, (tsk_id_t) record_index, &record);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_site_object(&record);
out:
    return ret;
}

static PyObject *
TreeSequence_get_mutation(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    tsk_mutation_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tsk_treeseq_get_num_mutations(
        self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tsk_treeseq_get_mutation(self->tree_sequence,
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
TreeSequence_get_individual(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    tsk_individual_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tsk_treeseq_get_num_individuals(self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tsk_treeseq_get_individual(self->tree_sequence, (size_t) record_index, &record);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_individual_object(&record);
out:
    return ret;
}

static PyObject *
TreeSequence_get_population(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    tsk_population_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tsk_treeseq_get_num_populations(self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tsk_treeseq_get_population(self->tree_sequence, (size_t) record_index, &record);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_population(&record);
out:
    return ret;
}

static PyObject *
TreeSequence_get_provenance(TreeSequence *self, PyObject *args)
{
    int err;
    PyObject *ret = NULL;
    Py_ssize_t record_index, num_records;
    tsk_provenance_t record;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "n", &record_index)) {
        goto out;
    }
    num_records = (Py_ssize_t) tsk_treeseq_get_num_provenances(self->tree_sequence);
    if (record_index < 0 || record_index >= num_records) {
        PyErr_SetString(PyExc_IndexError, "record index out of bounds");
        goto out;
    }
    err = tsk_treeseq_get_provenance(self->tree_sequence, (size_t) record_index, &record);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = make_provenance(&record);
out:
    return ret;
}

static PyObject *
TreeSequence_get_num_edges(TreeSequence *self, PyObject *args)
{
    PyObject *ret = NULL;
    size_t num_records;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_records = tsk_treeseq_get_num_edges(self->tree_sequence);
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
    num_records = tsk_treeseq_get_num_migrations(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_records);
out:
    return ret;
}

static PyObject *
TreeSequence_get_num_individuals(TreeSequence *self, PyObject *args)
{
    PyObject *ret = NULL;
    size_t num_records;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_records = tsk_treeseq_get_num_individuals(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_records);
out:
    return ret;
}

static PyObject *
TreeSequence_get_num_populations(TreeSequence *self, PyObject *args)
{
    PyObject *ret = NULL;
    size_t num_records;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_records = tsk_treeseq_get_num_populations(self->tree_sequence);
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
    num_trees = tsk_treeseq_get_num_trees(self->tree_sequence);
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
        tsk_treeseq_get_sequence_length(self->tree_sequence));
out:
    return ret;
}

static PyObject *
TreeSequence_get_file_uuid(TreeSequence  *self)
{
    PyObject *ret = NULL;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("s", tsk_treeseq_get_file_uuid(self->tree_sequence));
out:
    return ret;
}

static PyObject *
TreeSequence_get_num_samples(TreeSequence  *self)
{
    PyObject *ret = NULL;
    size_t num_samples;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_samples = tsk_treeseq_get_num_samples(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_samples);
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
    num_nodes = tsk_treeseq_get_num_nodes(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_nodes);
out:
    return ret;
}

static PyObject *
TreeSequence_get_samples(TreeSequence *self)
{
    PyObject *ret = NULL;
    tsk_id_t *samples;
    PyObject *py_samples = NULL;
    PyObject *py_int = NULL;
    size_t j, n;
    int err;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    n = tsk_treeseq_get_num_samples(self->tree_sequence);
    err = tsk_treeseq_get_samples(self->tree_sequence, &samples);
    if (err != 0) {
        handle_library_error(err);
    }
    py_samples = PyList_New(n);
    if (py_samples == NULL) {
        goto out;
    }
    for (j = 0; j < n; j++) {
        py_int = Py_BuildValue("i", (int) samples[j]);
        if (py_int == NULL) {
            Py_DECREF(py_samples);
            goto out;
        }
        PyList_SET_ITEM(py_samples, j, py_int);
    }
    ret = py_samples;
out:
    return ret;
}

static PyObject *
TreeSequence_get_pairwise_diversity(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    PyObject *py_samples = NULL;
    static char *kwlist[] = {"samples", NULL};
    tsk_id_t *samples = NULL;
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
    err = tsk_treeseq_get_pairwise_diversity(
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
TreeSequence_genealogical_nearest_neighbours(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    static char *kwlist[] = {"focal", "reference_sets", NULL};
    tsk_id_t **reference_sets = NULL;
    size_t *reference_set_size = NULL;
    PyObject *focal = NULL;
    PyObject *reference_sets_list = NULL;
    PyArrayObject *focal_array = NULL;
    PyArrayObject **reference_set_arrays = NULL;
    PyArrayObject *ret_array = NULL;
    npy_intp *shape, dims[2];
    size_t num_focal = 0;
    size_t num_reference_sets = 0;
    size_t j;
    int err;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO!", kwlist,
            &focal, &PyList_Type, &reference_sets_list)) {
        goto out;
    }

    /* We're releasing the GIL here so we need to make sure that the memory we
     * pass to the low-level code doesn't change while it's in use. This is
     * why we take copies of the input arrays. */
    focal_array = (PyArrayObject *) PyArray_FROMANY(focal, NPY_INT32, 1, 1,
            NPY_ARRAY_IN_ARRAY|NPY_ARRAY_ENSURECOPY);
    if (focal_array == NULL) {
        goto out;
    }
    shape = PyArray_DIMS(focal_array);
    num_focal = shape[0];
    num_reference_sets = PyList_Size(reference_sets_list);
    if (num_reference_sets == 0) {
        PyErr_SetString(PyExc_ValueError, "Must have at least one sample set");
        goto out;
    }
    reference_set_size = PyMem_Malloc(num_reference_sets * sizeof(*reference_set_size));
    reference_sets = PyMem_Malloc(num_reference_sets * sizeof(*reference_sets));
    reference_set_arrays = PyMem_Malloc(num_reference_sets * sizeof(*reference_set_arrays));
    if (reference_sets == NULL || reference_set_size == NULL || reference_set_arrays == NULL) {
        goto out;
    }
    memset(reference_set_arrays, 0, num_reference_sets * sizeof(*reference_set_arrays));
    for (j = 0; j < num_reference_sets; j++) {
        reference_set_arrays[j] = (PyArrayObject *) PyArray_FROMANY(
            PyList_GetItem(reference_sets_list, j), NPY_INT32, 1, 1,
            NPY_ARRAY_IN_ARRAY|NPY_ARRAY_ENSURECOPY);
        if (reference_set_arrays[j] == NULL) {
            goto out;
        }
        reference_sets[j] = PyArray_DATA(reference_set_arrays[j]);
        shape = PyArray_DIMS(reference_set_arrays[j]);
        reference_set_size[j] = shape[0];
    }

    /* Allocate the return array */
    dims[0] = num_focal;
    dims[1] = num_reference_sets;
    ret_array = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_FLOAT64);
    if (ret_array == NULL) {
        goto out;
    }

    Py_BEGIN_ALLOW_THREADS
    err = tsk_treeseq_genealogical_nearest_neighbours(self->tree_sequence,
        PyArray_DATA(focal_array), num_focal,
        reference_sets, reference_set_size, num_reference_sets,
        0, PyArray_DATA(ret_array));
    Py_END_ALLOW_THREADS
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }

    ret = (PyObject *) ret_array;
    ret_array = NULL;
out:
    if (reference_sets != NULL) {
        PyMem_Free(reference_sets);
    }
    if (reference_set_size != NULL) {
        PyMem_Free(reference_set_size);
    }
    if (reference_set_arrays != NULL) {
        for (j = 0; j < num_reference_sets; j++) {
            Py_XDECREF(reference_set_arrays[j]);
        }
        PyMem_Free(reference_set_arrays);
    }
    Py_XDECREF(focal_array);
    Py_XDECREF(ret_array);
    return ret;
}

static PyObject *
TreeSequence_mean_descendants(TreeSequence *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    static char *kwlist[] = {"reference_sets", NULL};
    tsk_id_t **reference_sets = NULL;
    size_t *reference_set_size = NULL;
    PyObject *reference_sets_list = NULL;
    PyArrayObject **reference_set_arrays = NULL;
    PyArrayObject *ret_array = NULL;
    npy_intp *shape, dims[2];
    size_t num_reference_sets = 0;
    size_t j;
    int err;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
            &PyList_Type, &reference_sets_list)) {
        goto out;
    }

    num_reference_sets = PyList_Size(reference_sets_list);
    if (num_reference_sets == 0) {
        PyErr_SetString(PyExc_ValueError, "Must have at least one sample set");
        goto out;
    }
    reference_set_size = PyMem_Malloc(num_reference_sets * sizeof(*reference_set_size));
    reference_sets = PyMem_Malloc(num_reference_sets * sizeof(*reference_sets));
    reference_set_arrays = PyMem_Malloc(num_reference_sets * sizeof(*reference_set_arrays));
    if (reference_sets == NULL || reference_set_size == NULL || reference_set_arrays == NULL) {
        goto out;
    }
    memset(reference_set_arrays, 0, num_reference_sets * sizeof(*reference_set_arrays));
    for (j = 0; j < num_reference_sets; j++) {
        /* We're releasing the GIL here so we need to make sure that the memory we
         * pass to the low-level code doesn't change while it's in use. This is
         * why we take copies of the input arrays. */
        reference_set_arrays[j] = (PyArrayObject *) PyArray_FROMANY(
            PyList_GetItem(reference_sets_list, j), NPY_INT32, 1, 1,
            NPY_ARRAY_IN_ARRAY|NPY_ARRAY_ENSURECOPY);
        if (reference_set_arrays[j] == NULL) {
            goto out;
        }
        reference_sets[j] = PyArray_DATA(reference_set_arrays[j]);
        shape = PyArray_DIMS(reference_set_arrays[j]);
        reference_set_size[j] = shape[0];
    }

    /* Allocate the return array */
    dims[0] = tsk_treeseq_get_num_nodes(self->tree_sequence);
    dims[1] = num_reference_sets;
    ret_array = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_FLOAT64);
    if (ret_array == NULL) {
        goto out;
    }

    Py_BEGIN_ALLOW_THREADS
    err = tsk_treeseq_mean_descendants(self->tree_sequence,
        reference_sets, reference_set_size, num_reference_sets,
        0, PyArray_DATA(ret_array));
    Py_END_ALLOW_THREADS
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }

    ret = (PyObject *) ret_array;
    ret_array = NULL;
out:
    if (reference_sets != NULL) {
        PyMem_Free(reference_sets);
    }
    if (reference_set_size != NULL) {
        PyMem_Free(reference_set_size);
    }
    if (reference_set_arrays != NULL) {
        for (j = 0; j < num_reference_sets; j++) {
            Py_XDECREF(reference_set_arrays[j]);
        }
        PyMem_Free(reference_set_arrays);
    }
    Py_XDECREF(ret_array);
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
    num_mutations = tsk_treeseq_get_num_mutations(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_mutations);
out:
    return ret;
}

static PyObject *
TreeSequence_get_num_sites(TreeSequence  *self)
{
    PyObject *ret = NULL;
    size_t num_sites;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_sites = tsk_treeseq_get_num_sites(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_sites);
out:
    return ret;
}

static PyObject *
TreeSequence_get_num_provenances(TreeSequence  *self)
{
    PyObject *ret = NULL;
    size_t num_provenances;

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_provenances = tsk_treeseq_get_num_provenances(self->tree_sequence);
    ret = Py_BuildValue("n", (Py_ssize_t) num_provenances);
out:
    return ret;
}

static PyObject *
TreeSequence_get_genotype_matrix(TreeSequence  *self)
{
    PyObject *ret = NULL;
    int err;
    size_t num_sites;
    size_t num_samples;
    npy_intp dims[2];
    PyArrayObject *genotype_matrix = NULL;
    tsk_vargen_t *vg = NULL;
    char *V;
    tsk_variant_t *variant;
    size_t j;

    /* TODO add option for 16 bit genotypes */

    if (TreeSequence_check_tree_sequence(self) != 0) {
        goto out;
    }
    num_sites = tsk_treeseq_get_num_sites(self->tree_sequence);
    num_samples = tsk_treeseq_get_num_samples(self->tree_sequence);
    dims[0] = num_sites;
    dims[1] = num_samples;

    genotype_matrix = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_UINT8);
    if (genotype_matrix == NULL) {
        goto out;
    }
    V = (char *) PyArray_DATA(genotype_matrix);
    vg = PyMem_Malloc(sizeof(tsk_vargen_t));
    if (vg == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tsk_vargen_alloc(vg, self->tree_sequence, NULL, 0, 0);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    j = 0;
    while ((err = tsk_vargen_next(vg, &variant)) == 1) {
        memcpy(V + (j * num_samples), variant->genotypes.u8, num_samples * sizeof(uint8_t));
        j++;
    }
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = (PyObject *) genotype_matrix;
    genotype_matrix = NULL;
out:
    if (vg != NULL) {
        tsk_vargen_free(vg);
        PyMem_Free(vg);
    }
    Py_XDECREF(genotype_matrix);
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
    {"get_node",
        (PyCFunction) TreeSequence_get_node, METH_VARARGS,
        "Returns the node record at the specified index."},
    {"get_edge",
        (PyCFunction) TreeSequence_get_edge, METH_VARARGS,
        "Returns the edge record at the specified index."},
    {"get_migration",
        (PyCFunction) TreeSequence_get_migration, METH_VARARGS,
        "Returns the migration record at the specified index."},
    {"get_site",
        (PyCFunction) TreeSequence_get_site, METH_VARARGS,
        "Returns the mutation type record at the specified index."},
    {"get_mutation",
        (PyCFunction) TreeSequence_get_mutation, METH_VARARGS,
        "Returns the mutation record at the specified index."},
    {"get_individual",
        (PyCFunction) TreeSequence_get_individual, METH_VARARGS,
        "Returns the individual record at the specified index."},
    {"get_population",
        (PyCFunction) TreeSequence_get_population, METH_VARARGS,
        "Returns the population record at the specified index."},
    {"get_provenance",
        (PyCFunction) TreeSequence_get_provenance, METH_VARARGS,
        "Returns the provenance record at the specified index."},
    {"get_num_edges", (PyCFunction) TreeSequence_get_num_edges,
        METH_NOARGS, "Returns the number of coalescence records." },
    {"get_num_migrations", (PyCFunction) TreeSequence_get_num_migrations,
        METH_NOARGS, "Returns the number of migration records." },
    {"get_num_populations", (PyCFunction) TreeSequence_get_num_populations,
        METH_NOARGS, "Returns the number of population records." },
    {"get_num_individuals", (PyCFunction) TreeSequence_get_num_individuals,
        METH_NOARGS, "Returns the number of individual records." },
    {"get_num_trees", (PyCFunction) TreeSequence_get_num_trees,
        METH_NOARGS, "Returns the number of trees in the tree sequence." },
    {"get_sequence_length", (PyCFunction) TreeSequence_get_sequence_length,
        METH_NOARGS, "Returns the sequence length in bases." },
    {"get_file_uuid", (PyCFunction) TreeSequence_get_file_uuid,
        METH_NOARGS, "Returns the UUID of the underlying file, if present." },
    {"get_num_sites", (PyCFunction) TreeSequence_get_num_sites,
        METH_NOARGS, "Returns the number of sites" },
    {"get_num_mutations", (PyCFunction) TreeSequence_get_num_mutations, METH_NOARGS,
        "Returns the number of mutations" },
    {"get_num_provenances", (PyCFunction) TreeSequence_get_num_provenances,
        METH_NOARGS, "Returns the number of provenances" },
    {"get_num_nodes", (PyCFunction) TreeSequence_get_num_nodes, METH_NOARGS,
        "Returns the number of unique nodes in the tree sequence." },
    {"get_num_samples", (PyCFunction) TreeSequence_get_num_samples, METH_NOARGS,
        "Returns the sample size" },
    {"get_samples", (PyCFunction) TreeSequence_get_samples, METH_NOARGS,
        "Returns the samples." },
    {"get_pairwise_diversity",
        (PyCFunction) TreeSequence_get_pairwise_diversity,
        METH_VARARGS|METH_KEYWORDS, "Returns the average pairwise diversity." },
    {"genealogical_nearest_neighbours",
        (PyCFunction) TreeSequence_genealogical_nearest_neighbours,
        METH_VARARGS|METH_KEYWORDS, "Returns the genealogical nearest neighbours statistic." },
    {"mean_descendants",
        (PyCFunction) TreeSequence_mean_descendants,
        METH_VARARGS|METH_KEYWORDS, "Returns the mean number of nodes descending from each node." },
    {"get_genotype_matrix", (PyCFunction) TreeSequence_get_genotype_matrix, METH_NOARGS,
        "Returns the genotypes matrix." },
    {NULL}  /* Sentinel */
};

static PyTypeObject TreeSequenceType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_tskit.TreeSequence",             /* tp_name */
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
        tsk_tree_free(self->sparse_tree);
        PyMem_Free(self->sparse_tree);
        self->sparse_tree = NULL;
    }
    Py_XDECREF(self->tree_sequence);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

/* TODO this API should be updated to remove the SparseTreeIterator object
 * and instead support the first(), last() etc methods. Until some seeking
 * function has been called, we should be in a state that errors if any
 * methods are called.
 *
 * The _free method below is also probably redundant now and should be
 * removed.
 */
static int
SparseTree_init(SparseTree *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"tree_sequence", "flags", "tracked_samples",
        NULL};
    PyObject *py_tracked_samples = NULL;
    TreeSequence *tree_sequence = NULL;
    tsk_id_t *tracked_samples = NULL;
    int flags = 0;
    uint32_t j, num_tracked_samples, num_nodes;
    PyObject *item;

    self->sparse_tree = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|iO!", kwlist,
            &TreeSequenceType, &tree_sequence,
            &flags, &PyList_Type, &py_tracked_samples)) {
        goto out;
    }
    self->tree_sequence = tree_sequence;
    Py_INCREF(self->tree_sequence);
    if (TreeSequence_check_tree_sequence(tree_sequence) != 0) {
        goto out;
    }
    num_nodes = tsk_treeseq_get_num_nodes(tree_sequence->tree_sequence);
    num_tracked_samples = 0;
    if (py_tracked_samples != NULL) {
        if (!(flags & TSK_SAMPLE_COUNTS)) {
            PyErr_SetString(PyExc_ValueError,
                "Cannot specified tracked_samples without count_samples flag");
            goto out;
        }
        num_tracked_samples = PyList_Size(py_tracked_samples);
    }
    tracked_samples = PyMem_Malloc(num_tracked_samples * sizeof(tsk_id_t));
    if (tracked_samples == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    for (j = 0; j < num_tracked_samples; j++) {
        item = PyList_GetItem(py_tracked_samples, j);
        if (!PyNumber_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "sample must be a number");
            goto out;
        }
        tracked_samples[j] = (tsk_id_t) PyLong_AsLong(item);
        if (tracked_samples[j] >= num_nodes) {
            PyErr_SetString(PyExc_ValueError, "samples must be valid nodes");
            goto out;
        }
    }
    self->sparse_tree = PyMem_Malloc(sizeof(tsk_tree_t));
    if (self->sparse_tree == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tsk_tree_alloc(self->sparse_tree, tree_sequence->tree_sequence,
           flags);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    if (!!(flags & TSK_SAMPLE_COUNTS)) {
        err = tsk_tree_set_tracked_samples(self->sparse_tree, num_tracked_samples,
                tracked_samples);
        if (err != 0) {
            handle_library_error(err);
            goto out;
        }
    }
    ret = 0;
out:
    if (tracked_samples != NULL) {
        PyMem_Free(tracked_samples);
    }
    return ret;
}

/* TODO this should be redundant; remove */
static PyObject *
SparseTree_free(SparseTree *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    /* This method is need because we have dangling references to
     * trees after a for loop and we can't run set_sites.
     */
    tsk_tree_free(self->sparse_tree);
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
    ret = Py_BuildValue("n", (Py_ssize_t) self->sparse_tree->tree_sequence->num_samples);
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
SparseTree_get_num_roots(SparseTree *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) tsk_tree_get_num_roots(self->sparse_tree));
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
SparseTree_get_left_root(SparseTree *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("i", (int) self->sparse_tree->left_root);
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
SparseTree_is_sample(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    int node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    ret = Py_BuildValue("i", tsk_tree_is_sample(self->sparse_tree, (tsk_id_t) node));
out:
    return ret;
}

static PyObject *
SparseTree_get_parent(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    tsk_id_t parent;
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
    tsk_node_t node;
    int node_id, err;

    if (SparseTree_get_node_argument(self, args, &node_id) != 0) {
        goto out;
    }
    err = tsk_treeseq_get_node(self->sparse_tree->tree_sequence, node_id, &node);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("i", (int) node.population);
out:
    return ret;
}

static PyObject *
SparseTree_get_time(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    double time;
    int node, err;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    err = tsk_tree_get_time(self->sparse_tree, node, &time);
    if (ret != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("d", time);
out:
    return ret;
}

static PyObject *
SparseTree_get_left_child(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    tsk_id_t child;
    int node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    child = self->sparse_tree->left_child[node];
    ret = Py_BuildValue("i", (int) child);
out:
    return ret;
}

static PyObject *
SparseTree_get_right_child(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    tsk_id_t child;
    int node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    child = self->sparse_tree->right_child[node];
    ret = Py_BuildValue("i", (int) child);
out:
    return ret;
}

static PyObject *
SparseTree_get_left_sib(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    tsk_id_t sib;
    int node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    sib = self->sparse_tree->left_sib[node];
    ret = Py_BuildValue("i", (int) sib);
out:
    return ret;
}

static PyObject *
SparseTree_get_right_sib(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    tsk_id_t sib;
    int node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    sib = self->sparse_tree->right_sib[node];
    ret = Py_BuildValue("i", (int) sib);
out:
    return ret;
}

static PyObject *
SparseTree_get_children(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    int node;
    tsk_id_t u;
    size_t j, num_children;
    tsk_id_t *children = NULL;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    num_children = 0;
    for (u = self->sparse_tree->left_child[node]; u != TSK_NULL_NODE;
            u = self->sparse_tree->right_sib[u]) {
        num_children++;
    }
    children = PyMem_Malloc(num_children * sizeof(tsk_id_t));
    if (children == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    j = 0;
    for (u = self->sparse_tree->left_child[node]; u != TSK_NULL_NODE;
            u = self->sparse_tree->right_sib[u]) {
        children[j] = u;
        j++;
    }
    ret = convert_node_id_list(children, num_children);
out:
    if (children != NULL) {
        PyMem_Free(children);
    }
    return ret;
}

static bool
SparseTree_check_sample_list(SparseTree *self)
{
    bool ret = tsk_tree_has_sample_lists(self->sparse_tree);
    if (! ret) {
        PyErr_SetString(PyExc_ValueError,
            "Sample lists not supported. Please set sample_lists=True.");
    }
    return ret;
}

static PyObject *
SparseTree_get_right_sample(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    tsk_id_t sample_index;
    int node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    if (!SparseTree_check_sample_list(self)) {
        goto out;
    }
    sample_index = self->sparse_tree->right_sample[node];
    ret = Py_BuildValue("i", (int) sample_index);
out:
    return ret;
}

static PyObject *
SparseTree_get_left_sample(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    tsk_id_t sample_index;
    int node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    if (!SparseTree_check_sample_list(self)) {
        goto out;
    }
    sample_index = self->sparse_tree->left_sample[node];
    ret = Py_BuildValue("i", (int) sample_index);
out:
    return ret;
}

static PyObject *
SparseTree_get_next_sample(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    tsk_id_t out_index;
    int in_index, num_samples;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "I", &in_index)) {
        goto out;
    }
    num_samples = (int) tsk_treeseq_get_num_samples(self->sparse_tree->tree_sequence);
    if (in_index < 0 || in_index >= num_samples) {
        PyErr_SetString(PyExc_ValueError, "Sample index out of bounds");
        goto out;
    }
    if (!SparseTree_check_sample_list(self)) {
        goto out;
    }
    out_index = self->sparse_tree->next_sample[in_index];
    ret = Py_BuildValue("i", (int) out_index);
out:
    return ret;
}

static PyObject *
SparseTree_get_mrca(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    int err;
    tsk_id_t mrca;
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
    err = tsk_tree_get_mrca(self->sparse_tree, (tsk_id_t) u,
            (tsk_id_t) v, &mrca);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("i", (int) mrca);
out:
    return ret;
}

static PyObject *
SparseTree_get_num_samples(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    size_t num_samples;
    int err, node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    err = tsk_tree_get_num_samples(self->sparse_tree, (tsk_id_t) node,
            &num_samples);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("I", (unsigned int) num_samples);
out:
    return ret;
}

static PyObject *
SparseTree_get_num_tracked_samples(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;
    size_t num_tracked_samples;
    int err, node;

    if (SparseTree_get_node_argument(self, args, &node) != 0) {
        goto out;
    }
    err = tsk_tree_get_num_tracked_samples(self->sparse_tree, (tsk_id_t) node,
            &num_tracked_samples);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = Py_BuildValue("I", (unsigned int) num_tracked_samples);
out:
    return ret;
}

static PyObject *
SparseTree_get_sites(SparseTree *self, PyObject *args)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = convert_sites(self->sparse_tree->sites, self->sparse_tree->sites_length);
out:
    return ret;
}

static PyObject *
SparseTree_get_num_sites(SparseTree  *self)
{
    PyObject *ret = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("n", (Py_ssize_t) self->sparse_tree->sites_length);
out:
    return ret;
}

static PyObject *
SparseTree_get_newick(SparseTree *self, PyObject *args, PyObject *kwds)
{
    PyObject *ret = NULL;
    static char *kwlist[] = {"root", "precision", NULL};
    int precision = 14;
    int root, err;
    size_t buffer_size;
    char *buffer = NULL;

    if (SparseTree_check_sparse_tree(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|i", kwlist, &root, &precision)) {
        goto out;
    }
    if (precision < 0 || precision > 16) {
        PyErr_SetString(PyExc_ValueError, "Precision must be between 0 and 16, inclusive");
        goto out;
    }
    buffer_size = tsk_treeseq_get_num_nodes(self->sparse_tree->tree_sequence);
    /* For every node, we have roughly precision bytes, plus bracketing and leading values.
     * This is a rough guess, so add 10 just to be on the safe side. We might need
     * to be more precise with this though if we have large time values.
     */
    buffer_size *= precision + 10;
    buffer = PyMem_Malloc(buffer_size);
    if (buffer == NULL) {
        PyErr_NoMemory();
    }
    err = tsk_convert_newick(self->sparse_tree, (tsk_id_t) root, precision, 0,
            buffer_size, buffer);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = PyBytes_FromString(buffer);
out:
    if (buffer != NULL) {
        PyMem_Free(buffer);
    }
    return ret;
}

static PyMemberDef SparseTree_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef SparseTree_methods[] = {
    {"free", (PyCFunction) SparseTree_free, METH_NOARGS,
            "Frees the underlying tree object." },
    {"get_sample_size", (PyCFunction) SparseTree_get_sample_size, METH_NOARGS,
            "Returns the number of samples in this tree." },
    {"get_num_nodes", (PyCFunction) SparseTree_get_num_nodes, METH_NOARGS,
            "Returns the number of nodes in the sparse tree." },
    {"get_num_roots", (PyCFunction) SparseTree_get_num_roots, METH_NOARGS,
            "Returns the number of roots in the sparse tree." },
    {"get_index", (PyCFunction) SparseTree_get_index, METH_NOARGS,
            "Returns the index this tree occupies within the tree sequence." },
    {"get_left_root", (PyCFunction) SparseTree_get_left_root, METH_NOARGS,
            "Returns the root of the tree." },
    {"get_left", (PyCFunction) SparseTree_get_left, METH_NOARGS,
            "Returns the left-most coordinate (inclusive)." },
    {"get_right", (PyCFunction) SparseTree_get_right, METH_NOARGS,
            "Returns the right-most coordinate (exclusive)." },
    {"get_sites", (PyCFunction) SparseTree_get_sites, METH_NOARGS,
            "Returns the list of sites on this tree." },
    {"get_flags", (PyCFunction) SparseTree_get_flags, METH_NOARGS,
            "Returns the value of the flags variable." },
    {"get_num_sites", (PyCFunction) SparseTree_get_num_sites, METH_NOARGS,
            "Returns the number of sites on this tree." },
    {"is_sample", (PyCFunction) SparseTree_is_sample, METH_VARARGS,
            "Returns True if the specified node is a sample." },
    {"get_parent", (PyCFunction) SparseTree_get_parent, METH_VARARGS,
            "Returns the parent of node u" },
    {"get_time", (PyCFunction) SparseTree_get_time, METH_VARARGS,
            "Returns the time of node u" },
    {"get_population", (PyCFunction) SparseTree_get_population, METH_VARARGS,
            "Returns the population of node u" },
    {"get_left_child", (PyCFunction) SparseTree_get_left_child, METH_VARARGS,
            "Returns the left-most child of node u" },
    {"get_right_child", (PyCFunction) SparseTree_get_right_child, METH_VARARGS,
            "Returns the right-most child of node u" },
    {"get_left_sib", (PyCFunction) SparseTree_get_left_sib, METH_VARARGS,
            "Returns the left-most sib of node u" },
    {"get_right_sib", (PyCFunction) SparseTree_get_right_sib, METH_VARARGS,
            "Returns the right-most sib of node u" },
    {"get_children", (PyCFunction) SparseTree_get_children, METH_VARARGS,
            "Returns the children of u in left-right order." },
    {"get_left_sample", (PyCFunction) SparseTree_get_left_sample, METH_VARARGS,
            "Returns the index of the left-most sample descending from u." },
    {"get_right_sample", (PyCFunction) SparseTree_get_right_sample, METH_VARARGS,
            "Returns the index of the right-most sample descending from u." },
    {"get_next_sample", (PyCFunction) SparseTree_get_next_sample, METH_VARARGS,
            "Returns the index of the next sample after the specified sample index." },
    {"get_mrca", (PyCFunction) SparseTree_get_mrca, METH_VARARGS,
            "Returns the MRCA of nodes u and v" },
    {"get_num_samples", (PyCFunction) SparseTree_get_num_samples, METH_VARARGS,
            "Returns the number of samples below node u." },
    {"get_num_tracked_samples", (PyCFunction) SparseTree_get_num_tracked_samples,
            METH_VARARGS,
            "Returns the number of tracked samples below node u." },
    {"get_newick", (PyCFunction) SparseTree_get_newick,
            METH_VARARGS|METH_KEYWORDS,
            "Returns the newick representation of this tree." },
    {NULL}  /* Sentinel */
};

static PyTypeObject SparseTreeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_tskit.SparseTree",             /* tp_name */
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
        tsk_diff_iter_free(self->tree_diff_iterator);
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
    self->tree_diff_iterator = PyMem_Malloc(sizeof(tsk_diff_iter_t));
    if (self->tree_diff_iterator == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->tree_diff_iterator, 0, sizeof(tsk_diff_iter_t));
    err = tsk_diff_iter_alloc(self->tree_diff_iterator,
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
    double left, right;
    size_t list_size, j;
    tsk_edge_list_t *records_out, *records_in, *record;

    if (TreeDiffIterator_check_state(self) != 0) {
        goto out;
    }
    err = tsk_diff_iter_next(self->tree_diff_iterator, &left, &right,
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
            value = Py_BuildValue("ddii", record->edge.left, record->edge.right,
                    record->edge.parent, record->edge.child);
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
            value = Py_BuildValue("ddii", record->edge.left, record->edge.right,
                    record->edge.parent, record->edge.child);
            if (value == NULL) {
                goto out;
            }
            PyList_SET_ITEM(in_list, j, value);
            record = record->next;
            j++;
        }
        ret = Py_BuildValue("(dd)OO", left, right, out_list, in_list);
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
    "_tskit.TreeDiffIterator",             /* tp_name */
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
        err = tsk_tree_first(self->sparse_tree->sparse_tree);
        self->first = 0;
    } else {
        err = tsk_tree_next(self->sparse_tree->sparse_tree);
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
    "_tskit.SparseTreeIterator",             /* tp_name */
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
 * VcfConverter
 *===================================================================
 */

static int
VcfConverter_check_state(VcfConverter *self)
{
    int ret = 0;
    if (self->tsk_vcf_converter == NULL) {
        PyErr_SetString(PyExc_SystemError, "converter not initialised");
        ret = -1;
    }
    return ret;
}

static void
VcfConverter_dealloc(VcfConverter* self)
{
    if (self->tsk_vcf_converter != NULL) {
        tsk_vcf_converter_free(self->tsk_vcf_converter);
        PyMem_Free(self->tsk_vcf_converter);
        self->tsk_vcf_converter = NULL;
    }
    Py_XDECREF(self->tree_sequence);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
VcfConverter_init(VcfConverter *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"tree_sequence", "ploidy", "contig_id", NULL};
    unsigned int ploidy = 1;
    const char *contig_id = "1";
    TreeSequence *tree_sequence;

    self->tsk_vcf_converter = NULL;
    self->tree_sequence = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|Is", kwlist,
            &TreeSequenceType, &tree_sequence, &ploidy, &contig_id)) {
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
    if (strlen(contig_id) == 0) {
        PyErr_SetString(PyExc_ValueError, "contig_id cannot be the empty string");
        goto out;
    }
    self->tsk_vcf_converter = PyMem_Malloc(sizeof(tsk_vcf_converter_t));
    if (self->tsk_vcf_converter == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = tsk_vcf_converter_alloc(self->tsk_vcf_converter,
            self->tree_sequence->tree_sequence, ploidy, contig_id);
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
    err = tsk_vcf_converter_next(self->tsk_vcf_converter, &record);
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
    err = tsk_vcf_converter_get_header(self->tsk_vcf_converter, &header);
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
    "_tskit.VcfConverter",             /* tp_name */
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
        tsk_hapgen_free(self->haplotype_generator);
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
    self->haplotype_generator = PyMem_Malloc(sizeof(tsk_hapgen_t));
    if (self->haplotype_generator == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->haplotype_generator, 0, sizeof(tsk_hapgen_t));
    err = tsk_hapgen_alloc(self->haplotype_generator,
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
    err = tsk_hapgen_get_haplotype(self->haplotype_generator,
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
    "_tskit.HaplotypeGenerator",             /* tp_name */
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
        tsk_vargen_free(self->variant_generator);
        PyMem_Free(self->variant_generator);
        self->variant_generator = NULL;
    }
    Py_XDECREF(self->tree_sequence);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
VariantGenerator_init(VariantGenerator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    int err;
    static char *kwlist[] = {"tree_sequence", "samples", NULL};
    TreeSequence *tree_sequence = NULL;
    PyObject *samples_input = Py_None;
    PyArrayObject *samples_array = NULL;
    tsk_id_t *samples = NULL;
    size_t num_samples = 0;
    npy_intp *shape;

    /* TODO add option for 16 bit genotypes */
    self->variant_generator = NULL;
    self->tree_sequence = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|O", kwlist,
            &TreeSequenceType, &tree_sequence, &samples_input)) {
        goto out;
    }
    self->tree_sequence = tree_sequence;
    Py_INCREF(self->tree_sequence);
    if (TreeSequence_check_tree_sequence(self->tree_sequence) != 0) {
        goto out;
    }
    if (samples_input != Py_None) {
        samples_array = (PyArrayObject *) PyArray_FROMANY(samples_input, NPY_INT32, 1, 1,
                NPY_ARRAY_IN_ARRAY);
        if (samples_array == NULL) {
            goto out;
        }
        shape = PyArray_DIMS(samples_array);
        num_samples = (size_t) shape[0];
        samples = PyArray_DATA(samples_array);
    }
    self->variant_generator = PyMem_Malloc(sizeof(tsk_vargen_t));
    if (self->variant_generator == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    /* Note: the vargen currently takes a copy of the samples list. If we wanted
     * to avoid this we would INCREF the samples array above and keep a reference
     * to in the object struct */
    err = tsk_vargen_alloc(self->variant_generator,
            self->tree_sequence->tree_sequence, samples, num_samples, 0);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    ret = 0;
out:
    Py_XDECREF(samples_array);
    return ret;
}

static PyObject *
VariantGenerator_next(VariantGenerator *self)
{
    PyObject *ret = NULL;
    tsk_variant_t *var;
    int err;

    if (VariantGenerator_check_state(self) != 0) {
        goto out;
    }
    err = tsk_vargen_next(self->variant_generator, &var);
    if (err < 0) {
        handle_library_error(err);
        goto out;
    }
    if (err == 1) {
        ret = make_variant(var, self->variant_generator->num_samples);
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
    "_tskit.VariantGenerator",             /* tp_name */
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
        tsk_ld_calc_free(self->ld_calc);
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
    self->ld_calc = PyMem_Malloc(sizeof(tsk_ld_calc_t));
    if (self->ld_calc == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    memset(self->ld_calc, 0, sizeof(tsk_ld_calc_t));
    err = tsk_ld_calc_alloc(self->ld_calc, self->tree_sequence->tree_sequence);
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
    err = tsk_ld_calc_get_r2(self->ld_calc, (size_t) a, (size_t) b, &r2);
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
    int direction = TSK_DIR_FORWARD;
    size_t num_r2_values = 0;
    int buffer_acquired = 0;

    if (LdCalculator_check_state(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "On|ind", kwlist,
            &dest, &source_index, &direction, &max_mutations, &max_distance)) {
        goto out;
    }
    if (direction != TSK_DIR_FORWARD && direction != TSK_DIR_REVERSE) {
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
    err = tsk_ld_calc_get_r2_array(
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
    "_tskit.LdCalculator",             /* tp_name */
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

static PyMethodDef tskit_methods[] = {
    {NULL}        /* Sentinel */
};

/* Initialisation code supports Python 2.x and 3.x. The framework uses the
 * recommended structure from http://docs.python.org/howto/cporting.html.
 * I've ignored the point about storing state in globals, as the examples
 * from the Python documentation still use this idiom.
 */

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef tskitmodule = {
    PyModuleDef_HEAD_INIT,
    "_tskit",   /* name of module */
    MODULE_DOC, /* module documentation, may be NULL */
    -1,
    tskit_methods,
    NULL, NULL, NULL, NULL
};

#define INITERROR return NULL

PyObject *
PyInit__tskit(void)

#else
#define INITERROR return

void
init_tskit(void)
#endif
{
    kas_version_t kastore_version;
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&tskitmodule);
#else
    PyObject *module = Py_InitModule3("_tskit", tskit_methods, MODULE_DOC);
#endif
    if (module == NULL) {
        INITERROR;
    }
    import_array();

    /* Set up the dynamic API for kastore. */
    kas_dynamic_api = (kas_funcptr *) PyCapsule_Import("_kastore._C_API", 0);
    if (kas_dynamic_api == NULL) {
        INITERROR;
    }
    kastore_version = kas_version();
    if (kastore_version.major != KAS_VERSION_MAJOR) {
        PyErr_SetString(PyExc_RuntimeError, "kastore C API major version mismatch");
        INITERROR;
    }
    /* If the minor version of kastore we compiled against is older than
     * the one we've just loaded it's OK. But it's not safe to run when we've
     * compiled against a newer version.
     */
    if (kastore_version.minor > KAS_VERSION_MINOR) {
        PyErr_SetString(PyExc_RuntimeError, "kastore C API minor version mismatch");
        INITERROR;
    }

    /* LightweightTableCollection type */
    LightweightTableCollectionType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&LightweightTableCollectionType) < 0) {
        INITERROR;
    }
    Py_INCREF(&LightweightTableCollectionType);
    PyModule_AddObject(module, "LightweightTableCollection",
            (PyObject *) &LightweightTableCollectionType);

    /* IndividualTable type */
    IndividualTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&IndividualTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&IndividualTableType);
    PyModule_AddObject(module, "IndividualTable", (PyObject *) &IndividualTableType);

    /* NodeTable type */
    NodeTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&NodeTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&NodeTableType);
    PyModule_AddObject(module, "NodeTable", (PyObject *) &NodeTableType);

    /* EdgeTable type */
    EdgeTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&EdgeTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&EdgeTableType);
    PyModule_AddObject(module, "EdgeTable", (PyObject *) &EdgeTableType);

    /* MigrationTable type */
    MigrationTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&MigrationTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&MigrationTableType);
    PyModule_AddObject(module, "MigrationTable", (PyObject *) &MigrationTableType);

    /* SiteTable type */
    SiteTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&SiteTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&SiteTableType);
    PyModule_AddObject(module, "SiteTable", (PyObject *) &SiteTableType);

    /* MutationTable type */
    MutationTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&MutationTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&MutationTableType);
    PyModule_AddObject(module, "MutationTable", (PyObject *) &MutationTableType);

    /* PopulationTable type */
    PopulationTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PopulationTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&PopulationTableType);
    PyModule_AddObject(module, "PopulationTable", (PyObject *) &PopulationTableType);

    /* ProvenanceTable type */
    ProvenanceTableType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&ProvenanceTableType) < 0) {
        INITERROR;
    }
    Py_INCREF(&ProvenanceTableType);
    PyModule_AddObject(module, "ProvenanceTable", (PyObject *) &ProvenanceTableType);

    /* TableCollectionTable type */
    TableCollectionType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&TableCollectionType) < 0) {
        INITERROR;
    }
    Py_INCREF(&TableCollectionType);
    PyModule_AddObject(module, "TableCollection", (PyObject *) &TableCollectionType);

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
    TskitLibraryError = PyErr_NewException("_tskit.LibraryError", NULL, NULL);
    Py_INCREF(TskitLibraryError);
    PyModule_AddObject(module, "LibraryError", TskitLibraryError);
    TskitFileFormatError = PyErr_NewException("_tskit.FileFormatError", NULL, NULL);
    Py_INCREF(TskitFileFormatError);
    PyModule_AddObject(module, "FileFormatError", TskitFileFormatError);
    TskitVersionTooNewError = PyErr_NewException("_tskit.VersionTooNewError", NULL, NULL);
    Py_INCREF(TskitVersionTooNewError);
    PyModule_AddObject(module, "VersionTooNewError", TskitVersionTooNewError);
    TskitVersionTooOldError = PyErr_NewException("_tskit.VersionTooOldError", NULL, NULL);
    Py_INCREF(TskitVersionTooOldError);
    PyModule_AddObject(module, "VersionTooOldError", TskitVersionTooOldError);

    /* Node flags */
    PyModule_AddIntConstant(module, "NODE_IS_SAMPLE", TSK_NODE_IS_SAMPLE);

    /* Tree flags */
    PyModule_AddIntConstant(module, "SAMPLE_COUNTS", TSK_SAMPLE_COUNTS);
    PyModule_AddIntConstant(module, "SAMPLE_LISTS", TSK_SAMPLE_LISTS);
    /* Directions */
    PyModule_AddIntConstant(module, "FORWARD", TSK_DIR_FORWARD);
    PyModule_AddIntConstant(module, "REVERSE", TSK_DIR_REVERSE);

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
