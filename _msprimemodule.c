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

#include "lib/msprime.h"

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


#if 0
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



static int
msprime_parse_event_classes(PyObject *py_events, event_class_t *events)
{
    int ret = -1;
    int j, size;
    double rate, u, r;
    PyObject *item, *value;
    size = PyList_Size(py_events);
    if (size == 0) {
        PyErr_SetString(MsprimeInputError, "must have > 0 events");
        goto out;
    }
    for (j = 0; j < size; j++) {
        item = PyList_GetItem(py_events, j);
        if (!PyDict_Check(item)) {
            PyErr_SetString(MsprimeInputError, "not a dictionary");
            goto out;
        }
        value = get_dict_number(item, "rate");
        if (value == NULL) {
            goto out;
        }
        rate = PyFloat_AsDouble(value);
        value = get_dict_number(item, "r");
        if (value == NULL) {
            goto out;
        }
        r = PyFloat_AsDouble(value);
        value = get_dict_number(item, "u");
        if (value == NULL) {
            goto out;
        }
        u = PyFloat_AsDouble(value);
        events[j].rate = rate;
        events[j].r = r;
        events[j].u = u;
    }
    ret = 0;
out:
    return ret;
}
#endif

/*===================================================================
 * Simulator
 *===================================================================
 */
#if 0
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
Simulator_parse_events(Simulator *self, PyObject *py_events)
{
    int ret = -1;
    int size;
    size = PyList_Size(py_events);
    if (size == 0) {
        PyErr_SetString(MsprimeInputError, "must have > 0 events");
        goto out;
    }
    self->sim->num_event_classes = size;
    self->sim->event_classes = PyMem_Malloc(size * sizeof(event_class_t));
    if (self->sim->event_classes == NULL) {
        ret = ERR_ALLOC_FAILED;
        goto out;
    }
    ret = msprime_parse_event_classes(py_events, self->sim->event_classes);
out:
    return ret;
}
#endif

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

#if 0
static int
Simulator_check_input(Simulator *self)
{
    int ret = -1;
    unsigned int j;
    sim_t *sim = self->sim;
    event_class_t *e;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (sim->dimension < 1 || sim->dimension > 2) {
        handle_input_error("dimension must be 1 or 2");
        goto out;
    }
    if (sim->dimension == 1 && sim->pixel_size != 2.0) {
        handle_input_error("pixel size must be 2.0 in 1D");
        goto out;
    }
    if (sim->simulate_pedigree < 0 || sim->simulate_pedigree > 1) {
        handle_input_error("simulate_pedigree must be 0 or 1");
        goto out;
    }
    if (sim->simulate_pedigree == 1 && sim->num_loci != 1) {
        handle_input_error("m must be 1 for pedigree simulation");
        goto out;
    }
    if (sim->torus_diameter <= 0.0) {
        handle_input_error("must have torus_edge > 0");
        goto out;
    }
    if (sim->num_loci == 0) {
        handle_input_error("must have num_loci > 0");
        goto out;
    }
    if (sim->num_parents == 0) {
        handle_input_error("must have num_parents > 0");
        goto out;
    }
    if (sim->max_population_size == 0) {
        handle_input_error("must have max_population_size > 0");
        goto out;
    }
    if (sim->max_occupancy == 0) {
        handle_input_error("must have max_occupancy > 0");
        goto out;
    }
    if (sim->recombination_probability < 0 ||
            sim->recombination_probability > 1) {
        handle_input_error("must have 0 <= recombination_probability <= 1");
        goto out;
    }
    if (sim->pixel_size <= 0 || sim->pixel_size > sim->torus_diameter / 4) {
        handle_input_error("must have 0 < pixel_size <= L/4 ");
        goto out;
    }
    if (fmod(sim->torus_diameter, sim->pixel_size) != 0.0) {
        handle_input_error("L/s must be an integer");
        goto out;
    }
    if (sim->num_event_classes == 0) {
        handle_input_error("at least one event class required");
        goto out;
    }
    for (j = 0; j < sim->num_event_classes; j++) {
        e = &sim->event_classes[j];
        if (e->r <= 0.0 || e->r > sim->torus_diameter / 4.0) {
            handle_input_error("must have 0 < r < L / 4");
            goto out;
        }
        if (e->u <= 0.0 || e->u >= 1.0) {
            handle_input_error("must have 0 < u < 1");
            goto out;
        }
        if (e->rate <= 0.0) {
            handle_input_error("must have 0 < rate < 1");
            goto out;
        }
    }
    ret = 0;
out:
    return ret;
}
#endif

static int
Simulator_init(Simulator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
#if 0
    int sim_ret;
    static char *kwlist[] = {"sample", "event_classes", "num_loci",
            "num_parents", "max_population_size", "max_occupancy",
            "dimension", "simulate_pedigree", "random_seed", "torus_diameter",
            "pixel_size", "recombination_probability", NULL};
    PyObject *sample, *events;
    sim_t *sim = PyMem_Malloc(sizeof(sim_t));
    self->sim = sim;
    if (self->sim == NULL) {
        goto out;
    }
    memset(self->sim, 0, sizeof(sim_t));
    sim->num_loci = 1;
    sim->num_parents = 2;
    sim->torus_diameter = 1000;
    sim->pixel_size = 2;
    sim->recombination_probability = 0.5;
    sim->random_seed = 1;
    sim->max_population_size = 1000;
    sim->max_occupancy = 10;
    sim->dimension = 2;
    sim->simulate_pedigree = 0;
    sim->max_time = DBL_MAX;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|IIIIIIkddd", kwlist,
            &PyList_Type, &sample,
            &PyList_Type, &events,
            &sim->num_loci, &sim->num_parents, &sim->max_population_size,
            &sim->max_occupancy, &sim->dimension, &sim->simulate_pedigree,
            &sim->random_seed, &sim->torus_diameter, &sim->pixel_size,
            &sim->recombination_probability)) {
        goto out;
    }
    if (Simulator_parse_sample(self, sample) != 0) {
        goto out;
    }
    if (Simulator_parse_events(self, events) != 0) {
        goto out;
    }
    if (Simulator_check_input(self) != 0) {
        goto out;
    }
    sim_ret = sim_alloc(self->sim);
    if (sim_ret != 0) {
        handle_library_error(sim_ret);
        goto out;
    }
    sim_ret = sim_initialise(self->sim);
    if (sim_ret != 0) {
        handle_library_error(sim_ret);
        goto out;
    }
    ret = 0;
out:
#endif
    return ret;
}

static PyMemberDef Simulator_members[] = {
    {NULL}  /* Sentinel */
};


#if 0

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
Simulator_get_num_parents(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->sim->num_parents);
out:
    return ret;
}

static PyObject *
Simulator_get_dimension(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->sim->dimension);
out:
    return ret;
}

static PyObject *
Simulator_get_simulate_pedigree(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->sim->simulate_pedigree);
out:
    return ret;
}

static PyObject *
Simulator_get_max_population_size(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->sim->max_population_size);
out:
    return ret;
}

static PyObject *
Simulator_get_max_occupancy(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("I", self->sim->max_occupancy);
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
Simulator_get_torus_diameter(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sim->torus_diameter);
out:
    return ret;
}

static PyObject *
Simulator_get_pixel_size(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sim->pixel_size);
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
Simulator_get_num_reproduction_events(Simulator  *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("K",
            (unsigned long long) self->sim->num_reproduction_events);
out:
    return ret;
}

static PyObject *
Simulator_get_recombination_probability(Simulator *self)
{
    PyObject *ret = NULL;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    ret = Py_BuildValue("d", self->sim->recombination_probability);
out:
    return ret;
}

static PyObject *
Simulator_get_event_classes(Simulator *self)
{
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *d = NULL;
    unsigned int j;
    event_class_t *e;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    l = PyList_New(self->sim->num_event_classes);
    if (l == NULL) {
        goto out;
    }
    for (j = 0; j < self->sim->num_event_classes; j++) {
        e = &self->sim->event_classes[j];
        d = Py_BuildValue("{s:d,s:d,s:d}", "r", e->r, "u", e->u,
                "rate", e->rate);
        if (d == NULL) {
            goto out;
        }
        if (PyList_SetItem(l, j, d) != 0) {
            goto out;
        }
    }
    ret = l;
    l = NULL;
out:
    Py_XDECREF(l);
    return ret;
}



static PyObject *
Simulator_individual_to_python(Simulator *self, individual_t *ind)
{
    PyObject *ret = NULL;
    PyObject *key, *value;
    int status;
    double *x = ind->location;
    avl_node_t *node;
    int_map_value_t *imv;
    PyObject *ancestry = NULL;
    PyObject *loc = NULL;
    if (self->sim->dimension == 1) {
        loc = Py_BuildValue("d", x[0]);
    } else {
        loc = Py_BuildValue("(d,d)", x[0], x[1]);
    }
    if (loc == NULL) {
        goto out;
    }
    if (self->sim->simulate_pedigree == 1) {
        ret = loc;
    } else {
        ancestry = PyDict_New();
        if (ancestry == NULL) {
            goto out;
        }
        for (node = ind->ancestry.head; node != NULL; node = node->next) {
            imv = (int_map_value_t *) node->item;
            key = Py_BuildValue("I", imv->key);
            if (key == NULL) {
                goto out;
            }
            value = Py_BuildValue("I", imv->value);
            if (value == NULL) {
                Py_DECREF(key);
                goto out;
            }
            status = PyDict_SetItem(ancestry, key, value);
            Py_DECREF(key);
            Py_DECREF(value);
            if (status != 0) {
                goto out;
            }
        }
        ret = PyTuple_Pack(2, loc, ancestry);
        if (ret == NULL) {
            goto out;
        }
    }
out:
    if (self->sim->simulate_pedigree == 0) {
        Py_XDECREF(loc);
        Py_XDECREF(ancestry);
    }
    return ret;
}

static PyObject *
Simulator_get_population(Simulator  *self)
{
    int err;
    unsigned int j;
    PyObject *ret = NULL;
    PyObject *l = NULL;
    PyObject *py_ind = NULL;
    avl_tree_t *pop = NULL;
    avl_node_t *node;
    uint64_t id;
    uintptr_t int_ptr;
    individual_t *ind;

    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    pop = PyMem_Malloc(sizeof(avl_tree_t));
    if (pop == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    err = sim_get_population(self->sim, pop);
    if (err != 0) {
        handle_library_error(err);
        goto out;
    }
    l = PyList_New(avl_count(pop));
    if (l == NULL) {
        goto out;
    }
    j = 0;
    for (node = pop->head; node != NULL; node = node->next) {
        id = *((uint64_t *) node->item);
        int_ptr = (uintptr_t) id;
        ind = (individual_t *) int_ptr;
        py_ind = Simulator_individual_to_python(self, ind);
        if (py_ind == NULL) {
            goto out;
        }
        if (PyList_SetItem(l, j, py_ind) != 0) {
            Py_DECREF(py_ind);
            goto out;
        }
        j++;
    }
    ret = l;
    l = NULL;
out:
    if (pop != NULL) {
        sim_free_population(self->sim, pop);
        PyMem_Free(pop);
    }
    Py_XDECREF(l);
    return ret;
}


static PyObject *
Simulator_get_history(Simulator  *self)
{
    PyObject *ret = NULL;
    PyObject *pi = NULL;
    PyObject *tau = NULL;
    PyObject *pi_locus, *tau_locus;
    unsigned int j, l, n;
    int err;
    sim_t *sim = self->sim;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (self->sim->simulate_pedigree == 1) {
        PyErr_SetString(PyExc_NotImplementedError,
                "Cannot get history for pedigree simulation");
        goto out;
    }

    pi = PyList_New(sim->num_loci);
    if (pi == NULL) {
        goto out;
    }
    tau = PyList_New(sim->num_loci);
    if (tau == NULL) {
        goto out;
    }
    n = 2 * sim->sample_size;
    for (l = 0; l < sim->num_loci; l++) {
        pi_locus = PyList_New(n);
        if (pi_locus == NULL) {
            goto out;
        }
        err = PyList_SetItem(pi, l, pi_locus);
        if (err < 0) {
            goto out;
        }
        tau_locus = PyList_New(n);
        if (tau_locus == NULL) {
            goto out;
        }
        err = PyList_SetItem(tau, l, tau_locus);
        if (err < 0) {
            goto out;
        }
        for (j = 0; j < n; j++) {
            err = PyList_SetItem(pi_locus, j, PyLong_FromLong(sim->pi[l][j]));
            if (err < 0) {
                goto out;
            }
            err = PyList_SetItem(tau_locus, j,
                    PyFloat_FromDouble(sim->tau[l][j]));
            if (err < 0) {
                goto out;
            }
        }
    }
    ret = Py_BuildValue("(O, O)", pi, tau);
out:
    Py_XDECREF(pi);
    Py_XDECREF(tau);

    return ret;
}

static PyObject *
Simulator_run(Simulator *self, PyObject *args)
{
    PyObject *ret = NULL;
    int status, not_done;
    uint64_t chunk = 8192;
    double max_time = DBL_MAX;
    if (Simulator_check_sim(self) != 0) {
        goto out;
    }
    if (!PyArg_ParseTuple(args, "|d", &max_time)) {
        goto out;
    }
    not_done = 1;
    self->sim->max_time = max_time;
    while (not_done) {
        status = sim_simulate(self->sim, chunk);
        if (status < 0) {
            handle_library_error(status);
            goto out;
        }
        not_done = status != 0;
        if (PyErr_CheckSignals() < 0) {
            goto out;
        }
    }
    /* return True if complete coalescence has occured */
    ret = self->sim->time < max_time ? Py_True : Py_False;
    Py_INCREF(ret);
out:
    return ret;
}

#endif
static PyMethodDef Simulator_methods[] = {
#if 0
    {"get_num_loci", (PyCFunction) Simulator_get_num_loci, METH_NOARGS,
            "Returns the number of loci" },
    {"get_num_parents", (PyCFunction) Simulator_get_num_parents, METH_NOARGS,
            "Returns the number of parents" },
    {"get_dimension", (PyCFunction) Simulator_get_dimension, METH_NOARGS,
            "Returns the dimension of the simulation." },
    {"get_simulate_pedigree", (PyCFunction) Simulator_get_simulate_pedigree, METH_NOARGS,
            "Returns 1 if we are simulating the pedigree; 0 otherwise." },
    {"get_max_population_size", (PyCFunction) Simulator_get_max_population_size,
            METH_NOARGS,
            "Returns the maximum size of the ancestral population"},
    {"get_max_occupancy", (PyCFunction) Simulator_get_max_occupancy,
            METH_NOARGS,
            "Returns the maximum occupancy of a single pixel"},
    {"get_random_seed", (PyCFunction) Simulator_get_random_seed, METH_NOARGS,
            "Returns the random seed" },
    {"get_torus_diameter", (PyCFunction) Simulator_get_torus_diameter,
            METH_NOARGS, "Returns the torus diameter" },
    {"get_pixel_size", (PyCFunction) Simulator_get_pixel_size, METH_NOARGS,
            "Returns the size of a pixel" },
    {"get_time", (PyCFunction) Simulator_get_time, METH_NOARGS,
            "Returns the current simulation time" },
    {"get_num_reproduction_events",
            (PyCFunction) Simulator_get_num_reproduction_events, METH_NOARGS,
            "Returns the number of reproduction events up to this point." },
    {"get_recombination_probability",
            (PyCFunction) Simulator_get_recombination_probability, METH_NOARGS,
            "Returns the probability of recombination between adjacent loci" },
    {"get_event_classes", (PyCFunction) Simulator_get_event_classes, METH_NOARGS,
            "Returns the event classes" },
    {"get_population", (PyCFunction) Simulator_get_population, METH_NOARGS,
            "Returns the state of the ancestral population" },
    {"get_history", (PyCFunction) Simulator_get_history, METH_NOARGS,
            "Returns the history of the sample as a tuple (pi, tau)" },
    {"run", (PyCFunction) Simulator_run, METH_VARARGS,
            "Simulates until at most the specified time. Returns True\
            if the required stopping conditions have been met and False \
            otherwise." },
#endif
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
    MsprimeInputError = PyErr_NewException("_msprime.InputError", NULL,
            NULL);
    Py_INCREF(MsprimeInputError);
    PyModule_AddObject(module, "InputError", MsprimeInputError);
    MsprimeLibraryError = PyErr_NewException("_msprime.LibraryError",
            NULL, NULL);
    Py_INCREF(MsprimeLibraryError);
    PyModule_AddObject(module, "LibraryError", MsprimeLibraryError);

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}


