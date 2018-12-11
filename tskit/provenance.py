"""
Common provenance methods used to determine the state and versions
of various dependencies and the OS.
"""
from __future__ import print_function
from __future__ import division

import platform
import json
import os.path

import jsonschema

import tskit.exceptions as exceptions

__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass


_environment = None


def _get_environment():
    env = {
        "libraries": {
            "kastore": {
                # Hard coding this for now as there is no way to get version
                # information from the kastore C API. See
                # https://github.com/tskit-dev/kastore/issues/41
                # We could import the kastore module here and use its version,
                # but this is not the same as the C code we have compiler against.
                "version": "0.1.0",
            }
        },
        "os": {
            "system": platform.system(),
            "node": platform.node(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
        },
        "python": {
            "implementation": platform.python_implementation(),
            "version": platform.python_version(),
        }
    }
    return env


def get_environment():
    """
    Returns a dictionary describing the environment in which tskit
    is currently running.
    """
    # Everything here is fixed so we cache it
    global _environment
    if _environment is None:
        _environment = _get_environment()
    return _environment


def get_provenance_dict(parameters=None):
    """
    Returns a dictionary encoding an execution of tskit conforming to the
    provenance schema.
    """
    document = {
        "schema_version": "1.0.0",
        "software": {
            "name": "tskit",
            "version": __version__,
        },
        "parameters": parameters,
        "environment": get_environment()
    }
    return document


def validate_provenance(provenance):
    """
    Validates the specified dict-like object against the tskit
    :ref:`provenance schema <sec_provenance>`. If the input does
    not represent a valid instance of the schema an exception is
    raised.

    :param dict provenance: The dictionary representing a JSON document
        to be validated against the schema.
    :raises: :class:`.ProvenanceValidationError`
    """
    base = os.path.dirname(__file__)
    schema_file = os.path.join(base, "provenance.schema.json")
    with open(schema_file) as f:
        schema = json.load(f)
    try:
        jsonschema.validate(provenance, schema)
    except jsonschema.exceptions.ValidationError as ve:
        raise exceptions.ProvenanceValidationError(str(ve))
