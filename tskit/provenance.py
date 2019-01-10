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
import _tskit

__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass


# NOTE: the APIs here are all preliminary. We should have a class that encapsulates
# all of the required functionality, including parsing and printing out provenance
# records. This will replace the current functions.


def get_environment(extra_libs=None, include_tskit=True):
    """
    Returns a dictionary describing the environment in which tskit
    is currently running.

    This API is tentative and will change in the future when a more
    comprehensive provenance API is implemented.
    """
    env = {
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
    libs = {
        "kastore": {
            "version": ".".join(map(str, _tskit.get_kastore_version()))
        }
    }
    if include_tskit:
        libs["tskit"] = {"version": __version__}
    if extra_libs is not None:
        libs.update(extra_libs)
    env["libraries"] = libs
    return env


def get_provenance_dict(parameters=None):
    """
    Returns a dictionary encoding an execution of tskit conforming to the
    provenance schema.
    """
    document = {
        "schema_version": "1.0.0",
        "software": {
            "name": "tskit",
            "version": __version__
        },
        "parameters": parameters,
        "environment": get_environment(include_tskit=False)
    }
    return document


# Cache the schema
_schema = None


def get_schema():
    """
    Returns the tskit provenance :ref:`provenance schema <sec_provenance>` as
    a dict.

    :return: The provenance schema.
    :rtype: dict
    """
    global _schema
    if _schema is None:
        base = os.path.dirname(__file__)
        schema_file = os.path.join(base, "provenance.schema.json")
        with open(schema_file) as f:
            _schema = json.load(f)
    # Return a copy to avoid issues with modifying the cached schema
    return dict(_schema)


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
    schema = get_schema()
    try:
        jsonschema.validate(provenance, schema)
    except jsonschema.exceptions.ValidationError as ve:
        raise exceptions.ProvenanceValidationError(str(ve))
