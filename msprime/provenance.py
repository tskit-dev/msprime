#
# Copyright (C) 2016-2017 University of Oxford
#
# This file is part of msprime.
#
# msprime is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# msprime is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with msprime.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Common provenance methods used to determine the state and versions
of various dependencies and the OS.
"""
import base64
import importlib
import json
import marshal
import logging
import types

import numpy
import tskit

import _msprime

__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass
logger = logging.getLogger(__name__)


def get_provenance_dict(parameters=None):
    """
    Returns a dictionary encoding an execution of msprime conforming to the
    tskit provenance schema.
    """
    document = {
        "schema_version": "1.0.0",
        "software": {
            "name": "msprime",
            "version": __version__,
        },
        "parameters": parameters,
        "environment": get_environment()
    }
    return document


def _get_environment():
    gsl_version = ".".join(map(str, _msprime.get_gsl_version()))
    libraries = {"gsl": {"version": gsl_version}}
    return tskit.provenance.get_environment(extra_libs=libraries)


_environment = None


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


class ProvenanceEncoderDecoder(json.JSONEncoder):
    """
    Extension of the `json` encoder that serializes arbitrary python class objects
    by calling an `asdict` method on the object that should provide the arguments
    needed to exactly recreate it. Special cases numpy arrays, classes and functions
    """

    def default(self, obj):
        if isinstance(obj, types.FunctionType):
            # Some provenance-relevant classes such as SimulationModelChange take
            # functional arguments. Note that we have a tight definition of
            # what we allow as a function here. We don't include builtins or
            # callable classes, we also block functions that refer to globals
            # or the outer scope for both complexity and security reasons.
            # We can deserialize these later with
            # FunctionType(fcode, fglobals, fname, fdefaults, fclosure)
            if len(obj.__code__.co_names) > 0:
                error = (f'Configuration function {obj.__code__.co_name} refers to'
                         f' global variables {obj.__code__.co_names} which are not'
                         f' currently serialized')
                logger.warning(error)
                return {'__error__': error}
            if obj.__closure__ is not None and len(obj.__closure__) > 0:
                error = (f'Configuration function {obj.__code__.co_name} refers to outer'
                         f' scope variables which are not currently serialized')
                logger.warning(error)
                return {'__error__': error}
            return {
                '__function__': base64.b64encode(
                    marshal.dumps(obj.__code__)).decode("utf-8"),
                'defaults': obj.__defaults__
            }

        if isinstance(obj, tskit.TreeSequence):
            # This is a special case in which we know that the input tree sequence
            # is strictly a subset of the current ts. Therefore, we can reproduce
            # the input by storing the offsets within tables.

            # Note: depends on current semantics of __iter__ in TableCollection.
            # See https://github.com/tskit-dev/tskit/issues/500
            offsets = {name: len(table) for name, table in obj.tables}
            return {'__from_ts__': offsets}

        elif isinstance(obj, numpy.ndarray):
            # The most failsafe way would be be to `base64.b64encode(obj.tostring())`
            # but as we expect only simple arrays we make this trade-off so that
            # the array is human readable
            return {
                '__ndarray__': obj.tolist(),
                'dtype': obj.dtype.str
            }

        elif isinstance(obj, numpy.number):
            # Some numbers come through as numpy types
            return {
                '__npgeneric__': str(obj),
                'dtype': obj.dtype.str,
            }

        try:
            ret = obj.asdict()
            cls = obj.__class__
            ret['__class__'] = f'{cls.__module__}.{cls.__name__}'
            return ret
        except AttributeError:
            raise TypeError(f'Object of type {obj.__class__.__name__} '
                            f'is not JSON serializable, provide an `asdict` method '
                            f'that gives the objects args')

    @staticmethod
    def decode(s, ts):
        def hook(obj):
            if '__function__' in obj:
                return types.FunctionType(
                    marshal.loads(base64.b64decode(obj['__function__'])),
                    {}, 'func', obj['defaults'], ())
            elif '__from_ts__' in obj:
                return subset_tree_sequence(ts, obj['__from_ts__'])
            elif '__ndarray__' in obj:
                return numpy.asarray(obj['__ndarray__'], dtype=obj['dtype'])
            elif '__npgeneric__' in obj:
                return numpy.array([obj['__npgeneric__']]).astype(obj['dtype'])[0]
            elif '__class__' in obj:
                module, cls = obj['__class__'].rsplit('.', 1)
                module = importlib.import_module(module)
                del obj['__class__']
                return getattr(module, cls)(**obj)
            return obj
        return json.JSONDecoder(object_hook=hook).decode(s)


def subset_tree_sequence(ts, offsets):
    tables = ts.dump_tables()
    # Note: depends on current semantics of __iter__ in TableCollection.
    # See https://github.com/tskit-dev/tskit/issues/500
    table_map = {name: table for name, table in tables}
    for name in ["individuals", "nodes", "edges", "provenances"]:
        table_map[name].truncate(offsets[name])
    for name in ["sites", "mutations", "migrations"]:
        if len(table_map[name]) != offsets[name]:
            raise ValueError(f"Table {name} in from_ts has been modified")
    return tables.tree_sequence()


def parse_provenance(ts, index=None):
    """
    Convert the specified provenance to tuple of the command used and dict suitable
    for calling the command again e.g `msprime.simulate(**parse_provenance(ts)[1])`
    """
    if index is None:
        index = ts.num_provenances - 1
    prov = ts.provenance(index).record
    ret = ProvenanceEncoderDecoder.decode(prov, ts)
    if ret['software']['name'] != 'msprime':
        raise ValueError(f'Only msprime provanances can be parsed,'
                         f' found {ret["software"]["name"]}')
    command = ret['parameters']['command']
    del ret['parameters']['command']
    return command, ret['parameters']


def _human_readable_size(size, decimal_places=2):
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            break
        size /= 1024.0
    return f"{size:.{decimal_places}f}{unit}"


def json_encode_provenance(provenance_dict, num_replicates=1):
    """
        Return a JSON representation of the provenance
    """
    prov = ProvenanceEncoderDecoder().encode(provenance_dict)
    if len(prov) > 2_097_152:
        logger.warning(f'The provenance information for the resulting tree sequence is'
                       f' {_human_readable_size(len(prov))}.'
                       f' This is nothing to worry about as provenance is a good thing'
                       f' to have, but if you want to save this memory/storage space'
                       f' you can disable provenance recording by setting'
                       f' record_provenance=False')
    return prov
