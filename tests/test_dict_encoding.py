"""
Test cases for the low-level dictionary encoding used to move
data around in C.
"""
import lwt_interface.dict_encoding_testlib

from msprime import _msprime

lwt_interface.dict_encoding_testlib.lwt_module = _msprime
# Bring the tests defined in dict_encoding_testlib into the current namespace
# so pytest will find and execute them.
from lwt_interface.dict_encoding_testlib import *  # noqa
