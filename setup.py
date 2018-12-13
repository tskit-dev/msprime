#
# Copyright (C) 2015-2018 University of Oxford
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
from __future__ import division
from __future__ import print_function

import subprocess
import platform
import os
import os.path
import sys
from warnings import warn

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


CONDA_PREFIX = os.getenv("MSP_CONDA_PREFIX", None)
IS_WINDOWS = platform.system() == "Windows"


class PathConfigurator(object):
    """
    A class to attempt configuration of the compile search paths
    on various platforms.
    """
    def __init__(self):
        self.include_dirs = []
        self.library_dirs = []
        try:
            self._configure_gsl()
        except OSError as e:
            warn("Error occured getting GSL path config: {}".format(e))
        # If the conda prefix is defined, then we are compiling in a conda
        # context. All include and lib paths should come from within this prefix.
        if CONDA_PREFIX is not None:
            prefix = CONDA_PREFIX
            if IS_WINDOWS:
                prefix = os.path.join(prefix, "Library")
            self.library_dirs.append(os.path.join(prefix, "lib"))
            self.include_dirs.append(os.path.join(prefix, "include"))

    def _run_command(self, args):
        return subprocess.check_output(args, universal_newlines=True)

    def _configure_gsl(self):
        output = self._run_command(["gsl-config", "--cflags"]).split()
        if len(output) > 0:
            token = output[0]
            self.include_dirs.append(token[2:])
        output = self._run_command(["gsl-config", "--libs"]).split()
        for token in output:
            if token.startswith("-L"):
                self.library_dirs.append(token[2:])


# Obscure magic required to allow numpy be used as an 'setup_requires'.
# Based on https://stackoverflow.com/questions/19919905
class local_build_ext(build_ext):
    def finalize_options(self):
        build_ext.finalize_options(self)
        if sys.version_info[0] >= 3:
            import builtins
        else:
            import __builtin__ as builtins
        # Prevent numpy from thinking it is still in its setup process:
        builtins.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())
        import kastore
        self.include_dirs.append(kastore.get_include())


libdir = "lib"
includes = [libdir, libdir + "/tskit"]

configurator = PathConfigurator()
msp_source_files = [
    "msprime.c", "fenwick.c", "avl.c", "util.c",
    "object_heap.c", "recomb_map.c", "mutgen.c"
]
tsk_source_files = [
    # TODO this will be removed once we move the tskit code out.
    "tskit/tsk_core.c",
    "tskit/tsk_tables.c",
    "tskit/tsk_trees.c",
    "tskit/tsk_genotypes.c",
    "tskit/tsk_stats.c",
    "tskit/tsk_convert.c",
]


# Now, setup the extension module. We have to do some quirky workarounds
# here so that we can get the current version number from setuptools_scm
# and also get this version provided as a compile time parameter to the
# extension module.
class DefineMacros(object):
    def __init__(self):
        self._msprime_version = None

    def __getitem__(self, index):
        if self._msprime_version is None:
            import setuptools_scm
            version = setuptools_scm.get_version()
            if IS_WINDOWS:
                self._msprime_version = '\\"{}\\"'.format(version)
            else:
                self._msprime_version = '"{}"'.format(version)

        defines = [
            # Define the library version
            # TODO: this is only used for the VCF converter to get the right version
            # in the header. We'll need something smarter in the future.
            ("TSK_LIBRARY_VERSION_STR", '{}'.format(self._msprime_version)),
            # Keeping this for now for compiling the C module.
            ("MSP_LIBRARY_VERSION_STR", '{}'.format(self._msprime_version)),
        ]
        if IS_WINDOWS:
            defines += [
                # These two are required for GSL to compile and link against the
                # conda-forge version.
                ("GSL_DLL", None), ("WIN32", None)]
        defines += [("KAS_DYNAMIC_API", None)]
        return defines[index]


libraries = ["gsl", "gslcblas"]
if IS_WINDOWS:
    # Needed for generating UUIDs
    libraries.append("Advapi32")

_msprime_module = Extension(
    '_msprime',
    sources=["_msprimemodule.c"] + [
        os.path.join(libdir, f) for f in msp_source_files + tsk_source_files],
    # Enable asserts by default.
    undef_macros=["NDEBUG"],
    extra_compile_args=["-std=c99"],
    libraries=libraries,
    define_macros=DefineMacros(),
    include_dirs=includes + configurator.include_dirs,
    library_dirs=configurator.library_dirs,
)

_tskit_module = Extension(
    '_tskit',
    sources=["_tskitmodule.c"] + [os.path.join(libdir, f) for f in tsk_source_files],
    # Enable asserts by default.
    undef_macros=["NDEBUG"],
    extra_compile_args=["-std=c99"],
    libraries=libraries,
    define_macros=DefineMacros(),
    include_dirs=includes + configurator.include_dirs,
    library_dirs=configurator.library_dirs,
)

numpy_ver = "numpy>=1.7"
kastore_ver = "kastore>=0.2.2"

with open("README.rst") as f:
    long_description = f.read()

setup(
    name="msprime",
    description="A fast and accurate coalescent simulator.",
    long_description=long_description,
    packages=["msprime", "tskit"],
    author="Jerome Kelleher",
    author_email="jerome.kelleher@well.ox.ac.uk",
    url="http://pypi.python.org/pypi/msprime",
    entry_points={
        'console_scripts': [
            'mspms=msprime.cli:mspms_main',
            'msp=msprime.cli:msp_main',
        ]
    },
    include_package_data=True,
    install_requires=[numpy_ver, kastore_ver, "h5py", "svgwrite", "six", "jsonschema"],
    ext_modules=[_msprime_module, _tskit_module],
    keywords=["Coalescent simulation", "ms"],
    license="GNU GPLv3+",
    platforms=["POSIX", "Windows", "MacOS X"],
    classifiers=[
        "Programming Language :: C",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: POSIX",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    setup_requires=[numpy_ver, kastore_ver, 'setuptools_scm'],
    use_scm_version={"write_to": "msprime/_version.py"},
    cmdclass={"build_ext": local_build_ext},
)
