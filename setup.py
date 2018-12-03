#
# Copyright (C) 2015-2017 University of Oxford
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

import platform
import os
import os.path
import glob

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext


IS_WINDOWS = platform.system() == "Windows"


# Obscure magic required to allow numpy be used as an 'setup_requires'.
class build_ext(_build_ext):
    def finalize_options(self):
        super(build_ext, self).finalize_options()
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())


# The above obscure magic doesn't seem to work on py2 and prevents the
# extension from building at all, so here's a nasty workaround:
libdir = "lib"
includes = [libdir]
try:
    import numpy
    includes.append(numpy.get_include())
except ImportError:
    pass

kastore_dir = os.path.join("kastore", "c")
gsl_dir = "gsl"

source_files = glob.glob("lib/gsl/rng/*.c")
source_files += glob.glob("lib/gsl/err/*.c")
source_files += glob.glob("lib/gsl/randist/*.c")
source_files += glob.glob("lib/gsl/sys/*.c")
source_files += [os.path.join("lib", "gsl", "specfunc", f) for f in [
    "gamma.c", "trig.c", "psi.c", "log.c", "zeta.c", "elementary.c", "exp.c"]]
source_files += [os.path.join("lib", "gsl", "complex", f) for f in [
    "math.c"]]

source_files += [os.path.join("lib", f) for f in [
    "msprime.c", "fenwick.c", "avl.c", "tree_sequence.c",
    "object_heap.c", "newick.c", "hapgen.c", "recomb_map.c", "mutgen.c",
    "vargen.c", "vcf.c", "ld.c", "tables.c", "util.c", "uuid.c",
    os.path.join(kastore_dir, "kastore.c")]]


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
            ("MSP_LIBRARY_VERSION_STR", '{}'.format(self._msprime_version)),
        ]
        return defines[index]


libraries = []
if IS_WINDOWS:
    # Needed for generating UUIDs
    libraries.append("Advapi32")

_msprime_module = Extension(
    '_msprime',
    sources=["_msprimemodule.c"] + source_files,
    # Enable asserts by default.
    undef_macros=["NDEBUG"],
    extra_compile_args=["-std=c99"],
    libraries=libraries,
    define_macros=DefineMacros(),
    include_dirs=includes + [
        os.path.join(libdir, kastore_dir),
        os.path.join(libdir, gsl_dir)]
)

with open("README.rst") as f:
    long_description = f.read()

setup(
    name="msprime",
    description="A fast and accurate coalescent simulator.",
    long_description=long_description,
    packages=["msprime"],
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
    install_requires=["numpy>=1.7.0", "h5py", "svgwrite", "six", "jsonschema"],
    ext_modules=[_msprime_module],
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
    setup_requires=['numpy', 'setuptools_scm'],
    use_scm_version={"write_to": "msprime/_version.py"},
)
