#
# Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

import re
import subprocess

# First, we try to use setuptools. If it's not available locally,
# we fall back on ez_setup.
try:
    from setuptools import setup, Extension
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup, Extension


class PathConfigurator(object):
    """
    A class to attempt configuration of the compile search paths
    on various platforms.
    """
    def __init__(self):
        # TODO: make some other guesses for this...
        self.include_dirs = []
        self.library_dirs = []
        self._check_hdf5_version()
        self._attempt_pkgconfig()

    def _check_hdf5_version(self):
        try:
            output = subprocess.check_output(["h5ls", "-V"]).split()
            version_str = output[2]
            version = list(map(int, version_str.split(b".")[:2]))
            if version < [1, 8]:
                # TODO is there a better exception to raise here?
                raise ValueError(
                    "hdf5 version {} found; we need 1.8.0 or greater".format(
                        version_str))
        except OSError as e:
            if e.errno == 2:
                print("Cannot find h5ls: is HDF5 installed?:")
            else:
                print("Error occured running h5ls:", e)

    def _run_pkgconfig(self, cmd):
        pkgconfig = "pkg-config"
        packages = ["gsl", "hdf5"]
        cmd = [pkgconfig] + cmd + packages
        output = subprocess.check_output(cmd).split()
        # Strip off the leading -I or -L
        return [arg[2:].decode() for arg in output]

    def _get_pkgconfig_list(self, option):
        ret = []
        try:
            ret = self._run_pkgconfig([option])
        except OSError as e:
            print("pkg-config error (not installed?):", e)
        except subprocess.CalledProcessError as e:
            print("pkg-config failed:", e)
        return ret

    def _attempt_pkgconfig(self):
        self.library_dirs = self._get_pkgconfig_list("--libs-only-L")
        self.include_dirs = self._get_pkgconfig_list("--cflags-only-I")


# Following the recommendations of PEP 396 we parse the version number
# out of the module.
def parse_version(module_file):
    """
    Parses the version string from the specified file.

    This implementation is ugly, but there doesn't seem to be a good way
    to do this in general at the moment.
    """
    f = open(module_file)
    s = f.read()
    f.close()
    match = re.findall("__version__ = '([^']+)'", s)
    return match[0]

# Now, setup the extension module.
msprime_version = parse_version("msprime/__init__.py")
configurator = PathConfigurator()
d = "lib/"
_msprime_module = Extension(
    '_msprime',
    sources=[
        "_msprimemodule.c", d + "msprime.c", d + "fenwick.c", d + "avl.c",
        d + "tree_sequence.c", d + "object_heap.c", d + "newick.c",
        d + "hapgen.c"],
    # Enable asserts by default.
    undef_macros=["NDEBUG"],
    # We define this macro to ensure we're using the v18 versions of the
    # HDF5 API and not earlier deprecated versions.
    define_macros=[
        ("H5_NO_DEPRECATED_SYMBOLS", None),
        ("MSP_LIBRARY_VERSION_STR", '"{}"'.format(msprime_version)),
    ],
    libraries=["gsl", "gslcblas", "hdf5"],
    include_dirs=[d] + configurator.include_dirs,
    library_dirs=configurator.library_dirs,
)

f = open("README.txt")
msprime_readme = f.read()
f.close()

setup(
    name="msprime",
    version=msprime_version,
    description="A fast and accurate coalescent simulator.",
    long_description=msprime_readme,
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
    install_requires=["svgwrite"],
    ext_modules=[_msprime_module],
    keywords=["Coalescent simulation", "ms"],
    license="GNU LGPLv3+",
    platforms=["POSIX"],
    classifiers=[
        "Programming Language :: C",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.1",
        "Programming Language :: Python :: 3.2",
        "Programming Language :: Python :: 3.3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        (
            "License :: OSI Approved :: GNU Lesser General Public License "
            "v3 or later (LGPLv3+)"
        ),
        "Operating System :: POSIX",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
