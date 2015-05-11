import re
import sys
try:
    from setuptools import setup, Extension
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup, Extension

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

f = open("README.txt")
msprime_readme = f.read()
f.close()
msprime_version = parse_version("msprime/__init__.py")

requirements = ["pkgconfig"]

# TODO proper pkgconfig setup. We need to try and see if GSL and HDF5
# are installed and then used pkgconfig to get the paths. If pkg-config
# isn't available, we make a guess.
import pkgconfig
pkg_info = pkgconfig.parse('gsl hdf5')
# On older systems, HDF5 is not supported by pkg-config. This however, we
# might still be able to compile and run just by making sure we link against
# hdf5. This works for Debian Wheezy.
if 'hdf5' not in pkg_info["libraries"]:
    pkg_info["libraries"].add('hdf5')

d = "lib/"
_msprime_module = Extension('_msprime',
    sources=[
        "_msprimemodule.c", d + "msprime.c", d + "fenwick.c", d + "avl.c",
        d + "tree_sequence.c", d + "object_heap.c"],
    # Enable asserts by default.
    undef_macros=['NDEBUG'],
    libraries=list(pkg_info["libraries"]),
    include_dirs = [d] + list(pkg_info["include_dirs"]),
    library_dirs = list(pkg_info["library_dirs"]),
)

setup(
    name="msprime",
    version=msprime_version,
    long_description=msprime_readme,
    packages=["msprime"],
    author="Jerome Kelleher",
    author_email="jerome.kelleher@well.ox.ac.uk",
    url="http://pypi.python.org/pypi/msprime",
    entry_points={
        'console_scripts': [
            'mspms=msprime.cli:msp_ms_main',
        ]
    },
    install_requires=requirements,
    ext_modules = [_msprime_module],
    keywords = ["Coalescent simulation", "ms"],
    license = "GNU LGPLv3+",
    platforms = ["POSIX"],
    classifiers = [
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
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: POSIX",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
