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

requirements = []

try:
    import pkgconfig
    pkg_info = pkgconfig.parse('gsl hdf5')
except ImportError:
    pkg_info = {"include_dirs":[], "library_dirs":[]}

d = "lib/"
_msprime_module = Extension('_msprime',
    sources=[
        "_msprimemodule.c", d + "msprime.c", d + "fenwick.c", d + "avl.c",
        d + "tree_sequence.c", d + "object_heap.c", d + "newick.c",
        d + "hapgen.c"],
    # Enable asserts by default.
    undef_macros=["NDEBUG"],
    # We define this macro to ensure we're using the v18 versions of the
    # HDF5 API and not earlier deprecated versions.
    define_macros=[("H5_NO_DEPRECATED_SYMBOLS", None)],
    libraries=["gsl", "gslcblas", "hdf5"],
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
    install_requires=[],
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
