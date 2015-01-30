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
v = sys.version_info[:2]
if v < (2, 7) or v == (3, 0) or v == (3, 1):
    requirements.append("argparse")

d = "lib/"
_msprime_module = Extension('_msprime',
    sources = ["_msprimemodule.c", d + "msprime.c", d + "fenwick.c",
            d + "tree_file.c", d + "hapgen.c", d + "newick.c",
            d + "avl.c"],
    libraries = ["gsl", "gslcblas"],
    # Enable asserts by default.
    undef_macros=['NDEBUG'],
    include_dirs = [d])

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
            'msp-ms=msprime.cli:msp_ms_main',
        ]
    },
    install_requires=requirements,
    ext_modules = [_msprime_module],
)
