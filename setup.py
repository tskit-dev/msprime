import os.path
import platform
import subprocess
from warnings import warn

import numpy
from setuptools import Extension
from setuptools import setup


CONDA_PREFIX = os.getenv("MSP_CONDA_PREFIX", None)
IS_WINDOWS = platform.system() == "Windows"


class PathConfigurator:
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
            warn(f"Error occured getting GSL path config: {e}")
        # If the conda prefix is defined, then we are compiling in a conda
        # context. All include and lib paths should come from within this prefix.
        if CONDA_PREFIX is not None:
            prefix = CONDA_PREFIX
            if IS_WINDOWS:
                prefix = os.path.join(prefix, "Library")
            self.library_dirs.append(os.path.join(prefix, "lib"))
            self.include_dirs.append(os.path.join(prefix, "include"))

    def _run_command(self, args):
        return subprocess.check_output(args, text=True)

    def _configure_gsl(self):
        output = self._run_command(["gsl-config", "--cflags"]).split()
        if len(output) > 0:
            token = output[0]
            self.include_dirs.append(token[2:])
        output = self._run_command(["gsl-config", "--libs"]).split()
        for token in output:
            if token.startswith("-L"):
                self.library_dirs.append(token[2:])


libdir = "lib"
tskroot = os.path.join("git-submodules", "tskit", "c")
tskdir = os.path.join(tskroot, "tskit")
kasdir = os.path.join(tskroot, "subprojects", "kastore")
includes = ["lwt_interface", libdir, tskroot, tskdir, kasdir]

configurator = PathConfigurator()
msp_source_files = [
    "msprime.c",
    "fenwick.c",
    "avl.c",
    "util.c",
    "object_heap.c",
    "rate_map.c",
    "mutgen.c",
    "likelihood.c",
]
tsk_source_files = ["core.c", "tables.c", "trees.c"]
kas_source_files = ["kastore.c"]

sources = (
    ["msprime/_msprimemodule.c"]
    + [os.path.join(libdir, f) for f in msp_source_files]
    + [os.path.join(tskdir, f) for f in tsk_source_files]
    + [os.path.join(kasdir, f) for f in kas_source_files]
)

libraries = ["gsl", "gslcblas"]
defines = []
if IS_WINDOWS:
    # Needed for generating UUIDs
    libraries.append("Advapi32")
    defines += [
        # These two are required for GSL to compile and link against the
        # conda-forge version.
        ("GSL_DLL", None),
        ("WIN32", None),
    ]

_msprime_module = Extension(
    "msprime._msprime",
    sources=sources,
    extra_compile_args=["-std=c99"],
    libraries=libraries,
    define_macros=defines,
    include_dirs=includes + configurator.include_dirs + [numpy.get_include()],
    library_dirs=configurator.library_dirs,
)


setup(
    # The package name along with all the other metadata is specified in setup.cfg
    # However, GitHub's dependency graph can't see the package unless we put this here.
    name="msprime",
    ext_modules=[_msprime_module],
)
