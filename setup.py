import os
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
            warn(f"Error occurred getting GSL path config: {e}")
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
        # Try vcpkg on Windows first
        if IS_WINDOWS:
            vcpkg_root = os.getenv("VCPKG_ROOT")
            if vcpkg_root is None:
                vcpkg_root = os.getenv("VCPKG_INSTALLATION_ROOT")
            if vcpkg_root and os.path.exists(vcpkg_root):
                gsl_include = os.path.join(
                    vcpkg_root, "installed", "x64-windows", "include"
                )
                gsl_lib = os.path.join(vcpkg_root, "installed", "x64-windows", "lib")
                if os.path.exists(gsl_include) and os.path.exists(gsl_lib):
                    self.include_dirs.append(gsl_include)
                    self.library_dirs.append(gsl_lib)
                    return

        # Fallback to gsl-config on Unix-like systems
        output = self._run_command(["gsl-config", "--cflags"]).split()
        if len(output) > 0:
            token = output[0]
            self.include_dirs.append(token[2:])
        output = self._run_command(["gsl-config", "--libs"]).split()
        for token in output:
            if token.startswith("-L"):
                self.library_dirs.append(token[2:])


def get_extension():
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
    tsk_source_files = ["core.c", "tables.c", "trees.c", "genotypes.c"]
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

    return Extension(
        "msprime._msprime",
        sources=sources,
        extra_compile_args=["-std=c99"],
        libraries=libraries,
        define_macros=defines,
        include_dirs=includes + configurator.include_dirs + [numpy.get_include()],
        library_dirs=configurator.library_dirs,
    )


def main():
    setup(
        ext_modules=[get_extension()],
    )


if __name__ == "__main__":
    main()
