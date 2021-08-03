(sec_installation)=

# Installation

There are three options for installing `msprime`:

1. {ref}`sec_installation_conda`: the recommended options for most users
2. {ref}`sec_installation_pip`: more flexibility and control for advanced users
3. {ref}`sec_installation_container`: to use
   `msprime` inside an extremely isolated environment

(sec_installation_conda)=

## Via Conda

Convenient [conda](<https://conda.io/docs/>) packages are available for Python
3.7+ on Linux, OSX and Windows.
These pre-built binary packages are built using
[conda-forge](<https://conda-forge.org/>).

### Quick Start

1. Install `conda` using [miniconda](<https://conda.io/miniconda.html>).
   Make sure you follow the instructions to fully activate your `conda`
   installation!
2. Set up the [conda-forge channel](<https://conda-forge.org/>) using
   `conda config --add channels conda-forge`.
3. Install msprime: `conda install msprime`.
4. Try it out: `msp --version`.

There are several different ways to obtain `conda`. Please see the
[anaconda installation documentation](<https://docs.anaconda.com/anaconda/install/>)
for full details.

### Full example

In this example we create a
[conda environment](<https://conda.io/docs/user-guide/tasks/manage-environments.html>)
and install `msprime` into it on an OSX machine.
We assume that `conda` has been installed  and bash shell is being used (Windows users will need to modify the
commands slightly).

```{code-block} none

$ conda create --name msprime-env
...

  environment location: /home/jean/miniconda3/envs/msprime-env

Proceed ([y]/n)? y

Preparing transaction: ...working... done
Verifying transaction: ...working... done
Executing transaction: ...working... done
#
# To activate this environment, use
#
#     $ conda activate msprime-env
#
# To deactivate an active environment, use
#
#     $ conda deactivate

$ conda activate msprime-env
(msprime-env) $ conda install -c conda-forge msprime
...

  added / updated specs:
    - msprime

The following NEW packages will be INSTALLED:
...

Proceed ([y]/n)? y

Downloading and Extracting Packages

Preparing transaction: ...working... done
Verifying transaction: ...working... done
Executing transaction: ...working... done

(msprime-env) $ python
Python 3.8.5
[GCC 7.3.0] :: Anaconda, Inc. on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import msprime
>>> ts = msprime.sim_ancestry(3)
>>> print(ts.draw_text())
3.18┊      10     ┊  
    ┊    ┏━━┻━━┓  ┊  
2.08┊    9     ┃  ┊  
    ┊  ┏━┻━┓   ┃  ┊  
0.80┊  ┃   ┃   8  ┊  
    ┊  ┃   ┃  ┏┻┓ ┊  
0.42┊  ┃   7  ┃ ┃ ┊  
    ┊  ┃  ┏┻┓ ┃ ┃ ┊  
0.34┊  6  ┃ ┃ ┃ ┃ ┊  
    ┊ ┏┻┓ ┃ ┃ ┃ ┃ ┊  
0.00┊ 0 4 1 3 2 5 ┊  
  0.00          1.00 
```

Please see the [conda documentation](<https://conda.io/docs/index.html>) for
more details on managing packages and environments.

(sec_installation_pip)=

## Via Pip

Installing using `pip` is more flexible than `conda` as it
can support more versions of Python and dependencies can be customised.

:::{warning}

If you installed Python using anaconda/miniconda, please install
msprime using `conda`. Link/load time errors will occur if you mix
a conda installed Python with system libraries.

:::

To install msprime via pip, first
{ref}`sec_installation_system_requirements`, then {ref}`sec_run_pip`.

(sec_installation_system_requirements)=

### check system requirements

Msprime requires Python 3.7+ and pip. Run

```{code-block} bash

$ python3 -m pip --version

```

to make sure you see Python 3.7 or greater.
If you do not, do installation {ref}`via conda <sec_installation_conda>`.

On most platforms, pip will install `msprime` pre-built binaries
without any additional requirements.
Follow {ref}`sec_pip_install_source` to install from sources instead pre-built
binaries.

(sec_run_pip)=

### run pip install

We can install `msprime` easily using pip:

```
$ python3 -m pip install msprime
```

(It is generally better to use `python3 -m pip` rather than call `pip`
directly since this allows you to control which installation of Python the
package is installed to.)

If you do not have root access to your machine, you can install
`msprime` into your local Python installation as follows:

```
$ python3 -m pip install msprime --user
```

To use the `mspms` program you must ensure
that the `~/.local/bin` directory is in your `PATH`, or
simply run it using:

```
$ ~/.local/bin/mspms
```

To uninstall `msprime`, simply run:

```
$ python3 -m pip uninstall msprime
```

(sec_pip_install_source)=

### pip install from source

Install from source if {ref}`sec_run_pip`  fails or you do
not want to use pre-built binaries.

There may be additional requirements depending on the specific platform.
In particular, installation of [GNU Scientific
Library](<http://www.gnu.org/software/gsl/>) (GSL) is sometimes needed:

```{eval-rst}
.. glossary:: 

    Debian/Ubuntu
        ::

            $ apt-get install python-dev libgsl0-dev

    Redhat/Fedora
        ::

            $ yum install gsl-devel

    FreeBSD
        ::

            $ pkg install gsl

    OS X
        We recommend using :ref:`sec_installation_conda` to install ``msprime`` on OS X.
        However, it is also possible to install using `Homebrew <http://brew.sh/>`_:
        ::

            $ brew update
            $ brew install gsl

    Windows
        Use :ref:`sec_installation_conda`, do not install via ``pip`` on Windows.

```

With GSL installed, install from source by doing:

```
python3 -m pip install msprime --no-binary msprime
```

(sec_installation_container)=

## Via Container

An [open container](<https://opencontainers.org/>) image (aka docker image) is built on
[Dockerhub](<https://hub.docker.com/r/tskit/msprime>) for each release of
msprime. Each image is [tagged](<https://hub.docker.com/r/tskit/msprime/tags>)
with the corresponding release. For example, for msprime release 0.7.5, the
corresponding image tag is `tskit/msprime:0.7.5`.

To run a container, you can use [docker](<https://www.docker.com/>),
[Singularity](<https://sylabs.io/singularity/>),
[podman](<https://podman.io/>) or similar tools supporting docker images.

### docker

`docker` requires root privilege to run a container:

```{code-block} bash

$ sudo docker run -it tskit/msprime:<release> mspms 10 1 -T

```

### podman

podman can run an msprime container without root privilege:

```{code-block} bash

$ podman run -it docker.io/tskit/msprime:<release> mspms 10 1 -T

```

### Singularity

A docker image can also be converted to a Singularity container and
then run without root privilege:

```{code-block} bash

$ singularity pull docker://tskit/msprime:<release> msprime-<release>.simg
$ singularity exec msprime-<release>.simg mspms 10 1 -T

```

::::::{note}

It is possible that your current environment may conflict with the environment in the singularity container.
There are two workarounds:

1. Ignore your home with the conflicting environment with `--contain` or `-H </new/path/to/home> -e`

```{code-block} bash

$ singularity shell --contain msprime-release-0.7.3.simg
Singularity: Invoking an interactive shell within container...

Singularity msprime-release-0.7.3.simg:~> python3
Python 3.6.8 (default, Jan 14 2019, 11:02:34)
[GCC 8.0.1 20180414 (experimental) [trunk revision 259383]] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import msprime
>>>

```

or use a different path as your home that does not have a conflicting environment

```{code-block} bash

$ singularity shell -H </new/path/to/home> -e msprime-release-0.7.3.simg
Singularity: Invoking an interactive shell within container...

Singularity msprime-release-0.7.3.simg:~/cnn_classify_demography> python3
Python 3.6.8 (default, Jan 14 2019, 11:02:34)
[GCC 8.0.1 20180414 (experimental) [trunk revision 259383]] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import msprime
>>>

```

2. In python get rid of your local path

```{code-block} bash

$ singularity shell msprime-release-0.7.3.simg
Singularity: Invoking an interactive shell within container...

Singularity msprime-release-0.7.3.simg:~> python3
Python 3.6.8 (default, Jan 14 2019, 11:02:34)
[GCC 8.0.1 20180414 (experimental) [trunk revision 259383]] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import sys
>>> for _path in sys.path:
...     if ".local" in _path:
...             sys.path.remove(_path)
...
>>> import msprime
>>>

```

::::::

For more information on Singularity, see <https://www.sylabs.io/guides/3.6/user-guide/>


