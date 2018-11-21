.. _sec_provenance:

==========
Provenance
==========

Every tree sequence has provenance information associated with it. The purpose of this
information is to improve `reproducibility <https://en.wikipedia.org/wiki/Reproducibility>`_:
given the provenance associated with a given tree sequence, it should be possible to
reproduce it. Provenance is split into three sections: the primary **software** used to
produce a tree sequence; the **parameters** provided to this software; and the computational
**environment** where the software was run.

This documentation serves two distinct purposes:

1. For developers using ``tskit`` in their own applications, it provides normative documentation
   for how provenance information should be stored.
2. For end-users of ``tskit``, it provides documentation to allows them to inspect and interpret
   the provenance information stored in ``.trees`` files.

Provenance information is encoded using `JSON <https://www.json.org/>`_.
To standardise the provenance information produced by different software and improve
interoperability we define a formal specification using `JSON Schema <http://json-schema.org/>`_.
The full schema is provided :ref:`below <sec_provenance_schema>`, which may be used to
automatically validate input. In the following we describe the intention of the various
sections in more detail.

This document defines specification version 1.0.0. Specification version numbers follow
`SemVer <https://semver.org/>`_ semantics.

.. _sec_provenance_example:

*******
Example
*******

To make things more concrete, let's consider an example:

.. code-block:: json

    {
      "schema_version": "1.0.0",
      "software": {
        "name": "msprime",
        "version": "0.6.1.dev123+ga252341.d20180820"
      },
      "parameters": {
        "sample_size": 5,
        "random_seed": 12345,
        "command": "simulate"
      },
      "environment": {
        "libraries": {
          "gsl": {
            "version": "2.1"
          },
          "kastore": {
            "version": "0.1.0"
          }
        },
        "python": {
          "version": "3.5.2",
          "implementation": "CPython"
        },
        "os": {
          "system": "Linux",
          "node": "powderfinger",
          "release": "4.15.0-29-generic",
          "version": "#31~16.04.1-Ubuntu SMP Wed Jul 18 08:54:04 UTC 2018",
          "machine": "x86_64"
        }
      }
    }

This information records the provenance for a very simple msprime simulation. The record is a JSON
object with three mandatory fields ("software", "parameters" and "environment")
which we discuss separately in the following sections.

.. _sec_provenance_software:

********
Software
********

Every tree sequence is produced by some piece of software. For example, this may be a
coalescent simulation produced by ``msprime``, a forwards-time simulation from ``SLiM``
or tree sequence inferred from data by ``tsinfer``. The software provenance is
intended to capture the details about this primary software.

================    ==============      ===========
Field               Type                Description
================    ==============      ===========
name                string              The name of the software.
version             string              The software version.
================    ==============      ===========


Note that libraries that the primary software links against are considered part of the
:ref:`sec_provenance_environment` and should be recorded there.

.. _sec_provenance_parameters:

**********
Parameters
**********

The parameters section of a provenance document records the input that was used to
produce a particular tree sequence. There are no requirements on what may be stored
within it, but we make some recommendations here on how to encode such information.

As a general principle, sufficient information should be recorded in the parameters
section to allow the output tree sequence to be reproduced exactly. There will be instances,
however, where this is not possible due to missing files, issues with numerical precision
and so on.

+++++++++++++++
API invocations
+++++++++++++++

Consider an API call like the following simple msprime simulation:

.. code-block:: python

    ts = msprime.simulate(sample_size=10, recombination_rate=2)

We recommend encoding the parameters provenance as follows (other fields omitted
for clarity):

.. code-block:: json

    {
      "parameters": {
        "command": "simulate",
        "sample_size": 10,
        "recombination_rate": 2,
        "random_seed": 123456789,
      }
    }

Specifically, we encode the name of the function using the ``command`` key and
the function parameters in the obvious way. Note that we include the ``random_seed``
here even though it was automatically generated.


+++++++++++++++
CLI invocations
+++++++++++++++

Consider the following invocation of a hypothetical command line program:

.. code-block:: bash

    $ supersim --sample-size=10 --do-some-stuff -O out.trees

We recommend encoding the parameters provenance as follows (other fields omitted
for clarity):

.. code-block:: json

    {
      "parameters": {
        "command": "supersim",
        "args": ["--sample-size=10", "--do-some-stuff", "-O", "out.trees"],
        "random_seed": 56789
      }
    }

Here we encode the name of the program using the ``command`` key
and its command line arguments as a list of strings in the ``args`` key. We
also include the automatically generated random seed in the parameters list.

If parameters that affect the output tree sequence are derived from environment
variables these should also be recorded.

.. _sec_provenance_environment:

***********
Environment
***********

The environment section captures details about the computational environment in
which the software was executed. Two optional fields are defined: ``os``
and ``libraries``. We recommend including any additional relevant platform
information here; for example, if using Python store the interpreter information
as shown in the example above.

++++++++++++++++
Operating system
++++++++++++++++

The ``os`` section records details about the operating system on which the
software was executed. This section is optional and has no required internal
structure. We recommend the following structure based on the output of the
POSIX `uname <http://pubs.opengroup.org/onlinepubs/009695399/functions/uname.html>`_
function:

.. code-block:: json

    {
      "environment": {
        "os": {
          "system": "Linux",
          "node": "powderfinger",
          "release": "4.15.0-29-generic",
          "version": "#31~16.04.1-Ubuntu SMP Wed Jul 18 08:54:04 UTC 2018",
          "machine": "x86_64"
        }
    }


+++++++++
Libraries
+++++++++

The ``libraries`` section captures information about important libraries that the
primary software links against. There is no required structure.


.. _sec_provenance_schema:

***********
Full schema
***********

This schema is formally defined using `JSON Schema <http://json-schema.org/>`_ and
given in full here. Developers writing provenance information to ``.trees`` files
should validate the output JSON against this schema.

.. literalinclude:: ../msprime/provenance.schema.json
    :language: json
