---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---


(sec_logging)=

# Logging

Msprime uses the Python {mod}`logging` infrastructure to help debugging complex
simulations. Messages at the INFO level are high-level information about the
state of the current simulation, such as switches in simulation model, etc.
While it is straightforward to set up logging messages using the built-in
Python methods, the [daiquiri](<https://daiquiri.readthedocs.io/en/latest/>)
library is a bit more convenient. For example,

% Note: the output is very ugly here but there doesn't seem to be much we
% can do about until streams can be "coalesced" upstream:
% https://github.com/executablebooks/jupyter-book/issues/973
% We could try to change the buffering on stderr, I guess, but it seems
% fiddly.

```{code-cell}
import daiquiri
import msprime

daiquiri.setup(level="INFO")
ts = msprime.sim_ancestry(
    10,
    population_size=1000,
    model=[
        msprime.DiscreteTimeWrightFisher(duration=100),
        msprime.StandardCoalescent(),
    ],
    random_seed=1234
)
```

When running larger simulations and trying to figure out when
they might finish, it can be helpful to use the DEBUG logging output.
For example:

```{code-cell}

daiquiri.setup(level="DEBUG")
ts = msprime.sim_ancestry(
    10 ** 5,
    population_size=10000,
    recombination_rate=2e-8,
    sequence_length=1e6,
    random_seed=32
)
```

In this example we run a reasonably large simulation and turn on
the DEBUG output. This will then periodically (every 10,000 simulation
events) print out the current time in the simulation, and the
number of extant ancestral lineages.

:::{warning}
The format of these logging messages is not fixed and may change
arbitrarily in the future. If you need to obtain the information within
them, please open an issue on GitHub so that we can provide a documented
API for this.
:::
