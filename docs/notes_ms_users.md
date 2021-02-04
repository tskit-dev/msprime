(sec_notes_ms_users)=
# Notes for ms users

```{eval-rst}
.. todo:: This is copied from the old api.rst page and needs some updating.
```

The simulation model in `msprime` closely follows the classical `ms`
program. Unlike `ms`, however, time is measured in generations rather than
in units of {math}`4 N_e` generations, i.e., "coalescent units".
This means that when simulating a population with diploid effective size {math}`N_e`,
the mean time to coalescence between two samples
in an `msprime` simulation will be around {math}`2 N_e`,
while in an `ms` simulation, the mean time will be around {math}`0.5`.
Internally, `msprime` uses the same algorithm as `ms`,
and so the `Ne` parameter to the {func}`.simulate` function
still acts as a time scaling, and can be set to `0.5` to match many theoretical results,
or to `0.25` to match `ms`. Population sizes for each
subpopulation and for past demographic events are also defined as absolute values, **not**
scaled by `Ne`. All migration rates and growth rates are also per generation.

:::{warning}

This parameterisation of recombination, mutation and
migration rates is different to {program}`ms`, which states these
rates over the entire region and in coalescent time units. The
motivation for this is to allow the user: 1) to change the size of the simulated
region without having to rescale the recombination, gene conversion, and mutation rates,
and 2) to directly state times and rates in units of
generations. However, the `mspms` command line application is
fully {program}`ms` compatible.
If recombination and gene conversion are combined the gene conversion
rate in {program}`ms` is determined by the ratio {math}`f`, which corresponds to
setting {math}`g = f r`. In `msprime` the gene conversion rate {math}`g` is
set independently and does not depend on the recombination rate. However,
`mspms` mimics the {program}`ms` behaviour.

:::
