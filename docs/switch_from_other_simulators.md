(sec_switch_from_other_simulators)=

# Switching from other simulators

:::{seealso}
Msprime outputs results as a succinct tree sequence using
the  [tskit](https://tskit.dev) library.
Please see {ref}`tutorials:sec_tskit_getting_started` for information
on how to use this data structure to compute statistics and
(if necessary) export to other formats.
:::

## Notes for ms users

Msprime began as an efficient reimplementation of the classical ``ms``
program, and therefore largely follows the same underlying models
as ``ms``. If you wish to use ``msprime`` as a direct replacement
for ``ms``, please use the {ref}`mspms<sec_mspms>` program. This
aims to be 100% compatible with ``ms``, and should be substantially
faster when simulating large regions.

However, simulations using ``mspms`` will still be hampered by the
text output of ``ms``, which is very inefficient for large
simulations. For larger simulations the [tskit](https://tskit.dev/tskit)
output produced by ``msprime``'s Python API will be many times
smaller and faster to process than text based formats.

The Python APIs for {ref}`describing demographic models<sec_demography>`,
running {ref}`ancestry<sec_ancestry>` and
{ref}`mutation<sec_mutations>` simulations are extensively documented.
There are a few basic differences to ``ms`` that are worth keeping
in mind:

- Time is measured in generations rather than coalescent units.
- Population sizes are defined as absolute values, not scaled relative to Ne.
- All rates are absolute rates per generation (i.e., not scaled)
- Recombination rates are per unit of sequence length (not total rates
  along the genome)
- Gene conversion rates are absolute, and do not depend on the recombination
  rate.
- Mutations are at finite discrete sites by default, following the Jukes-Cantor
  nucleotide model (but infinite sites 0/1 mutations at continuous locations
  can be produced, if required).


## Notes for discoal, msms and macs users

Different coalescent simulators may use different scaling factors for
simulation parameters, making it confusing to users who would like to switch
from simulators such as ``discoal`` ([Kern,
2016](https://academic.oup.com/bioinformatics/article/32/24/3839/2525668)),
``msms`` ([Ewing,
2010](https://academic.oup.com/bioinformatics/article/26/16/2064/216505)),
and ``macs`` ([Chen, 2009](https://genome.cshlp.org/content/19/1/136.long))
to ``msprime``. Here, we compared a few parameters and arguments commonly
used in these simulators to mitigate this pain. As ``msprime`` python API is
more efficient and expressive than ``mspms`` interface, the following
comparisons is done between the {func}`.sim_ancestry` and
{func}`.sim_mutations` interface of ``msprime``, and the command-line options
of the other simulators.

To better compare them, let us define a few terms:

- `r`: the recombination rate per generation per base pair.
- `u`: the mutation rate per generation per base pair.
- `g`: time in generations before present.
- `nsam`: sample size or the number of chromosomes simulated.
- `N`: the effective population size (of diploids).
- `L`: the length of the simulated chromosome in base pairs.
- `s`: the selection coefficient is `s` for genotype homozygous for the
selected allele (`AA`), `hs` for (`Aa`). Currently, `h` is always 0.5 in
`msprime`.
- `selpos`: the position (in base pairs) of the site under selection.

### Scaled parameters

|                       | Msprime    | Discoal    | Msms        | Macs   |
|:----------------------|:-----------|:-----------|:------------|:-------|
| Sample size           | nsam/2     | nsam       | nsam        | nsam   |
| Recombination rate    | r          | 4Nr(L-1)   | 4Nr(L-1)    | 4Nr    |
| mutation rate         | u          | 4NuL       | 4NuL        | 4Nu    |
| Selection coefficient | AA: s      | AA: 2Ns    | AA: 2Ns     | -      |
|                       | Aa: 0.5s   | Aa:  Ns    | Aa: 2Nhs    | -      |
| Selection position    | selpos     | selopos/L  | selpos/L    | -      |
| Time                  | g          | g/(4N)     | g/(4N)      | g/(2N) |

Additional comments:

- `msprime` can directly use `nsam` as sample size parameter if `ploidy` is
specified to 1. See {ref}`sec_ancestry_ploidy` section for detailed
explanation about `ploidy`.
- Time in `macs` is also calculated as g/(4N), see
[Browning, 2015](https://www.cell.com/ajhg/fulltext/S0002-9297(15)00288-8).

### Example 1

In **Example 1**, we simulate 100 samples with a sample size of 300
(chromosomes), a chromosome length of 1000,000 bp, recombination and mutation
rates of 1e-8 per generation per base pair, and an effective population size
of 1000.

#### Msprime

```python
import msprime

ts = msprime.sim_ancestry(
    samples=150,  # Default ploidy is 2, so 300 sampled chromosomes
    population_size=1000,
    recombination_rate=1e-8,
    sequence_length=1e6,
    discrete_genome=False,
    model='Hudson',
)
mutated_ts = msprime.sim_mutations(
    ts,
    rate=1e-8,
    discrete_genome=False,
    model=msprime.BinaryMutationModel(),
)
```

:::{note}
``msprime`` separates {ref}`ancestry<sec_ancestry>` and {ref}`mutation<sec_mutations>`
simulations, although the two steps are usually run together in other
simulators. See the {ref}`sec_quickstart` section for an introduction to
simulating ancestry and mutations in ``msprime``.
:::

:::{seealso}
We run a single replicate simulation here. Please
see the {ref}`sec_randomness` section for more information on
how to run replicate simulations in msprime.
:::

#### Discoal
```sh
# discoal {nsam} {num_repeats} {L} -t {4*N*u*L} -r {4*N*r*(L-1)}

$ discoal 300 100 1000000 -t 40 -r 40
```

#### Msms

```sh
# java -jar msms.jar {nsam} {num_repeats} \
#   -t {4*N*u*L} -r {4*N*r*(L-1)} {L}

$ java -Xmx1G -jar msms.jar 300 100 -t 40 -r 40 1000000
```

#### MaCS

```sh
# macs {nsam} {L} -t {4*N*u} -r {4*N*r}

$ macs 300 1000000 -t 0.00004 -r 0.00004
```

### Example 2

In **Example 2**, we add a selective sweep at the midpoint of the chromosome
(position 500,000) ended at 80 generations ago. The end frequency of the
selected allele is 0.9; the selection coefficient is 0.2 for genotype `AA` and
0.1 for `Aa`. This type of simulation can be done in ``msprime``, ``discoal``
and ``msms``.


#### Msprime

```python
import msprime

N = 1000  # effective population size

model_list = [
    msprime.StandardCoalescent(duration=80),  # From generation 0 to 80
    msprime.SweepGenicSelection(              # From generation 80 to
        position=500000,  # selpos            # selection start time (random)
        start_frequency=1 / (2 * N),
        end_frequency=0.9,
        s=0.2,  # s for AA
        dt=1.0 / (40 * N)
    ),
    msprime.StandardCoalescent()  # From selection start time to coalescence
]

ts = msprime.sim_ancestry(
    samples=150,  # Default ploidy is 2, so 300 sampled chromosomes
    population_size=N,
    recombination_rate=1e-8,
    sequence_length=1e6,
    discrete_genome=False,
    model=model_list,
)
mutated_ts = msprime.sim_mutations(
    ts,
    rate=1e-8,
    discrete_genome=False,
    model=msprime.BinaryMutationModel(),
)
```

The {class}`.SweepGenicSelection` class provides an interface to define the
selective sweep model. More examples can be found in the
{ref}`sec_ancestry_models_selective_sweeps` section.

:::{note}
``msprime`` ancestry `model` parameters can be specified as a single model
(**Example 1**), or a list of models (**Example 2**). Please
see the {ref}`sec_ancestry_models` section for more details.
:::

:::{seealso}
We run a single replicate simulation here. Please
see the {ref}`sec_randomness` section for more information on
how to run replicate simulations in msprime.
:::

#### Discoal

```sh
# discoal {nsam} {num_repeats} {L} -t {4*N*u*L} -r {4*N*r*(L-1)} \
#   -ws {sel_end_time/(4N)} -a {s * 0.5} \
#   -x {selpos/L} -c {end_frequency} -N {N}

$ discoal 300 100 1000000 -t 40 -r 40 \
   -ws 0.02 -a 200 -x 0.5 -c 0.9 -N 1000
```

#### Msms

```sh
# java -jar msms.jar {nsam} {num_repeats} -t {4*N*u*L} -r {4*N*r*(L-1)} {L} \
#   -SaA {2*N*(s*0.5)} -SAA {2*N*s} \
#   -SF {sel_end_time/(4N)} {end_frequency} \
#   -Sp {selpos/L} -N {N}

$ java -Xmx1G -jar msms.jar 300 100 -t 40 -r 40 1000000 \
  -SaA 200 -SAA 400 -SF 0.02 0.9 -Sp 0.5 -N 1000
```
