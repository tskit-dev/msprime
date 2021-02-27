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

```{code-cell}
:tags: [remove-cell]

    import msprime
    from IPython.display import SVG
    import numpy as np
    # Doing this to make the notebook outputs deterministic. 
    # DO NOT DO THIS IN YOUR CODE
    msprime.core.set_seed_rng_seed(42)
```

(sec_mutations)=

# Mutations


:::{warning}

This documentation is under heavy construction. Please note
any outstanding TODOs before opening issues.

:::

---

## Quick reference


{func}`.sim_mutations`
: Add simulated mutations to a tree sequence

**Models**

{class}`.BinaryMutationModel`             
: Binary mutation model with two flip-flopping alleles: "0" and "1".

{class}`.JC69MutationModel`               
: Jukes & Cantor model ('69), equal probability of transitions between nucleotides

{class}`.HKYMutationModel`                
: Hasegawa, Kishino & Yano model ('85), different probabilities for transitions and transversions

{class}`.F84MutationModel`                
: Felsenstein model ('84), different probabilities for transitions and transversions

{class}`.GTRMutationModel`                
: Generalised Time-Reversible nucleotide mutation model

{class}`.BLOSUM62MutationModel`                
: The BLOSUM62 model of time-reversible amino acid mutation

{class}`.PAMMutationModel`              
: The PAM model of time-reversible amino acid mutation

{class}`.MatrixMutationModel`             
: Superclass of mutation models with a finite set of states

{class}`.InfiniteAllelesMutationModel`    
: A generic infinite-alleles mutation model          

{class}`.SLiMMutationModel`               
: An infinite-alleles model producing SLiM-style mutations

---

{func}`.sim_mutations` adds mutations to a tree sequence simulated with
{func}`.sim_ancestry`, or any other program that outputs tree sequence topologies.
We can specify a uniform mutation rate or a
{class}`.RateMap`, whether or not mutations should only occur at 
{ref}`discrete coordinates <sec_mutations_discrete>`,
a {ref}`specific time span within which to add mutations <sec_mutations_time_span>`,
as well as how to handle any
{ref}`existing mutations in a tree sequence <sec_mutations_existing>`.

One of the most powerful features of `msprime` is the flexibility of specifying
the model under which mutations are simulated. Select from among
the many predefined {ref}`models <sec_mutations_models>`, such as the 
{ref}`infinite alleles model <sec_mutations_mutation_infinite_alleles>`, or
{ref}`define your own <sec_mutations_matrix_mutation_models_details>`.


(sec_mutations_rate)=

## Specifying rates

The `rate` parameter of {func}`.sim_mutations` determines the rate of mutations per unit
of sequence length per generation. A rate must be specified for mutations to be generated.

```{code-cell}
ts = msprime.sim_ancestry(2, random_seed=1)
mts = msprime.sim_mutations(ts, rate=1, random_seed=1)
mts.num_mutations
```

Running {func}`.sim_mutations` on this simulated tree of length 1 with `rate=1`
produces 7 mutations.

Specifying higher rates lead to proportionately more mutations:

```{code-cell}
ts = msprime.sim_ancestry(2, random_seed=1)
mts = msprime.sim_mutations(ts, rate=5, random_seed=1)
mts.num_mutations
```

It's also possible to provide a {class}`.RateMap` which specifies variable
mutation rates over different stretches of the sequence. 

In the following example, the mutation rate between positions 0 and 2 is 0.5,
between positions 2 and 5 is 0, and between positions 5 and 10 is 0.1:

```{code-cell}
ts = msprime.sim_ancestry(4, sequence_length=7, random_seed=1)
pos = [0, 2, 5, 7]
rate = [0.5, 0, 0.1]
ratemap = msprime.RateMap(pos, rate)
mts = msprime.sim_mutations(ts, rate=ratemap, random_seed=10)
SVG(mts.draw_svg(node_labels={}, size=(400, 300)))
```

As we can see from the output, there are three mutations at site 0 (position 0),
three mutations at site 1 (position 1), no mutations between positions 2 and 5 
(where the mutation rate is zero) and one mutation at site 2 (position 5).
This illustrates that as the rate increases, the probability of recurrent mutation
(multiple mutations at a site) also increases.

:::{note}

If using a {class}`.RateMap` with {func}`.sim_mutations`, be sure that the final
position of the {class}`.RateMap` is the same as the
`sequence_length` of the tree sequence you're adding mutations to!

:::


(sec_mutations_randomness)=

## Controlling randomness

As described in {ref}`sec_ancestry_random_seed`, ``msprime`` uses a
[psuedorandom number generator](<https://en.wikipedia.org/wiki/Pseudorandom_number_generator>)
from the [GNU Scientific Library](<https://www.gnu.org/software/gsl/doc/html/rng.html>) for
random number generation. Passing an integer to the ``random_seed`` parameter 
defines a trajectory for a call of {func}`.sim_mutations`, making output deterministic.

```{code-cell}
ts = msprime.sim_ancestry(10, random_seed=1)
mts_1 = msprime.sim_mutations(ts, rate=1, random_seed=7)
mts_2 = msprime.sim_mutations(ts, rate=1, random_seed=8)
mts_1.tables.equals(mts_2.tables, ignore_timestamps=True)
```

Here, we check for equality between the TableCollections of the `mts_1` and `mts_2` tree
sequences, excluding the differing timestamps in Provenance tables. The result is
`False` because we used different seeds.

```{code-cell}
mts_3 = msprime.sim_mutations(ts, rate=1, random_seed=7)
mts_1.tables.equals(mts_3.tables, ignore_timestamps=True)
```

When we use the same seed, the resulting tree sequence is identical. 


(sec_mutations_discrete)=

## Discrete or continuous

As in {func}`.sim_ancestry` (see {ref}`sec_ancestry_discrete_genome`),
the `discrete_genome` parameter controls whether mutations should be placed at
discrete, integer coordinates or continously at floating point positions.
`discrete_genome=True` is the default:

```{code-cell}
ts = msprime.sim_ancestry(4, sequence_length=5, random_seed=1) 
mts = msprime.sim_mutations(ts, rate=0.1, random_seed=5)
mts.tables.sites.position
```

Specifying `discrete_genome=False` places mutations at floating point positions
continously along the genome.

```{code-cell}
mts = msprime.sim_mutations(ts, rate=0.1, random_seed=5, discrete_genome=False)
mts.tables.sites.position
```

Note that using `discrete_genome=False` is equivalent to specifying an infinite sites
model.


(sec_mutations_time_span)=

## Restricting time span

The time range where mutations can be placed is specified using the `start_time`
and `end_time` parameters. For instance, the following only allows mutations to occur
earlier than time ``1`` in the tree:

```{code-cell}
ts = msprime.sim_ancestry(2, random_seed=1)
mts = msprime.sim_mutations(ts, rate=2, start_time=1, random_seed=1)
mts.tables.mutations.time
```

Note, however, that the child node of the edge where the mutation occured can be younger
than `start_time`.

```{code-cell}
mts.tables.nodes.time[mts.tables.mutations.node]
```

It is also possible to use multiple calls of {func}`.sim_mutations` with `start_time`
and `end_time` specified to simulate differing mutation rates at different
periods of time. For instance, the following code simulates mutations at a low rate
prior to time ``1`` and at a higher rate more recently.

```{code-cell}
ts = msprime.sim_ancestry(2, random_seed=1)
mts = msprime.sim_mutations(ts, rate=0.1, start_time=1, random_seed=1)
print(mts.tables.mutations.time)
mts = msprime.sim_mutations(mts, rate=4, end_time=1, random_seed=1)
print(mts.tables.mutations.time)
```

As explained in {ref}`the following section <sec_mutations_existing>`, reversing the order of these two
lines will result in an error, as older mutations must be added first.

(sec_mutations_existing)=

## Existing mutations

If adding mutations to a tree sequence which already contains them, the `keep` parameter
controls whether existing mutations are kept or discarded (the default is `keep=True`).
For instance, in final code block in {ref}`sec_mutations_time_span`, mutations were
progressively added to a simulated tree sequence, beginninng with the oldest time period.

While it is simple to add younger mutations to a tree sequence which already contains
mutations, adding a mutation ancestral to an existing mutation will result in
an error unless `add_ancestral=True` is specified, as in the following code block:

```{code-cell}
ts = msprime.sim_ancestry(5, random_seed=1)
mts = msprime.sim_mutations(ts, rate=1, random_seed=1)
mts = msprime.sim_mutations(mts, rate=0.1, random_seed=5, add_ancestral=True)
```

:::{note}

Note that specifying `add_ancestral=True` means that "silent transitions" may occur.
A silent transition is a mutation which causes a nucleotide to transition to 
itself, such as "A -> A".

:::

Adding ancestral mutations can be problematic, so `add_ancestral` should be specified
with care. For instance, if mutations are added to
a tree sequence under two different models, impossible transitions may result. See
{func}`.sim_mutations` for further explanation.


(sec_mutations_models)=

## Models

Mutation models are specified using the `model` parameter to
{func}`.sim_mutations`. This parameter can either take the form of a
string describing the model (e.g. `model="jc69"`) or an instance of a
model definition class (e.g `model=msprime.JC69MutationModel()`).
Here are the available models; they are documented in more detail below.

- {class}`.BinaryMutationModel`: Basic binary mutation model with two flip-flopping alleles: "0" and "1".
- {class}`.JC69MutationModel`: Jukes & Cantor model ('69), equal probability of transitions between nucleotides
- {class}`.HKYMutationModel`: Hasegawa, Kishino & Yano model ('85), different probabilities for transitions and transversions
- {class}`.F84MutationModel`: Felsenstein model ('84), different probabilities for transitions and transversions
- {class}`.GTRMutationModel`: The Generalised Time-Reversible nucleotide mutation model, a general parameterization of a time-reversible mutation process
- {class}`.BLOSUM62MutationModel`: The BLOSUM62 model of time-reversible amino acid mutation
- {class}`.PAMMutationModel`: The PAM model of time-reversible amino acid mutation
- {class}`.MatrixMutationModel`: Superclass of the specific mutation models with a finite set of states
- {class}`.InfiniteAllelesMutationModel`: A generic infinite-alleles mutation model
- {class}`.SLiMMutationModel`: An infinite-alleles model of mutation producing SLiM-style mutations

(sec_mutations_matrix_mutations_models)=

### Matrix Mutation Models

These classes are defined by an alphabet of possible alleles (`alleles`); an array of
probabilities that determines how likely each allele is to be the root, ancestral allele
(`root_distribution`); and a transition matrix specifying the probability for each allele
to mutate to every other allele. Each class has specific values of these parameters to
create the specific model. For your own custom model these parameters can be set using
{class}`msprime.MatrixMutationModel`. For more detail about how mutations are simulated
in these models see {ref}`sec_mutations_matrix_mutation_models_details`.

(sec_mutations_matrix_mutation_models_details)=

### Mutation Matrix Models Details

Mutation matrix models are specified by three things: an alphabet,
a root distribution, and a transition matrix.
These leave one free parameter: an overall mutation rate,
specified by the mutation `rate` in the call to {func}`.sim_mutations`.
Concisely,
the underlying model of mutation is a continuous-time Markov chain on the alphabet,
started by a draw from `root_distribution`, and
with instantaneous transition rate from `i` to `j` that is equal to
`rate` multiplied by `transition_matrix[i,j]`.
The `root distribution` and every row in the `transition_matrix`
must give *probabilities*, i.e., they must be nonnegative numbers summing to 1.
For the precise interpretation of these parameters
(especially when the transition matrix has nonzero entries on the diagonal)
see {ref}`sec_mutations_matrix_mutation_theory`.

You can define your own, but you probably don't need to:
there are several mutation matrix models already implemented in `msprime`,
using binary (0/1), nucleotide, or amino acid alphabets:

### Defining your own finite-sites model

If you want to define your own {class}`.MatrixMutationModel`, you have a good
deal of freedom. For instance, here's a "decomposition/growth/disturbance"
mutation model, where the only possible transitions are 🎄 to 🔥, 🔥 to 💩, and
💩 to 🎄, with the first transition happening at one-fifth the rate of the
other two:

```{code-cell} python

alleles = ["💩", "🎄", "🔥"]
model = msprime.MatrixMutationModel(
    alleles,
    root_distribution = [1.0, 0.0, 0.0],
    transition_matrix = [[0.0, 1.0, 0.0],
                         [0.0, 0.8, 0.2],
                         [1.0, 0.0, 0.0]]
)
ts = msprime.sim_ancestry(6, population_size=10, random_seed=2, sequence_length=7)
mts = msprime.sim_mutations(ts, rate=2, random_seed=1, model=model)

```

We have simulated from this model at rate 2, so the overall rate of mutation
from 💩 to 🎄 and 🔥 to 💩 is 2, and from 🎄 to 🔥 is {math}`2 \times 0.2 = 0.4`.
As a result, roughly 5/7th of the states will be 🎄, with the remainder
divided evenly between 💩 and 🔥. Here is the resulting "genotype matrix":

```{code-cell}

for v in mts.variants():
   print("".join(v.alleles[k] for k in v.genotypes))


```

(sec_mutations_matrix_mutation_theory)=

### Parameterization of Matrix Mutation Models

Mutation matrix models are specified by three things: an alphabet,
a root distribution, and a transition matrix.
These leave one free parameter: an overall mutation rate,
specified by the mutation `rate` in the call to {func}`.sim_mutations`.
Concisely,
the underlying model of mutation is a continuous-time Markov chain on the alphabet,
started by a draw from `root_distribution`, and
with instantaneous transition rate from `i` to `j` that is equal to
`rate` multiplied by `transition_matrix[i,j]`.
The `root distribution` and every row in the `transition_matrix`
must give *probabilities*, i.e., they must be nonnegative numbers summing to 1.

To interpret these parameters,
it helps to know how the underlying mutational process is implemented.
First, "possible" mutations are placed on the tree,
with a mean density equal to the `rate`, per unit of time and sequence length.
If `discrete_genome=False` then this is an infinite-sites model,
so each possible mutation occurs at a distinct location.
If `discrete_genome=True` (the default setting) then at each integer position,
each branch of the tree at that position gets a Poisson number of mutations
with mean equal to `rate` multiplied by the length of the branch.
Next, each site that has a possible mutation is assigned an ancestral state,
i.e., the allele at the root of the tree at that position,
by drawing an allele from the probabilities in the `root_distribution`.
Now, each possible mutation is examined, moving down the tree.
For each, a derived state is chosen using the probabilities given in the
row of the `transition_matrix` that corresponds to the "parental state",
i.e., the allele that this mutation will replace.
Importantly, if the chosen allele is the *same* as the parental allele,
no mutation is recorded (that's why they were called "possible mutations").
And, any site at which no mutations are recorded is not recorded either.

This arrangement is necessary to fully specify Markov models of mutation,
with a free "mutation rate" parameter.
However, there are some surprising consequences.
For instance, the distribution of ancestral alleles, across all sites,
is *not* necessarily equal to the root distribution.
This is because the root distribution gives the distribution of
"ancestral" alleles across the entire sequence,
but we only see the ancestral alleles at *mutated* sites,
and some alleles may have a higher mutation rate than others.
For instance, if we have

```{code-block} python

import numpy as np

alleles = ["A", "C", "G", "T"]
root_distribution = np.array([0.25, 0.25, 0.25, 0.25])
transition_matrix = np.array([
   [0.25, 0.25, 0.25, 0.25],
   [ 0.3,  0.0,  0.4,  0.3],
   [ 0.3,  0.4,  0.0,  0.3],
   [0.25, 0.25, 0.25, 0.25]
])

```

then A and T alleles have a 25% lower mutation rate than do C and G alleles,
since 25% of the time that we consider mutating them, we leave them unchanged.
From the properties of the Poisson distribution,
the probability that a tree of total length {math}`T`
has no mutations at a given discrete site is {math}`\exp(-rT)`,
if mutations are put down at a rate of {math}`r`.
Suppose that a single tree of total length {math}`T = 1.5`
extends over many discrete sites,
and that mutations are placed on it at rate {math}`r = 2.0`.
Every site that is assigned a "C" or "G" ancestral allele is retained,
but of those sites that are assigned an "A" or "T",
some are not recorded in the resulting tree sequence.
The expected proportions of the ancestral states
across all sites is proportional to the root distribution
multiplied by the probability that at least one mutation is retained on the tree.
In this situation it can be computed as follows:

```{code-block} python

r = 2.0
T = 1.5
prob_mut = 1.0 - np.diag(transition_matrix)
ancestral_state_distribution = root_distribution * (1 - np.exp(- r * T * prob_mut))
ancestral_state_distribution /= sum(ancestral_state_distribution)

```

Two more facts about Markov chains are useful to interpret the statistics
of these mutation models.
First, suppose we have tabulated all mutations, and so for each pair of alleles
{math}`i` and {math}`j` we have the proportion of mutations that caused an {math}`i \to j` change.
If allele {math}`i` mutates to a different allele, the chance it mutates to allele {math}`j`
is proportional to `transition_matrix[i,j]` but excluding the diagonal (no-change) entry,
so is equal to `transition_matrix[i,j] / (1 - transition_matrix[i,i])`.
Second, suppose that an ancestor carries allele {math}`i` at a given position.
The probability that her descendant some time {math}`t` in the future carries allele {math}`j` at that position
is given by a matrix exponential of
the scaled [infinitestimal rate matrix](<https://en.wikipedia.org/wiki/Transition_rate_matrix>) of the Markov chain,
which can be computed as follows:

```{code-block} python
import scipy

Q = (transition_matrix - np.eye(len(alleles)))
Pt = scipy.linalg.expm(t * rate * Q)[i,j]

```

If the top of a branch of length {math}`t` has allele {math}`i`,
the bottom of the branch has allele {math}`j` with probability `Pt[i,j]`.

(sec_mutations_mutation_infinite_alleles)=

### Infinite Alleles Mutation Models

You can also use a model of *infinite alleles* mutation: where each new mutation produces a unique,
never-before-seen allele. The underlying mutation model just assigns the derived state
to be a new integer every time a new mutation appears.
By default these integers start at zero, but a different starting point can be chosen,
with the `start_allele` parameter.
It does this globally across all mutations, so that the first assigned allele will be `start_allele`,
and if `n` alleles are assigned in total (across ancestral and derived states),
these will be the next `n-1` integers.
Many theoretical results are derived based on this mutation model (e.g., Ewens' sampling formula).

For instance, here we'll simulate with the infinite alleles model on a single tree,
and print the resulting tree, labeling each mutation with its derived state:

```{code-cell} python

ts = msprime.sim_ancestry(6, random_seed=2, sequence_length=1)
model = msprime.InfiniteAllelesMutationModel()
mts = msprime.sim_mutations(ts, rate=2, random_seed=1, model=model)
t = mts.first()
ml = {m.id: m.derived_state for m in mts.mutations()}
SVG(t.draw_svg(mutation_labels=ml, node_labels={}, size=(400, 300)))

```

Apparently, there were 20 mutations at this site, but the alleles present in the population are
"13" (in five copies), "17" (in two copies), and one copy each of "14", "15", "19", and "20".
Note that all other mutations on the tree are not observed in the population as they have
been are "overwritten" by subsequent mutations.


:::{warning}

Neither this nor the next infinite alleles mutation model check to see if the alleles
they produce already exist at the mutated sites. So, if you are using these
models to add mutations to an already-mutated tree sequence, it is up to you
to set the starting allele appropriately, and to make sure the results make sense!

:::

(sec_mutations_mutation_slim_mutations)=

### SLiM mutations

A special class of infinite alleles model is provided for use with [SLiM](<https://messerlab.org/slim/>),
to agree with the underlying mutation model in SLiM.
As with the InfiniteAlleles model, it assigns each new mutation a unique integer,
by keeping track of the `next_id` and incrementing it each time a new mutation appears.

This differs from the {class}`.InfiniteAllelesMutationmodel` because mutations
in SLiM can "stack": new mutations can add to the existing state, rather than
replacing the previous state. So, derived states are comma-separated lists of
mutation IDs, and the ancestral state is always the empty string. For instance,
if a new mutation with ID 5 occurs at a site, and then later another mutation
appears with ID 64, the sequence of alleles moving along this line of descent
would be `""`, then `"5"`, and finally `"5,64"`. Furthermore, the mutation
model adds SLiM metadata to each mutation, which records, among other things,
the SLiM mutation type of each mutation, and the selection coefficient (which
is always 0.0, since adding mutations in this way only makes sense if they are
neutral). For this reason, the model has one required parameter: the `type`
of the mutation, a nonnegative integer. If, for instance, you specify
`type=1`, then the mutations in SLiM will be of type `m1`. For more
information, and for how to modify the metadata (e.g., changing the selection
coefficients), see
[the pyslim documentation](<https://pyslim.readthedocs.io/en/latest/>).
For instance,

```{code-cell} python

model = msprime.SLiMMutationModel(type=1)
mts = msprime.sim_mutations(
    ts, rate=1, random_seed=1, model=model, add_ancestral=True)
t = mts.first()
ml = {m.id: m.derived_state for m in mts.mutations()}
SVG(t.draw_svg(mutation_labels=ml, node_labels={}, size=(400, 300)))

```

These resulting alleles show how derived states are built.

The behavior of this mutation model when used to add mutations to a previously mutated
tree sequence can be subtle. Let's look at a simple example.
Here, we first lay down mutations of type 1, starting from ID 0:

```{code-cell} python

model_1 = msprime.SLiMMutationModel(type=1)
mts_1 = msprime.sim_mutations(ts, rate=0.5, random_seed=2, model=model_1)
t = mts_1.first()
ml = {m.id: m.derived_state for m in mts_1.mutations()}
SVG(t.draw_svg(mutation_labels=ml, node_labels={}, size=(400, 300)))

```

Next, we lay down mutations of type 2.
These we assign starting from ID 100,
to make it easy to see which are which:
in general just need to make sure that we start at an ID greater than any
previously assigned.
Note that without the `add_ancestral=True` parameter
this would cause an error, because we are adding mutations above existing ones.

```{code-cell} python

model_2 = msprime.SLiMMutationModel(type=2, next_id=100)
mts = msprime.sim_mutations(
    mts_1, rate=0.5, random_seed=3, model=model_2, add_ancestral=True, keep=True)
t = mts.first()
ml = {m.id: m.derived_state for m in mts.mutations()}
SVG(t.draw_svg(mutation_labels=ml, node_labels={}, size=(400, 300)))

```


Note what has happened here: on the top branch on the right side of the tree,
with the first model we added two mutations: first a mutation with ID `0`,
then a mutation with ID `3`.
Then, with the second model, we added two more mutations to this same branch,
with IDs `100` and `102`, between these two mutations.
These were added to mutation `0`, obtaining alleles `0,100` and `0,100,102`.
But then, moving down the branch, we come upon the mutation with ID `3`.
This was already present in the tree sequence, so its derived state is not modified:
`0,3`. We can rationalize this, post-hoc, by saying that the type 1 mutation `3`
has "erased" the type 2 mutations `100` and `102`.
If you want a different arrangment,
you can go back and edit the derived states (and metadata) as you like.

