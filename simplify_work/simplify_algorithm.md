As in the msprime implementation,
we always keep 

- $P$ : a list of *populations*, which are lists of *ancestors*, 
    which are linked lists of (left, right, node) ancestry segments, and
- $C$ : the coalescence records that will be output at the end.

We don't need the Fenwick tree $L$.

Suppose we will simplify a tree sequence of $N$ samples
by reducing it to the subset describing the ancestry of samples
$\{i_1, \ldots, i_n\}$.
Do this translates between three things:

1. node IDs in the *input* tree sequence
2. ancestor IDs in the *internal* state of the algorithm
3. node IDs in the *output* tree sequence

We just keep track of the conversion from *input* to *internal*;
the algorithm does appropriate *output*.
This means we don't end up with the conversion from *input* to *output* node IDs,
but we don't actually need this, except for the samples.
Currently, we store this conversion in a dict, $A$,
so that for an input node ID $u$, the value of $A[u]$
gives the "internal ancestor",
i.e., the index of the first ancestry segment corresponding to that ancestor.
This means that we need also to be able to find ancestral segments based on index.

**Initialization:**

We first let $P$ be a list of null ancestors,
but then set $P[i_k] = (0, m, k)$ for each $k$.


**Update:**

We iterate through parents, ordered by increasing time-in-the-past.
Suppose all the edgesets corresponding to the next parent, $u$, is
$$
\{
(\ell_1, r_1, u, (v_{11}, \ldots, v_{1k_1})),
\ldots
(\ell_n, r_n, u, (v_{n1}, \ldots, v_{nk_n})) 
\} ,
$$
and are ordered by position $\ell_1 \le r_1 \le \ell_2 \cdots \le \ell_n \le r_n$.

We then need to do the following:

* REMOVE the segments $(\ell_i, r_i)$ from all entries in $P$ corresponding to $(v_{i1}, \ldots, v_{ik})$, and
* MERGE those segments into a new entry in $P$ with parent $u$.

The index of the first such segment (corresponding to $(\ell_1, r_1)$) we then store in $A[u]$.

*Note* that we do all edges corresponding to a given parent at once
so that every MERGE operation that produces an coalescence produces a new, unique ancestor in the *output* records.

*Note* that $u$ may or may not end up in the output NodeTable,
as it may not contain any coalescences or indeed possibly any ancestry of the samples;
whether it is or not is handled in the coalescence output during the MERGE step.

*REMOVE:* 
step through each of the linked lists of ancestral segments for every offspring $v$
and delete/shorten overlapping segments; the removed sections are stored in $H$.


*MERGE:* 
existing algorithm `merge_ancestors` applied to $H$;
modified only to return the index of the first segment of $u$.


**Notes:**

msprime doesn't associate unique IDs with ancestors, so we have to maintain this mapping.
It's currently done in a dictionary, but should figure out the right way to do it.
