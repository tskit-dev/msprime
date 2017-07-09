As in the msprime implementation,
we always keep 

- $P$ : a list of *populations*, which are lists of *ancestors*, 
    which are linked lists of (left, right, node) ancestry segments, and
- $C$ : the coalescence records that will be output at the end.

We don't need the two trees $S$ and $L$.

Suppose we will simplify a tree sequence of $N$ samples
by reducing it to the subset describing the ancestry of samples
$\{i_1, \ldots, i_n\}$.

**Initialization:**

We first let $P$ be a list of null ancestors,
but then set $P[i_k] = (0, m, k)$ for each $k$.

**Update:**

Iterate through edgeset records of the old tree sequence.
Suppose the next is
$$
(\ell, r, u, (v_1, \ldots, v_k)) .
$$

We then need to

* REMOVE the segment $(\ell, r)$ from all entries in $P$ corresponding to $(v_1, \ldots, v_k)$, and
* MERGE those segments into a new entry in $P$ for $u$.

*REMOVE:* 
must step through each of the appropriate linked lists
and delete overlapping segments; storing these in $H$


*MERGE:* existing algorithm, applied to $H$


**Notes:**

msprime doesn't associate unique IDs with ancestors, so we have to maintain this mapping.
It's currently done in a dictionary, but should figure out the right way to do it.
