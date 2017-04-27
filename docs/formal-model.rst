.. _sec-file-format:

======================================================
Formal (and Algorithmic) Definition of a Tree Sequence
======================================================

Here we write down the formal requirements for a set of tables to give
a valid ARG and the algorithms to work properly on them.

Copied from elsewhere, must update:

1. Offspring must be born after their parents (and hence, no loops).
2. The set of intervals on which individual $a$ is a child must be disjoint, for every $a$.
3. The set of intervals on which individual $a$ is a parent must be disjoint, for every $a$.
4. All records with the same parent must occur at the same time.
5. The samples must be numbered 0,...,n-1, and the smallest non-sampled label must be an internal node.
6. The list of offspring in a coalescence record must be sorted.
7. Records must be sorted in nondecreasing time order.
8. Node times must be strictly greater than zero.

The first two disallow time travel and multiple inheritance;
the third and fifth are algorithmic requirements; 
and the fourth implies that we measure branch lengths in clock time
(and hence get ultrametric trees).



*********
NodeTable
*********

************
EdgesetTable
************


*************
MutationTable
*************
