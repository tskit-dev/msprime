.. _sec-file-format:

======================================================
Formal (and Algorithmic) Definition of a Tree Sequence
======================================================

This document records the formal requirements for a set of tables to give a
valid ARG and the algorithms to work properly on them.

To disallow time travel and multiple inheritance:

1. Offspring must be born after their parents (and hence, no loops).
2. The set of intervals on which individual $a$ is a child must be disjoint,
   for every $a$.

and for algorithmic reasons:

3. The set of intervals on which individual $a$ is a parent must be disjoint,
   for every $a$.
4. The list of offspring in a coalescence record must be sorted.
5. Records must be sorted in nondecreasing time order.
6. Node times must be strictly greater than zero.

Note that each node time is equal to the (birth) time of the corresponding parent.
This implies that time is measured in clock time (not meioses).


*********
NodeTable
*********

(insert documentation for NodeTable here)

************
EdgesetTable
************

(insert documentation for EdgesetTable here)

*************
MutationTable
*************

(insert documentation for MutationTable here)

