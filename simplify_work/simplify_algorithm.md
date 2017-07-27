# Tree sequence simplification

This is an algorithm that takes a set of tables (Nodes and Edgesets, at least)
and a list of "sampled" node IDs in these tables, and produces a new set of tables
that contains all information relevant to the history of the sampled IDs, but no more.
It does this in such a way that the sampled node IDs have IDs 0, ..., n-1 in the new set of tables.

Suppose we will simplify a tree sequence of $N$ samples
by reducing it to the subset describing the ancestry of samples
$\{i_1, \ldots, i_n\}$.
Do this entails an identification between three things.

1. node IDs in the *input* tree sequence (call these *input IDs* or *input nodes*),
2. ancestors in the *internal* state of the algorithm (call these *ancestors*), and
3. node IDs in the *output* tree sequence (call these *output IDs* or *output nodes*)

There is a one-to-one mapping between *input nodes* and *ancestors*,
but not every *output node* has a corresponding ancestor: those that do not contain a coalescence
event in the reduced tree sequence are not represented.
Furthermore, each output node corresponds only to a *subset* of the ancestral segments
in one particular ancestor: those segments that have coalesced in that ancestor.
Nevertheless, the output to input ID map could be recorded explicitly, but we do not need to.

## Data structures

- an *ancestral segment* records `(left, right, node, prev, next)`, meaning it:
    * covers the region `[left, right)`
    * is a segment of ancestry inherited by *output node ID* `node`
    * is preceded by the ancestral segment `prev` in this ancestor (`NULL` if it is first)
    * and is followed by the ancestral segment `next` in this ancestor (`NULL` if it is last)

- an *ancestor* is a sorted, (linked) list of nonoverlapping ancestral
  segments, referred to by (a pointer to) the first one.

**Input:**

- `in_N` : the input node table
- `in_E` : the input edgeset table
- `samples` : the list of input nodes to preserve
- `sequence_length` : the total length of the chromosome

**Internal state:**

- `ancestor_map` : a vector of length equal to the number of nodes in `in_N` of pointers to ancestors.
    * `ancestor_map[u]` is the ancestor corresponding to input ID `u`
    * we need this so that we know which ancestors should be acted on when processing the input edgesets.

- `merge_queue` : a heap of ancestors, sorted on `left` endpoint
    * these are destined to be merged into a new `ancestor`
    * we process `in_E` in such a way as to go pull *all* ancestral segments correponding to each *input ID* out at once,
      storing them here before processing coalescence events and producing the new ancestor

- `overlap_count` : a heap of `(position, count)` pairs sorted on `position`
    * `overlap_count[pos]` is the number of remaining ancestral segments
      covering the region from position `pos` until the next position to the right
    * we need this to know when we are all done processing a given segment 
      (i.e., have reached the root locally)

**Output:**

- `out_N` : the output node table
    * `nrow(out_N)` is the number of rows in `out_N`, and hence the *next* output ID to be assigned.

- `out_E` : the output edgeset table

## Algorithm

*Overview:*
We process all input edgesets corresponding to a given (input) parent together,
removing them from the children and merging them together to make the new ancestor.
Each segment of ancestry is labeled by which *output* node label it corresponds to.
Any segment over which there was a coalescence event in a given ancestor
is labeled by the *output* node label of that ancestor,
but if an ancestor does not contain any such segments,
she will not be assigned an output label.

Output nodes will be recorded in order, backwards in time:

* `record_node(input_id)` copies the `input_id`th row of `in_N` to the end of `out_N`
* and so after this occurs, `input_id` will correspond to `nrow(out_N) - 1`


**Initialization:**

For each `0 <= k < len(samples)`:

1. `samples[k]` gets output label `k`: 
    create a new ancestral segment `a = (0, sequence_length, k, NULL, NULL)`.
2. Record this: `ancestor_map[samples[k]] = a`, and
3. `record_node(samples[k])`.


**Loop:**

Since edgesets are ordered by parent time and then by left endpoint,
we can iterate down `in_E` until an edgeset with a new parent is reached,
accumulating ancestral segments in `merge_queue`,
at which point these are merged.

Each edgeset is of the form `(left, right, parent, children)`.
For each `child` in each edgeset:

- if `ancestor_map[child]` is `NULL` then go to the next child.
- otherwise, `remove_ancestry(left, right, child)`, which moves (sub-)segments to `merge_queue`.

Upon encountering an edgeset with a different parent, we

- `merge_ancestors(parent)`, which empties `merge_queue`,

before moving on to the next parent.

These two operations are:

`remove_ancestry(left, right, child)`:

At the end, `y` will be the rightmost segment remaining in `child` before the removed segment,
`x` will be the leftmost one after, and `w` will be the next segment to be output.

0. Let `x` be the first ancestral segment in `child`, and set `y = w = NULL`.

1. While `x.right < left`, set `y = x` and `x = x.next`.

2. If `x.left < left` then we must split `x`: set `y` to be the part of x before `left`

    and set `x` to be a new segment that is the remaining part.
3. While `x.left < right`, output segments:

    * allocate `next_w` to be a new segment with `(x.left, min(x.right, right), x.node, w, None)`
    * if `w` is `NULL` then insert `next_w` into `merge_queue`.
    * if `x.right <= right`, set `x = x.next` and move on
    * otherwise, set `x.left = right`

4. Wrap up, setting `y.next = x` and `x.prev = y` if these are not `NULL`;
    and updating `ancestor_map[child] = x` if `y` is `NULL`.


*MERGE:* 

`merge_ancestors(parent)` will take all ancestors in `merge_queue`
and merge them into a single, new ancestor who is recorded in `ancestor_map` under `parent`.
(Note that these are not ancestors who have a corresponding input ID,
but they are linked lists of ancestral segments.)
In the process, any overlapping segments produce a coalescence event, so:

    * copy the `parent`th row of `in_N` to the bottom of `out_N`
    * add a new edgeset to `out_E` with parent `nrow(out_N)-1` and children given
      by the node labels on the contributing ancestral segments
    * reduce `overlap_count` on this segment appropriately
    * if `overlap_count` tells us no more coalescences on this segment are possible,
      then do not add the merged segment to the new ancestor
    * otherwise, add the segment, with node label `nrow(out_N)-1`.


**Updating mutations:**

We can update the `node` entries in the mutation table as we go along, then go back and cull out unused `sites` afterwards.

At each input edgeset, the input node ID associated with a particular line of ancestry changes,
so we need to deal with any associated mutations.  To do this, we just need to assign any mutations associated
with any input nodes encountered with the corresponding output node ID.  Since each ancestry
segment carries with it the output node ID, this is easy.
If
- there is no coalescence, we just find corresponding mutations, assign the output node ID, and output.
- there is coalescence and there are remaining overlapping segments, assign the (new) output node ID, and output
- there is coalescence and no remaining overlapping segments, then find the node's state and assign this to the site's ancestral state.

