# Updating the UI to remove lots of `get_`s

Process:

1. wrap old functions in new API functions, without documentations
2. ... use them
3. update tests
4. swap implementation into the new functions, and make old ones deprecated in the docs
5. add deprecation warnings (but don't turn on till later versions)

# Making these changes

1. A vim regex to fix methods with non-trivial arguments:
```
'<,'>s,\(\s\+\)\(def\s\)get_\(.*(\)self\,\s\(.*)\):,\1\2\3self\, \4:\r\1\1return self.get_\3\4\r\r\1\2get_\3self\, \4:,gc
```
NOTE: this will not properly treat functions with default args --- these must be updated by hand

2. A regex to apply the `@property` properly to methods with out args other than self:
```
'<,'>s,\(\s\+\)\(def\s\)get_\(.*(\)self):,\1@property\r\1\2\3self):\r\1\1return self.get_\3)\r\r\1\2get_\3self):,gc
```

## Initial targets

NB: below is the initial set of targets for the changes. To list all possible, use
`git grep -p 'def \(get.*\)(' msprime`.

### In `SparseTree`

```
msprime/trees.py=class SparseTree(object):
msprime/trees.py:    def get_root(self):
msprime/trees.py:    def get_index(self):
msprime/trees.py:    def get_interval(self):
msprime/trees.py:    def get_length(self):
msprime/trees.py:    def get_sample_size(self):
msprime/trees.py:    def get_num_mutations(self):
msprime/trees.py:    def get_parent_dict(self):
msprime/trees.py:    def get_time_dict(self):
msprime/trees.py:    def get_total_branch_length(self):
msprime/trees.py:    def get_branch_length(self, u):
msprime/trees.py:    def get_mrca(self, u, v):
msprime/trees.py:    def get_tmrca(self, u, v):
msprime/trees.py:    def get_parent(self, u):
msprime/trees.py:    def get_children(self, u):
msprime/trees.py:    def get_time(self, u):
msprime/trees.py:    def get_population(self, u):
msprime/trees.py:    def get_num_leaves(self, u):
msprime/trees.py:    def get_num_tracked_leaves(self, u):
```

### In `TreeSequence`

```
msprime/trees.py=class TreeSequence(object):
msprime/trees.py:    def get_ll_tree_sequence(self):
msprime/trees.py:    def get_provenance(self):
msprime/trees.py:    def get_sample_size(self):
msprime/trees.py:    def get_sequence_length(self):
msprime/trees.py:    def get_num_records(self):
msprime/trees.py:    def get_num_trees(self):
msprime/trees.py:    def get_num_mutations(self):
msprime/trees.py:    def get_num_nodes(self):
msprime/trees.py:    def get_pairwise_diversity(self, samples=None):
msprime/trees.py:    def get_time(self, sample):
msprime/trees.py:    def get_population(self, sample):
msprime/trees.py:    def get_samples(self, population_id=None):
```

