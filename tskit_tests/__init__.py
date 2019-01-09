from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import base64

# TODO remove this code and refactor elsewhere.

from .simplify import *  # NOQA

import tskit


class PythonTree(object):
    """
    Presents the same interface as the Tree object for testing. This
    is tightly coupled with the PythonTreeSequence object below which updates
    the internal structures during iteration.
    """
    def __init__(self, num_nodes):
        self.num_nodes = num_nodes
        self.parent = [tskit.NULL for _ in range(num_nodes)]
        self.left_child = [tskit.NULL for _ in range(num_nodes)]
        self.right_child = [tskit.NULL for _ in range(num_nodes)]
        self.left_sib = [tskit.NULL for _ in range(num_nodes)]
        self.right_sib = [tskit.NULL for _ in range(num_nodes)]
        self.above_sample = [False for _ in range(num_nodes)]
        self.is_sample = [False for _ in range(num_nodes)]
        self.left = 0
        self.right = 0
        self.root = 0
        self.index = -1
        self.left_root = -1
        # We need a sites function, so this name is taken.
        self.site_list = []

    @classmethod
    def from_tree(cls, tree):
        ret = PythonTree(tree.num_nodes)
        ret.left, ret.right = tree.get_interval()
        ret.site_list = list(tree.sites())
        ret.index = tree.get_index()
        ret.left_root = tree.left_root
        ret.tree = tree
        for u in range(ret.num_nodes):
            ret.parent[u] = tree.parent(u)
            ret.left_child[u] = tree.left_child(u)
            ret.right_child[u] = tree.right_child(u)
            ret.left_sib[u] = tree.left_sib(u)
            ret.right_sib[u] = tree.right_sib(u)
        assert ret == tree
        return ret

    @property
    def roots(self):
        u = self.left_root
        roots = []
        while u != tskit.NULL:
            roots.append(u)
            u = self.right_sib[u]
        return roots

    def children(self, u):
        v = self.left_child[u]
        ret = []
        while v != tskit.NULL:
            ret.append(v)
            v = self.right_sib[v]
        return ret

    def _preorder_nodes(self, u, l):
        l.append(u)
        for c in self.children(u):
            self._preorder_nodes(c, l)

    def _postorder_nodes(self, u, l):
        for c in self.children(u):
            self._postorder_nodes(c, l)
        l.append(u)

    def _inorder_nodes(self, u, l):
        children = self.children(u)
        if len(children) > 0:
            mid = len(children) // 2
            for v in children[:mid]:
                self._inorder_nodes(v, l)
            l.append(u)
            for v in children[mid:]:
                self._inorder_nodes(v, l)
        else:
            l.append(u)

    def _levelorder_nodes(self, u, l, level):
        l[level].append(u) if level < len(l) else l.append([u])
        for c in self.children(u):
            self._levelorder_nodes(c, l, level + 1)

    def nodes(self, root=None, order="preorder"):
        roots = [root]
        if root is None:
            roots = self.roots
        for u in roots:
            node_list = []
            if order == "preorder":
                self._preorder_nodes(u, node_list)
            elif order == "inorder":
                self._inorder_nodes(u, node_list)
            elif order == "postorder":
                self._postorder_nodes(u, node_list)
            elif order == "levelorder" or order == "breadthfirst":
                # Returns nodes in their respective levels
                # Nested list comprehension flattens node_list in order
                self._levelorder_nodes(u, node_list, 0)
                node_list = iter([i for level in node_list for i in level])
            else:
                raise ValueError("order not supported")
            for v in node_list:
                yield v

    def get_interval(self):
        return self.left, self.right

    def get_parent(self, node):
        return self.parent[node]

    def get_children(self, node):
        return self.children[node]

    def get_index(self):
        return self.index

    def get_parent_dict(self):
        d = {
            u: self.parent[u] for u in range(self.num_nodes)
            if self.parent[u] != tskit.NULL}
        return d

    def sites(self):
        return iter(self.site_list)

    def __eq__(self, other):
        return (
            self.get_parent_dict() == other.get_parent_dict() and
            self.get_interval() == other.get_interval() and
            self.roots == other.roots and
            self.get_index() == other.get_index() and
            list(self.sites()) == list(other.sites()))

    def __ne__(self, other):
        return not self.__eq__(other)

    def newick(self, root=None, precision=16, node_labels=None):
        if node_labels is None:
            node_labels = {u: str(u + 1) for u in self.tree.leaves()}
        if root is None:
            root = self.left_root
        return self._build_newick(root, precision, node_labels) + ";"

    def _build_newick(self, node, precision, node_labels):
        label = node_labels.get(node, "")
        if self.left_child[node] == tskit.NULL:
            s = label
        else:
            s = "("
            for child in self.children(node):
                branch_length = self.tree.branch_length(child)
                subtree = self._build_newick(child, precision, node_labels)
                s += subtree + ":{0:.{1}f},".format(branch_length, precision)
            s = s[:-1] + label + ")"
        return s


class PythonTreeSequence(object):
    """
    A python implementation of the TreeSequence object.
    """
    def __init__(self, tree_sequence, breakpoints=None):
        self._tree_sequence = tree_sequence
        self._num_samples = tree_sequence.get_num_samples()
        self._breakpoints = breakpoints
        self._sites = []

        def make_mutation(id_):
            site, node, derived_state, parent, metadata = tree_sequence.get_mutation(id_)
            return tskit.Mutation(
                id_=id_, site=site, node=node, derived_state=derived_state,
                parent=parent, metadata=metadata)
        for j in range(tree_sequence.get_num_sites()):
            pos, ancestral_state, ll_mutations, id_, metadata = tree_sequence.get_site(j)
            self._sites.append(tskit.Site(
                id_=id_, position=pos, ancestral_state=ancestral_state,
                mutations=[make_mutation(ll_mut) for ll_mut in ll_mutations],
                metadata=metadata))

    def edge_diffs(self):
        M = self._tree_sequence.get_num_edges()
        sequence_length = self._tree_sequence.get_sequence_length()
        edges = [tskit.Edge(*self._tree_sequence.get_edge(j)) for j in range(M)]
        time = [self._tree_sequence.get_node(edge.parent)[1] for edge in edges]
        in_order = sorted(range(M), key=lambda j: (
            edges[j].left, time[j], edges[j].parent, edges[j].child))
        out_order = sorted(range(M), key=lambda j: (
            edges[j].right, -time[j], -edges[j].parent, -edges[j].child))
        j = 0
        k = 0
        left = 0
        while j < M or left < sequence_length:
            e_out = []
            e_in = []
            while k < M and edges[out_order[k]].right == left:
                h = out_order[k]
                e_out.append(edges[h])
                k += 1
            while j < M and edges[in_order[j]].left == left:
                h = in_order[j]
                e_in.append(edges[h])
                j += 1
            right = sequence_length
            if j < M:
                right = min(right, edges[in_order[j]].left)
            if k < M:
                right = min(right, edges[out_order[k]].right)
            yield (left, right), e_out, e_in
            left = right

    def trees(self):
        M = self._tree_sequence.get_num_edges()
        sequence_length = self._tree_sequence.get_sequence_length()
        edges = [
            tskit.Edge(*self._tree_sequence.get_edge(j)) for j in range(M)]
        t = [
            self._tree_sequence.get_node(j)[1]
            for j in range(self._tree_sequence.get_num_nodes())]
        in_order = sorted(
            range(M), key=lambda j: (
                edges[j].left, t[edges[j].parent], edges[j].parent, edges[j].child))
        out_order = sorted(
            range(M), key=lambda j: (
                edges[j].right, -t[edges[j].parent], -edges[j].parent, -edges[j].child))
        j = 0
        k = 0
        N = self._tree_sequence.get_num_nodes()
        st = PythonTree(N)

        samples = list(self._tree_sequence.get_samples())
        for l in range(len(samples)):
            if l < len(samples) - 1:
                st.right_sib[samples[l]] = samples[l + 1]
            if l > 0:
                st.left_sib[samples[l]] = samples[l - 1]
            st.above_sample[samples[l]] = True
            st.is_sample[samples[l]] = True

        st.left_root = tskit.NULL
        if len(samples) > 0:
            st.left_root = samples[0]

        u = st.left_root
        roots = []
        while u != -1:
            roots.append(u)
            v = st.right_sib[u]
            if v != -1:
                assert st.left_sib[v] == u
            u = v

        st.left = 0
        while j < M or st.left < sequence_length:
            while k < M and edges[out_order[k]].right == st.left:
                p = edges[out_order[k]].parent
                c = edges[out_order[k]].child
                k += 1

                lsib = st.left_sib[c]
                rsib = st.right_sib[c]
                if lsib == tskit.NULL:
                    st.left_child[p] = rsib
                else:
                    st.right_sib[lsib] = rsib
                if rsib == tskit.NULL:
                    st.right_child[p] = lsib
                else:
                    st.left_sib[rsib] = lsib
                st.parent[c] = tskit.NULL
                st.left_sib[c] = tskit.NULL
                st.right_sib[c] = tskit.NULL

                # If c is not above a sample then we have nothing to do as we
                # cannot affect the status of any roots.
                if st.above_sample[c]:
                    # Compute the new above sample status for the nodes from
                    # p up to root.
                    v = p
                    above_sample = False
                    while v != tskit.NULL and not above_sample:
                        above_sample = st.is_sample[v]
                        u = st.left_child[v]
                        while u != tskit.NULL:
                            above_sample = above_sample or st.above_sample[u]
                            u = st.right_sib[u]
                        st.above_sample[v] = above_sample
                        root = v
                        v = st.parent[v]

                    if not above_sample:
                        # root is no longer above samples. Remove it from the root list.
                        lroot = st.left_sib[root]
                        rroot = st.right_sib[root]
                        st.left_root = tskit.NULL
                        if lroot != tskit.NULL:
                            st.right_sib[lroot] = rroot
                            st.left_root = lroot
                        if rroot != tskit.NULL:
                            st.left_sib[rroot] = lroot
                            st.left_root = rroot
                        st.left_sib[root] = tskit.NULL
                        st.right_sib[root] = tskit.NULL

                    # Add c to the root list.
                    # print("Insert ", c, "into root list")
                    if st.left_root != tskit.NULL:
                        lroot = st.left_sib[st.left_root]
                        if lroot != tskit.NULL:
                            st.right_sib[lroot] = c
                        st.left_sib[c] = lroot
                        st.left_sib[st.left_root] = c
                    st.right_sib[c] = st.left_root
                    st.left_root = c

            while j < M and edges[in_order[j]].left == st.left:
                p = edges[in_order[j]].parent
                c = edges[in_order[j]].child
                j += 1

                # print("insert ", c, "->", p)
                st.parent[c] = p
                u = st.right_child[p]
                lsib = st.left_sib[c]
                rsib = st.right_sib[c]
                if u == tskit.NULL:
                    st.left_child[p] = c
                    st.left_sib[c] = tskit.NULL
                    st.right_sib[c] = tskit.NULL
                else:
                    st.right_sib[u] = c
                    st.left_sib[c] = u
                    st.right_sib[c] = tskit.NULL
                st.right_child[p] = c

                if st.above_sample[c]:
                    v = p
                    above_sample = False
                    while v != tskit.NULL and not above_sample:
                        above_sample = st.above_sample[v]
                        st.above_sample[v] = st.above_sample[v] or st.above_sample[c]
                        root = v
                        v = st.parent[v]
                    # print("root = ", root, st.above_sample[root])

                    if not above_sample:
                        # Replace c with root in root list.
                        # print("replacing", root, "with ", c ," in root list")
                        if lsib != tskit.NULL:
                            st.right_sib[lsib] = root
                        if rsib != tskit.NULL:
                            st.left_sib[rsib] = root
                        st.left_sib[root] = lsib
                        st.right_sib[root] = rsib
                        st.left_root = root
                    else:
                        # Remove c from root list.
                        # print("remove ", c ," from root list")
                        st.left_root = tskit.NULL
                        if lsib != tskit.NULL:
                            st.right_sib[lsib] = rsib
                            st.left_root = lsib
                        if rsib != tskit.NULL:
                            st.left_sib[rsib] = lsib
                            st.left_root = rsib

            st.right = sequence_length
            if j < M:
                st.right = min(st.right, edges[in_order[j]].left)
            if k < M:
                st.right = min(st.right, edges[out_order[k]].right)
            assert st.left_root != tskit.NULL
            while st.left_sib[st.left_root] != tskit.NULL:
                st.left_root = st.left_sib[st.left_root]
            st.index += 1
            # Add in all the sites
            st.site_list = [
                site for site in self._sites if st.left <= site.position < st.right]
            yield st
            st.left = st.right


class MRCACalculator(object):
    """
    Class to that allows us to compute the nearest common ancestor of arbitrary
    nodes in an oriented forest.

    This is an implementation of Schieber and Vishkin's nearest common ancestor
    algorithm from TAOCP volume 4A, pg.164-167 [K11]_. Preprocesses the
    input tree into a sideways heap in O(n) time and processes queries for the
    nearest common ancestor between an arbitary pair of nodes in O(1) time.

    :param oriented_forest: the input oriented forest
    :type oriented_forest: list of integers
    """
    LAMBDA = 0

    def __init__(self, oriented_forest):
        # We turn this oriened forest into a 1 based array by adding 1
        # to everything
        converted = [0] + [x + 1 for x in oriented_forest]
        self.__preprocess(converted)

    def __preprocess(self, oriented_forest):
        """
        Preprocess the oriented forest, so that we can answer mrca queries
        in constant time.
        """
        n = len(oriented_forest)
        child = [self.LAMBDA for i in range(n)]
        parent = [self.LAMBDA for i in range(n)]
        sib = [self.LAMBDA for i in range(n)]
        self.__lambda = [0 for i in range(n)]
        self.__pi = [0 for i in range(n)]
        self.__tau = [0 for i in range(n)]
        self.__beta = [0 for i in range(n)]
        self.__alpha = [0 for i in range(n)]
        for u in range(n):
            v = oriented_forest[u]
            sib[u] = child[v]
            child[v] = u
            parent[u] = v
        p = child[self.LAMBDA]
        n = 0
        self.__lambda[0] = -1
        while p != self.LAMBDA:
            notDone = True
            while notDone:
                n += 1
                self.__pi[p] = n
                self.__tau[n] = self.LAMBDA
                self.__lambda[n] = 1 + self.__lambda[n >> 1]
                if child[p] != self.LAMBDA:
                    p = child[p]
                else:
                    notDone = False
            self.__beta[p] = n
            notDone = True
            while notDone:
                self.__tau[self.__beta[p]] = parent[p]
                if sib[p] != self.LAMBDA:
                    p = sib[p]
                    notDone = False
                else:
                    p = parent[p]
                    if p != self.LAMBDA:
                        h = self.__lambda[n & -self.__pi[p]]
                        self.__beta[p] = ((n >> h) | 1) << h
                    else:
                        notDone = False
        # Begin the second traversal
        self.__lambda[0] = self.__lambda[n]
        self.__pi[self.LAMBDA] = 0
        self.__beta[self.LAMBDA] = 0
        self.__alpha[self.LAMBDA] = 0
        p = child[self.LAMBDA]
        while p != self.LAMBDA:
            notDone = True
            while notDone:
                a = (
                    self.__alpha[parent[p]] |
                    (self.__beta[p] & -self.__beta[p])
                )
                self.__alpha[p] = a
                if child[p] != self.LAMBDA:
                    p = child[p]
                else:
                    notDone = False
            notDone = True
            while notDone:
                if sib[p] != self.LAMBDA:
                    p = sib[p]
                    notDone = False
                else:
                    p = parent[p]
                    notDone = p != self.LAMBDA

    def get_mrca(self, x, y):
        """
        Returns the most recent common ancestor of the nodes x and y,
        or -1 if the nodes belong to different trees.

        :param x: the first node
        :param y: the second node
        :return: the MRCA of nodes x and y
        """
        # WE need to rescale here because SV expects 1-based arrays.
        return self._sv_mrca(x + 1, y + 1) - 1

    def _sv_mrca(self, x, y):
        if self.__beta[x] <= self.__beta[y]:
            h = self.__lambda[self.__beta[y] & -self.__beta[x]]
        else:
            h = self.__lambda[self.__beta[x] & -self.__beta[y]]
        k = self.__alpha[x] & self.__alpha[y] & -(1 << h)
        h = self.__lambda[k & -k]
        j = ((self.__beta[x] >> h) | 1) << h
        if j == self.__beta[x]:
            xhat = x
        else:
            ell = self.__lambda[self.__alpha[x] & ((1 << h) - 1)]
            xhat = self.__tau[((self.__beta[x] >> ell) | 1) << ell]
        if j == self.__beta[y]:
            yhat = y
        else:
            ell = self.__lambda[self.__alpha[y] & ((1 << h) - 1)]
            yhat = self.__tau[((self.__beta[y] >> ell) | 1) << ell]
        if self.__pi[xhat] <= self.__pi[yhat]:
            z = xhat
        else:
            z = yhat
        return z


def base64_encode(metadata):
    """
    Returns the specified metadata bytes object encoded as an ASCII-safe
    string.
    """
    return base64.b64encode(metadata).decode('utf8')
