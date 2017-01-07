import msprime


def path_length(tr,x,y):
    L = 0
    mrca = tr.mrca(x,y)
    for u in x,y:
        while u != mrca:
            L += tr.branch_length(u)
            u = tr.parent(u)
    return L

def branch_length_diversity(ts,X,Y):
    '''
    Computes average pairwise diversity between a random choice from x
    and a random choice from y.
    '''
    S = 0
    for tr in ts.trees():
        SS = 0
        for x in X:
            for y in Y:
                SS += path_length(tr,x,y)
        S += SS*tr.length
    return S/(ts.sequence_length*len(X)*len(Y))

def branch_length_Y(ts,x,y,z):
    S = 0
    for tr in ts.trees():
        xy_mrca = tr.mrca(x,y)
        xz_mrca = tr.mrca(x,z)
        yz_mrca = tr.mrca(y,z)
        if xy_mrca == xz_mrca:
            #   /\
            #  / /\
            # x y  z
            S += path_length(tr,x,yz_mrca)*tr.length
        elif xy_mrca == yz_mrca:
            #   /\
            #  / /\
            # y x  z
            S += path_length(tr,x,xz_mrca)*tr.length
        elif xz_mrca == yz_mrca:
            #   /\
            #  / /\
            # z x  y
            S += path_length(tr,x,xy_mrca)*tr.length
    return S/ts.sequence_length

def branch_length_f4(ts,A,B,C,D):
    S = 0
    for tr in ts.trees():
        SS = 0
        for a in A:
            for b in B:
                for c in C:
                    for d in D:
                        SS += path_length(tr,tr.mrca(a,c),tr.mrca(b,d))
                        SS -= path_length(tr,tr.mrca(a,d),tr.mrca(b,c))
        S += SS*tr.length
    return S/(ts.sequence_length*len(A)*len(B)*len(C)*len(D))

def branch_stats_node_iter(ts,leaf_sets,weight_fun,method='length'):
    '''
    Here leaf_sets is a list of lists of leaves, and weight_fun is a function
    whose argument is a list of integers of the same length as leaf_sets
    that returns a number.  Each branch in a tree is weighted by weight_fun(x),
    where x[i] is the number of leaves in leaf_sets[i] below that
    branch.  This finds the sum of all counted branches for each tree,
    and averages this across the tree sequence ts, weighted by genomic length.

    If method='mutations' instead, then branch lengths will be measured in
    numbers of mutations instead of time.

    This version is inefficient as it iterates over all nodes in each tree.
    '''
    return branch_stats_vector_node_iter(ts,leaf_sets,lambda x: [weight_fun(x)],method)[0]

def branch_stats(ts,leaf_sets,weight_fun,method='length'):
    '''
    Here leaf_sets is a list of lists of leaves, and weight_fun is a function
    whose argument is a list of integers of the same length as leaf_sets
    that returns a boolean.  A branch in a tree is weighted by weight_fun(x),
    where x[i] is the number of leaves in leaf_sets[i] below that
    branch.  This finds the sum of all counted branches for each tree,
    and averages this across the tree sequence ts, weighted by genomic length.

    Doesn't do method='mutations'.
    '''
    return branch_stats_vector(ts,leaf_sets,lambda x: [weight_fun(x)],method)[0]

def branch_stats_vector_node_iter(ts,leaf_sets,weight_fun,method='length'):
    '''
    Here leaf_sets is a list of lists of leaves, and weight_fun is a function
    whose argument is a list of integers of the same length as leaf_sets
    that returns a list of numbers; there will be one output for each element.
    For each value, each branch in a tree is weighted by weight_fun(x),
    where x[i] is the number of leaves in leaf_sets[i] below that
    branch.  This finds the sum of all counted branches for each tree,
    and averages this across the tree sequence ts, weighted by genomic length.

    If method='mutations' instead, then branch lengths will be measured in
    numbers of mutations instead of time.

    This version is inefficient as it iterates over all nodes in each tree.
    '''
    tr_its = [ ts.trees(tracked_leaves=x,leaf_counts=True,leaf_lists=True) for x in leaf_sets ]
    n_out = len(weight_fun([0 for a in leaf_sets]))
    S = [ 0.0 for j in range(n_out) ]
    for k in range(ts.num_trees):
        trs = [ next(x) for x in tr_its ]
        root = trs[0].root
        tr_len = trs[0].length
        if method=='length':
            for node in trs[0].nodes():
                if node != root:
                    x = [ tr.num_tracked_leaves(node) for tr in trs ]
                    w = weight_fun(x)
                    for j in range(n_out):
                        S[j] += w[j] * trs[0].branch_length(node) * tr_len
        elif method=='mutations':
            count_nodes = dict([ 
                (node,weight_fun([ tr.num_tracked_leaves(node) for tr in trs ])) 
                for node in trs[0].nodes() if node != root ])
            # print(count_nodes)
            for mut in trs[0].mutations():
                # print(mut)
                for j in range(n_out):
                    S[j] += count_nodes[mut.node][j]
        else:
            raise(TypeError("Unknown method "+method))
    for j in range(n_out):
        S[j] /= ts.get_sequence_length()
    return S

def branch_stats_vector(ts,leaf_sets,weight_fun,method='length'):
    '''
    Here leaf_sets is a list of lists of leaves, and weight_fun is a function
    whose argument is a list of integers of the same length as leaf_sets
    that returns a boolean.  A branch in a tree is weighted by weight_fun(x),
    where x[i] is the number of leaves in leaf_sets[i] below that
    branch.  This finds the sum of all counted branches for each tree,
    and averages this across the tree sequence ts, weighted by genomic length.

    Doesn't do method='mutations'.
    '''
    # initialize
    num_leaf_sets = len(leaf_sets)
    n_out = len(weight_fun([0 for a in range(num_leaf_sets)]))
    # print("n_out:",n_out)
    S = [ 0.0 for j in range(n_out) ]
    L = [ 0.0 for j in range(n_out) ]
    N = ts.num_nodes
    X = [ [ int(u in a) for a in leaf_sets] for u in range(N) ]
    # we will essentially construct the tree
    pi = [-1 for j in range(N)]
    node_time = [0.0 for u in range(N)]
    for length,records_out,records_in in ts.diffs():
        for sign,records in ((-1,records_out), (+1,records_in)):
            for node,children,time in records:
                # print("Record (",sign,"):",node,children,time)
                # print("\t",X, "-->", L)
                if sign==+1:
                    node_time[node] = time
                dx = [0 for k in range(num_leaf_sets)]
                for child in children:
                    if sign==+1:
                        pi[child] = node
                    for k in range(num_leaf_sets):
                        dx[k] += sign * X[child][k]
                    w = weight_fun(X[child])
                    dt = (node_time[pi[child]] - node_time[child])
                    for j in range(n_out):
                        L[j] += sign * dt * w[j]
                    # print("\t\tchild:",child,sign,weight_fun(X[child]),node_time[pi[child]],node_time[child],"-->",L)
                    if sign==-1:
                        pi[child] = -1
                old_w = weight_fun(X[node])
                for k in range(num_leaf_sets):
                    X[node][k] += dx[k]
                if pi[node] != -1:
                    w = weight_fun(X[node])
                    dt = (node_time[pi[node]] - node_time[node])
                    for j in range(n_out):
                        L[j] += dt * (w[j]-old_w[j])
                    # print("\t\tnode:",node,weight_fun(X[node]),old_w,"-->",L)
                # propagate change up the tree
                u = pi[node]
                if u != -1:
                    next_u = pi[u]
                    while u != -1:
                        old_w = weight_fun(X[u])
                        for k in range(num_leaf_sets):
                            X[u][k] += dx[k]
                        # need to update X for the root,
                        # but the root does not have a branch length
                        if next_u != -1:
                            w = weight_fun(X[u])
                            dt = (node_time[pi[u]] - node_time[u])
                            for j in range(n_out):
                                L[j] += dt*(w[j] - old_w[j])
                            # print("\t\tanc:",u,weight_fun(X[u]),old_w,"-->",L)
                        u = next_u
                        next_u = pi[next_u]
                # print("\t",X, "-->", L)
        # print("next tree:",L,length)
        for j in range(n_out):
            S[j] += L[j] * length
    for j in range(n_out):
        S[j] /= ts.get_sequence_length()
    return S


