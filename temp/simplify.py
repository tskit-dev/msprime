import msprime
import math

def add_stepfuns(funs, add_fun, parent, time, empty=-1):
    """
    "add" a collection of stepfunctions
    using add_fun(values, left, right, time)
    """
    nfuns = len(funs)
    if nfuns == 1:
        return funs[0]
    breaks = [0.0]
    values = []
    # records which interval we are on for each function
    inds = [0 for _ in funs]
    # next breakpoint for each function (they all start at 0.0)
    nexts = [x[0][1] for x in funs]
    # current values for each function
    current = [x[1][0] for x in funs]
    not_done = [True for k in range(nfuns)]
    last_break = 0.0
    while any(not_done):
        next_break = min(nexts)
        breaks.append(next_break)
        values.append(add_fun(current, last_break, next_break, parent, time))
        for k in range(nfuns):
            if nexts[k] == next_break:
                if inds[k] == len(funs[k][1]) - 1:
                    not_done[k] = False
                else:
                    inds[k] += 1
                    nexts[k] = funs[k][0][inds[k]+1]
                    current[k] = funs[k][1][inds[k]]
        last_break = next_break
    return (breaks, values)

def intervals_to_stepfun(intervals, empty=-1):
    """
    Given a list of (left,right,value) tuples, with left >= 0,
    return a tuple (breaks,values), where the function takes value[k] 
    on interval [breaks[k],breaks[k+1]), with breaks[0] = 0.0 always.  
    Any interval not specified in intervals takes the value 'empty'
    (including the last semi-infinite interval).
    """
    # print("i to st", intervals)
    breaks = [0.0]
    values = []
    for left, right, value in intervals:
        if left > breaks[-1]:
            breaks.append(left)
            values.append(empty)
        breaks.append(right)
        values.append(value)
    breaks.append(math.inf)
    values.append(empty)
    return (breaks,values)

def stepfun_to_intervals(stepfun, empty=-1):
    out = []
    for k in range(len(stepfun[1])):
        if stepfun[1][k] is not empty:
            out.append((stepfun[0][k], stepfun[0][k+1], stepfun[1][k]))
    return out

def remove_paint(painting, left, right):
    """
    Remove (modifying in place) and return the subset of the painting 
    lying within the interval (left, right).  A painting is a list of 
    nonoverlapping tuples of (left, right, color) in order by left endpoint.
    """
    new_painting = []
    k = 0
    # print("remove_paint:", painting, left, right)
    while k < len(painting):
        seg_left, seg_right, seg_color = painting[k]
        # print("  :", k, seg_left, seg_right, seg_color)
        if right < seg_left:
            break
        if left < seg_right and right > seg_left:
            out_left = max(left,seg_left)
            out_right = min(right, seg_right)
            new_painting += [(out_left, out_right, seg_color)]
            if seg_left < out_left:
                painting[k] = (seg_left, out_left, seg_color)
                k += 1
            else:
                del painting[k]
            if seg_right > out_right:
                painting.insert(k, (out_right, seg_right, seg_color))
                break
        else:
            k += 1
    # print("...", painting)
    # print("-->", new_painting)
    return new_painting

def simplify(ts, samples):
    """
    Steps through records, painting each according to who descends from them.
    When two colors combine, they get a new color and we record a record
    in the simplified tables.
    """
    nodes = msprime.NodeTable()
    edgesets = msprime.EdgesetTable()
    # dict of old labels -> new labels
    node_colors = {}

    def add_colors(values, left, right, parent, time, empty=-1):
        """
        Return a new color if there is more than one nonempty one,
        storing these events as we go.
        """
        colors = list(set(values).difference([empty]))
        # print(colors)
        if len(colors) == 0:
            return empty
        elif len(colors) == 1:
            return colors[0]
        else:
            if parent in node_colors:
                new_color = node_colors[parent]
            else:
                new_color = nodes.num_rows
                node_colors[parent] = new_color
                nodes.add_row(time=time)
            colors.sort()
            edgesets.add_row(left=left, right=right, parent=new_color, children=tuple(colors))
            return new_color

    # the 'painting' is a dict of paintings indexed by nodes in the input
    painting = {}
    for k in samples:
        new_color = nodes.num_rows
        painting[k] = [(0, ts.sequence_length, new_color)]
        nodes.add_row(time=0.0, flags=1)
        node_colors[k] = new_color
    orig_nodes = msprime.NodeTable()
    ts.dump_tables(nodes=orig_nodes)
    for edge in ts.edgesets():
        # print(edge)
        # here need to remove all painted segments overlapping this segment from children
        # and insert merged painting into parent
        merge_segments = [remove_paint(painting[child], edge.left, edge.right) 
                            for child in edge.children if child in painting]
        # print("merging:", merge_segments)
        merge_stepfuns = [intervals_to_stepfun(x) for x in merge_segments]
        # for f in merge_stepfuns:
        #     print("- ", f)
        merged_stepfun = add_stepfuns(merge_stepfuns, add_fun=add_colors, 
                                      parent=edge.parent,
                                      time=orig_nodes.time[edge.parent])
        # print("-->", merged_stepfun)
        merged_segments = stepfun_to_intervals(merged_stepfun)
        # print("------>", merged_segments)
        if edge.parent not in painting:
            painting[edge.parent] = []
        painting[edge.parent].extend(merged_segments)
    # print(nodes)
    # print(edgesets)
    new_ts = msprime.load_tables(nodes=nodes, edgesets=edgesets)
    return new_ts


###
# attempt #1 - using intervals only
###

def add_paint_intervals(painting, additions, pot):
    """
    Step through the intervals in the two paintings:
    * removing everything from old
    * if old does not overlap with anything in new, add old to new.
    * if old does overlap, add a new interval with a new color from the pot to new
    """
    k = 0
    paint_left = painting[k].left
    add_k = 0
    while add_k < len(additions):
        while k < len(painting) and painting[k].right <= additions[add_k].left:
            # [-paint-]
            #           [--add--]
            k += 1
        if k == len(painting):
            # nothing overlapping left in painting
            painting.extend(additions)
            break
        if additions[add_k].right <= painting[k].left:
            # [-add-]
            #         [--paint--]
            painting.insert(k,(additions[add_k].left, additions[add_k].right, additions[add_k].color))
            k += 1
            add_k += 1
            continue
        if additions[add_k].left > painting[k].left:
            #     [-add-?
            # [--paint--?
            painting.insert(k+1, (additions[add_k].left, painting[k].right, painting[k].color))
            painting[k].right = additions[add_k].left
            k += 1
            continue
        elif additions[add_k].left < painting[k].left:
            # [----add---?
            #   [--paint-?
            painting.insert(k, (additions[add_k].left, painting[k].left, additions[add_k].color))
            k += 1
            additions[add_k].left = painting[k].left
            # DO NOT continue as we haven't updated additions[add_k]
        # ok now we are finally in this situation:
        # [--add---?
        # [--paint-?
        merge_color = next(pot)
        merge_right = min(painting[k].right, additions[add_k].right)
        # new record to output!
        merges.append((painting[k].left, merge_right, merge_color))
        if painting[k].right > additions[add_k].right:
            # [--add---]
            # [--paint----]
            painting[k].left = additions[add_k].right
            painting.insert(k, (additions[add_k].left, additions[add_k].right, next(pot)))
            k += 1
            add_k += 1
            continue
        if painting[k].right > add_right:
            # [-----add-----]
            # [--paint--]
            k += 1
            additions[add_k].left = painting[k].right
            continue
