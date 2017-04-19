import msprime
import math

def add_paintings(paintings, add_fun, parent, time, empty=-1):
    """
    'Add' a collection of paintings using add_fun(colors, left, right, parent,
    time) to produce the new color.
    """
    if len(paintings) == 1:
        return paintings[0]
    # remove the empty ones, they are just trouble
    paints = [p for p in paintings if len(p) > 0]
    # print(paints)
    new_paint = []
    if len(paints) == 0:
        return new_paint
    np = len(paints)
    not_done = [True for _ in paints]
    # records which interval we are on for each painting
    inds = [0 for _ in paints]
    # and whether we are inside the interval or in the space before it
    inside = [(x[0][0] == 0.0) for x in paints]
    # next breakpoint for each function (they all start at 0.0)
    nexts = [(x[0][1] if z else x[0][0]) for x,z in zip(paints, inside)]
    # current colors for each painting
    current = [(x[0][2] if z else empty) for x,z in zip(paints, inside)]
    last_break = 0.0
    while any(not_done):
        next_break = min(nexts)
        this_color = add_fun(current, left=last_break, right=next_break,
                             parent=parent, time=time)
        # print(last_break, "---", next_break, "-->", this_color)
        # print("  ", not_done, inds, inside, nexts, current)
        if this_color is not empty:
            new_paint.append([last_break, next_break, this_color])
        for k in range(np):
            if nexts[k] == next_break:
                if not inside[k]:
                    nexts[k] = paints[k][inds[k]][1]
                    inside[k] = True
                    current[k] = paints[k][inds[k]][2]
                elif inds[k] == len(paints[k]) - 1:
                    not_done[k] = False
                    nexts[k] = math.inf
                    current[k] = empty
                else:
                    inds[k] += 1
                    next_left = paints[k][inds[k]][0]
                    if next_left == next_break:
                        nexts[k] = paints[k][inds[k]][1]
                        current[k] = paints[k][inds[k]][2]
                    else:
                        inside[k] = False
                        nexts[k] = next_left
                        current[k] = empty
        last_break = next_break
    # print(new_paint)
    return new_paint

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
        # print("merging:", merge_segments)
        merge_segments = [remove_paint(painting[child], edge.left, edge.right) 
                            for child in edge.children if child in painting]
        merged_segments = add_paintings(merge_segments, add_fun=add_colors,
                                        parent=edge.parent,
                                        time=orig_nodes.time[edge.parent])
        # print("------>", merged_segments)
        if edge.parent not in painting:
            painting[edge.parent] = []
        painting[edge.parent].extend(merged_segments)
    # print(nodes)
    # print(edgesets)
    new_ts = msprime.load_tables(nodes=nodes, edgesets=edgesets)
    return new_ts

