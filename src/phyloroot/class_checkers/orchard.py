from copy import deepcopy


def get_leaves_and_root(network):
    """Returns the root and the set of leaves of the network"""
    leaves = set()
    for node in network:
        if network.out_degree(node) == 0:
            leaves.add(node)
        if network.in_degree(node) == 0:
            root = node
    return (root, leaves)


def is_orchard(network):
    """A function that returns True iff the given network is orchard."""
    reducing_network = deepcopy(network)
    root, leaves = get_leaves_and_root(network)

    if reducing_network.out_degree(root) == 2:
        reducing_network.add_edge(-2, root)
    checked_all_leaves = False
    while not checked_all_leaves:
        # if we get through
        reduction_done = False
        for leaf in leaves:
            pair = is_second_in_pair(reducing_network, leaf)
            if pair:
                reduced = reduce_pair(reducing_network, *pair)
                if reduced == "C":
                    leaves.remove(pair[0])
                reduction_done = True
                break
        if len(reducing_network.edges) == 1:
            return True
        # if there was no reduction after checking all the leaves,
        # then we can reduce the network no further
        checked_all_leaves = not reduction_done
    return False


def reduce_pair(network, x, y):
    """Reduces the pair (x,y) in network if it is reducible in network.

    Changes the network in place, but also returns the reduced network.
    """
    for node in network.predecessors(x):
        px = node
    for node in network.predecessors(y):
        py = node
    if px == py:
        for node in network.predecessors(py):
            ppy = node
        network.remove_edges_from([(ppy, py), (py, x), (py, y)])
        network.add_edge(ppy, y)
        return "C"
    if py in network.predecessors(px):
        for node in network.predecessors(py):
            ppy = node
        for node in network.predecessors(px):
            if node != py:
                ppx = node
        network.remove_edges_from([(ppy, py), (py, px), (py, y), (ppx, px), (px, x)])
        network.add_edges_from([(ppy, y), (ppx, x)])
        return "RC"
    return False


def is_second_in_pair(network, x):
    """Returns a reducible pair of network with x as second element if it exists, returns False otherwise."""
    for node in network.predecessors(x):
        px = node
    cpx = False
    for node in network.successors(px):
        if node != x:
            cpx = node
    if type(cpx) == bool and not cpx:
        return False
    if network.out_degree(cpx) == 0:
        return (cpx, x)
    if network.out_degree(cpx) == 1:
        for node in network.successors(cpx):
            if network.out_degree(node) == 0:
                return (node, x)
    return False
