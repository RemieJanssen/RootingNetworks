from phyloroot.class_checkers import is_network
from phyloroot.chain_reduction import reduce_chains
from phyloroot.rooting_exponential import c_orientation_exponential, c_orientation_exponential_cycle_basis

def get_leaf_adjacent_to_node(network, node):
    """Returns the leaf node adjacent to a node in the network, if it exists.

    Args:
        network (nx.Graph): The network
        node (int): The node in the network

    Returns:
        int or None: The leaf node adjacent to the node, or None if no such node exists.
    """
    for nb in network.neighbors(node):
        if network.degree(nb) == 1:
            return nb
    return None

def expand_reduced_rootings_for_short_side(
    reduced_network, reduced_rootings, reduced_side
):
    """Expands the rootings of a reduced network to the original network for sides of length <= ell+2

    Args:
        reduced_network (nx.Graph): The reduced network
        reduced_rootings (dict(tuple(int))): A dict with all valid C-root-edges of the reduced network as keys, and a tuple of
            reticulation nodes of a valid C-orientation with this root edge of `reduced_network`.
        reduced_side (list(int)): The list of nodes in the reduced side

    Returns:
        dict(tuple(int)): A dict with all valid C-root-edges on the side of the original network corresponding to `reduced_side`
        as keys, and a tuple of reticulation nodes of a valid C-orientation with this root edge of `network`.
    """

    all_side_edges = [
        (reduced_side[i], reduced_side[i + 1]) for i in range(len(reduced_side) - 1)
    ]
    all_side_leaf_edges = [
        (node, nb)
        for node in reduced_side[1:-1]
        for nb in reduced_network.neighbors(node)
        if reduced_network.degree(nb) == 1
    ]
    rootings = dict()
    for edge in all_side_edges + all_side_leaf_edges:
        rootingAtEdge = reduced_rootings.get(edge)
        if not rootingAtEdge:
            rootingAtEdge = reduced_rootings.get((edge[1], edge[0]))
        if rootingAtEdge:
            rootings[edge] = rootingAtEdge
    return rootings


def expand_reduced_rootings_for_long_side(
    network, reduced_network, reduced_rootings, reduced_side, side, ell
):
    """Expands the rootings of a reduced network to the original network for sides of length > ell+2

    Uses the fact that the rootings of the network are determined by the rootings at the leaf edges of the reduced side.

    Args:
        network (nx.Graph): The original network
        reduced_network (nx.Graph): The reduced network
        reduced_rootings (dict(tuple(int))): A dict with all valid C-root-edges of the reduced network as keys, and a tuple of
            reticulation nodes of a valid C-orientation with this root edge of `reduced_network`.
        reduced_side (list(int)): The list of nodes in the reduced side
        side (list(int)): The list of nodes in the original side
        ell (int): The length to which chains are reduced.

    Returns:
        dict(tuple(int)): A dict with all valid C-root-edges on the side of the original network corresponding to `reduced_side`
        as keys, and a tuple of reticulation nodes of a valid C-orientation with this root edge of `network`.
    """

    rootings = dict()
    # Infer rootings from leaf edge rootability.
    n = len(side) - 2  # number of leaves on the side, called n in the paper
    for index, current in enumerate(reduced_side[1:-1]):
        i = index + 1
        # Find the leaf edge attached to current node
        leaf = get_leaf_adjacent_to_node(reduced_network, current)
        # Check if the leaf edge (current,leaf) is a root-edge in the reduced network
        rooting_at_edge = reduced_rootings.get((current, leaf))
        if not rooting_at_edge:
            rooting_at_edge = reduced_rootings.get((leaf, current))
        if not rooting_at_edge:
            continue

        # If the network is rootable at the leaf edge, infer rootings of the original network
        # First the leaf edges on the side
        for j in range(i, n - (ell - i) + 1):
            # Find the leaf edge at position j
            side_node = side[j]
            leaf_edge = (side_node, get_leaf_adjacent_to_node(network, side_node))
            # Now we extend to an orientation of the original network
            # be careful with inferring orientation: there can be a reticulation on the root side.
            # In that case, shift the reticulation to the relative position to the root as in the reduced network
            new_rooting = []
            for reticulation in rooting_at_edge:
                if reticulation in reduced_side[1:-1]:
                    # Here, we have found a reticulation node reticNode on the reduced side containing the root.
                    # We find its position on the reduced side, and give it the same relative position to the root on the original network
                    reduced_index = reduced_side.index(reticulation)
                    # Current position of the root is j
                    # In the reduced side, the root is at position i, and the retic at position reduced_index, a difference of reduced_index-i.
                    new_rooting += [side[j + reduced_index - i]]
                else:
                    # If the reticulation node is not on the root side, we can take this node as the reticulation (by placing the leaves back where they were removed)
                    new_rooting += [reticulation]
            rootings[leaf_edge] = tuple(new_rooting)
        # Now the internal edges of the side
        for j in range(i - 1, n - (ell - i) + 1):
            # Do something similar as for the leaf edges, to find the right reticulation node on the root side
            new_rooting = []
            for reticulation in rooting_at_edge:
                if reticulation in reduced_side[1:-1]:
                    reduced_index = reduced_side.index(reticulation)
                    # In the reduced side, the root is at position i, and the reticulation at position reduced_index, a difference of reduced_index-i.
                    # We correct for the fact that a new leaf is introduced between positions j and j+1 when we root at an internal edge:
                    relative_position = reduced_index - i
                    if reduced_index < i:
                        relative_position += 1
                    new_rooting += [side[j + relative_position]]
                else:
                    # If the reticulation node is not on the root side, we can take this node as the reticulation (by placing the leaves back where they were removed)
                    new_rooting += [reticulation]
            rootings[(side[j], side[j + 1])] = tuple(new_rooting)
    return rootings


def c_orientation_fpt_reticulation_number(network, ell, class_checker=is_network, use_cycle_basis=True):
    """Solves C-orientation for the given network and class
    This uses the FPT-time algorithm with the reticulation number as parameter
    (Algorithm 3) from the paper "Orienting Undirected Phylogenetic Networks".

    The class must be `ell`-chain reducible, blob-determined, and leaf-addable
    for the algorithm to be correct.

    Assumes the network has no pendant subtrees.

    Args:
        network (nx.Graph): The network that is to be oriented
        ell (int): The length to which chains are reduced.
        class_checker (function: nx.DiGraph -> Bool, optional):
            A function that determines whether a network is in a certain class
        use_cycle_basis (bool, optional):
            Whether to use the cycle basis method to restrict the possible reticulations

    Returns:
        dict(tuple(int)): A dict with all valid C-root-edges as keys, and a tuple of
        reticulation nodes of a valid C-orientation with this root edge of `network`.
    """

    if ell <= 2 and len(network.edges) == len(network.nodes):
        # In this case, the chain reduction would result in parallel edges, so keep the chain length 3 in this case.
        ell = 3
    reduced_network, sidesDict = reduce_chains(network, ell)
    if use_cycle_basis:
        reduced_rootings = c_orientation_exponential_cycle_basis(reduced_network, class_checker)
    else:
        reduced_rootings = c_orientation_exponential(reduced_network, class_checker)
    rootings = dict()
    for reduced_side, side in sidesDict.items():
        if len(side) <= ell + 2:
            rootings = {
                **rootings,
                **expand_reduced_rootings_for_short_side(
                    reduced_network, reduced_rootings, reduced_side
                ),
            }
        else:
            rootings = {
                **rootings,
                **expand_reduced_rootings_for_long_side(
                    network, reduced_network, reduced_rootings, reduced_side, side, ell
                ),
            }
    return rootings