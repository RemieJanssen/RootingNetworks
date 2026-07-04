import networkx as nx
import itertools

from phyloroot.class_checkers import (
    is_network,
    is_tree,
)
from phyloroot.constrained_orientation import constrained_orientation_binary
from phyloroot.chain_reduction import reduce_chains


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


def root_at_edge(network, root_edge, class_checker=is_network):
    """Subroutine for c_orientation_exponential.
    Checks C-orientation for a given root-edge.

    Args:
        network (nx.Graph): The network that is to be oriented
        root_edge (tupe/list(int, int)): An edge of `network`
        class_checker (function: nx.DiGraph -> Bool, optional):
            A function that determines whether a network is in a certain class

    Returns:
        tuple(int) or bool: The tuple of reticulation nodes of a valid
           C-orientation with root edge `root_edge` of `network` if it exists,
           and False otherwise
    """
    noOfReticulations = len(network.edges) - len(network.nodes) + 1
    nonLeafNodes = []
    for node in network.nodes:
        if network.degree(node) > 1:
            nonLeafNodes += [node]
    for reticulations in itertools.combinations(nonLeafNodes, noOfReticulations):
        result = constrained_orientation_binary(network, root_edge, set(reticulations))
        if result and class_checker(result):
            return reticulations
    return False


def c_orientation_exponential(network, class_checker=is_network):
    """Solves C-orientation for the given network and class
    This uses the txponential-time algorithm (Algorithm 2) from
    the paper "Orienting Undirected Phylogenetic Networks".

    Args:
        network (nx.Graph): The network that is to be oriented
        class_checker (function: nx.DiGraph -> Bool, optional):
            A function that determines whether a network is in a certain class

    Returns:
        dict(tuple(int): A dict with all valid C-root-edges as keys, and a tuple of
        reticulation nodes of a valid C-orientation with this root edge of `network`.
    """
    rootings = dict()
    for e in network.edges:
        result = root_at_edge(network, e, class_checker)
        if result:
            rootings[e] = result
    return rootings


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


def c_orientation_fpt_reticulation_number(network, ell, class_checker=is_network):
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

    Returns:
        dict(tuple(int)): A dict with all valid C-root-edges as keys, and a tuple of
        reticulation nodes of a valid C-orientation with this root edge of `network`.
    """

    if ell <= 2 and len(network.edges) == len(network.nodes):
        # In this case, the chain reduction would result in parallel edges, so keep the chain length 3 in this case.
        ell = 3
    reduced_network, sidesDict = reduce_chains(network, ell)
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


def get_blobs_and_leaf_edges(network):
    """Returns a list of blobs and a list of leaf edges for the given network.

    Args:
        network (nx.Graph): The network
    Returns:
        tuple(list(nx.Graph), list(tuple(int))): a list of tuples, where each tuple contains a blob (as a subgraph) 
            and a list of leaf edges (as tuples of nodes with the leaf as the second element of the pair) for that blob.
    """
    blobs_and_leaf_edges = []
    for c in nx.biconnected_components(network):
        blob = network.subgraph(c).copy()
        leaf_edges = set()
        for node in blob.nodes:
            if blob.degree(node) == 2:
                for nb in network.neighbors(node):
                    if nb not in blob:
                        leaf_edges.add((node, nb))
        blob.add_edges_from(leaf_edges)
        blobs_and_leaf_edges.append((blob, leaf_edges))
    return blobs_and_leaf_edges


def c_orientation_fpt_level(network, ell, ClassChecker=is_network):
    """Solves C-orientation for the given network and class
    This uses the FPT-time algorithm with the level as parameter
    (Algorithm 4) from the paper "Orienting Undirected Phylogenetic Networks".

    The class must be `ell`-chain reducible, blob-determined, and leaf-addable
    for the algorithm to be correct.

    Args:
        network (nx.Graph): The network that is to be oriented
        ell (int): The length to which chains are reduced.
        class_checker (function: nx.DiGraph -> Bool, optional):
            A function that determines whether a network is in a certain class

    Returns:
        dict(tuple(int): A dict with all valid C-root-edges as keys, and a tuple of
        reticulation nodes of a valid C-orientation with this root edge of `network`.
    """

    blobs_and_leaf_edges = get_blobs_and_leaf_edges(network)
    # Prepare the partially directed network that we will condense into T_CN
    partially_oriented_network = network.to_directed()
    # Empty list to store all orientations for all blobs
    blob_orientations = []
    # For each blob, compute the orientations
    for blob, blob_leaf_edges in blobs_and_leaf_edges:
        if len(blob) <= 2:
            # For trivial biconnected components, add a trivial list to the blob orientations, 
            # so that indices still match between blobOrientations, and blobs
            blob_orientations += [[]]
            continue
        # Find all rootings of the blob
        blob_rootings = c_orientation_fpt_reticulation_number(
            blob, ell, ClassChecker
        )
        blob_orientations += [blob_rootings]
        if not blob_rootings:
            # If there is no orientation for this blob, then there is no orientation for the network.
            return False
        # Partially orient at the leaves, according to where the blob can be rooted
        for leaf_edge in blob_leaf_edges:
            if (
                leaf_edge in blob_rootings
                or (leaf_edge[1], leaf_edge[0]) in blob_rootings
            ):
                continue
            # If now both arcs are gone, there is no rooting of the original network, so we may return False
            if not partially_oriented_network.has_edge(*leaf_edge):
                return False
            partially_oriented_network.remove_edge(leaf_edge[1], leaf_edge[0])
    # Create T_CN by condensing partiallyOrientedNetwork
    T_CN = nx.condensation(partially_oriented_network)
    # Find the root of T_CN if it exists; if it does not, the network is not C-orientable
    rootComponent = is_tree(T_CN)
    if not type(rootComponent) == int:
        return False

    # Go through all edges to find all orientations
    rootings = dict()
    root_component_nodes = T_CN.nodes(data=True)[rootComponent]["members"]
    for rootEdge in network.edges:
        reticulations = []
        edges_to_continue_at = False
        # Check if the edge is in the rootComponent
        if not (rootEdge[0] in root_component_nodes and rootEdge[1] in root_component_nodes):
            continue
        # Find a blob containing the root edge
        for i, (blob, blob_leaf_edges) in enumerate(blobs_and_leaf_edges):
            if blob.has_edge(*rootEdge):
                break

        # check if the blob containing the potential root edge can be rooted at this edge, and if so, continue to root the other blobs
        if len(blob) == 2:
            edges_to_continue_at = set([rootEdge, (rootEdge[1], rootEdge[0])])
        elif rootEdge in blob_orientations[i]:
            reticulations += blob_orientations[i][rootEdge]
            edges_to_continue_at = blob_leaf_edges
        elif (rootEdge[1], rootEdge[0]) in blob_orientations[i]:
            reticulations += blob_orientations[i][(rootEdge[1], rootEdge[0])]
            edges_to_continue_at = blob_leaf_edges
        else:
            # If the blob containing the potential root edge cannot be rooted at this edge, continue to the next edge
            continue

        # Continue to root all other blobs, by moving away from the blob with the root.
        # edges_to_continue_at keeps a list of edges along which we still have to move away from the root
        while edges_to_continue_at:
            edge = edges_to_continue_at.pop()
            # find the blob this edge points to
            for i, (blob, blob_leaf_edges) in enumerate(blobs_and_leaf_edges):
                if edge[1] not in blob:
                    continue
                # Continue at trivial biconnected components with an endpoint edge[1] (but not edge[0], to prevent cycling in the algorithm)
                if len(blob) == 2 and edge[0] not in blob:
                    otherNode = False
                    for v in blob:
                        if v != edge[1]:
                            otherNode = v
                    edges_to_continue_at.add((edge[1], otherNode))
                    break
                # Continue at blobs that contain edge[1] in the interior (so the degree of edge[1] in the blob is not 1)
                if len(blob) > 2 and blob.degree(edge[1]) != 1:
                    edges_to_continue_at |= blob_leaf_edges - set(
                        [(edge[1], edge[0])]
                    )
                    if edge in blob_orientations[i]:
                        reticulations += blob_orientations[i][edge]
                    else:
                        reticulations += blob_orientations[i][
                            (edge[1], edge[0])
                        ]
                    break
        rootings[rootEdge] = reticulations
    return rootings
