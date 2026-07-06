import networkx as nx

from phyloroot.class_checkers import is_network
from phyloroot.rooting_fpt_reticulation_number import c_orientation_fpt_reticulation_number

def is_tree(digraph):
    """Checks if a directed graph is a rooted tree. If so, it returns the root node"""
    if len(digraph.edges) != len(digraph.nodes) - 1 or not nx.is_weakly_connected(
        digraph
    ):
        return False
    root_found = False
    root = False
    for node in digraph:
        if digraph.in_degree(node) == 0:
            if root_found:
                return False
            root_found = True
            root = node
        if digraph.in_degree(node) > 1:
            return False
    return root


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


def construct_full_orientation(blobs_and_leaf_edges, blob_orientations, root_edge):
    """Constructs a full orientation of the network from the orientations of the blobs and a given root edge.

    Starts at the blob containing the root edge, and moves away from the root to orient all other blobs.

    Args:
        blobs_and_leaf_edges (list(tuple(nx.Graph, list(tuple(int))))): a list of tuples, where each tuple contains a blob (as a subgraph) and a list of leaf edges (as tuples of nodes with the leaf as the second element of the pair) for that blob.
        blob_orientations (list(dict(tuple(int), list(int)))): a list of dictionaries, where each dictionary contains the orientations of a blob (as keys) and the corresponding reticulation nodes (as values).
        root_edge (tuple(int)): The root edge of the network.
    
    Returns:
        list(int) or bool: A list of reticulation nodes of a valid C-orientation with the given root edge of the network if it exists, and False otherwise.
    """
    # Find a blob containing the root edge
    for i, (blob, blob_leaf_edges) in enumerate(blobs_and_leaf_edges):
        if blob.has_edge(*root_edge):
            break

    reticulations = []
    edges_to_continue_at = False
    # check if the blob containing the potential root edge can be rooted at this edge, and if so, continue to root the other blobs
    if len(blob) == 2:
        edges_to_continue_at = set([root_edge, (root_edge[1], root_edge[0])])
    elif root_edge in blob_orientations[i]:
        reticulations += blob_orientations[i][root_edge]
        edges_to_continue_at = blob_leaf_edges
    elif (root_edge[1], root_edge[0]) in blob_orientations[i]:
        reticulations += blob_orientations[i][(root_edge[1], root_edge[0])]
        edges_to_continue_at = blob_leaf_edges
    else:
        # If the blob containing the potential root edge cannot be rooted at this edge, continue to the next edge
        return False

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
                other_node = False
                for v in blob:
                    if v != edge[1]:
                        other_node = v
                edges_to_continue_at.add((edge[1], other_node))
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
    return reticulations

        
def c_orientation_fpt_level(network, ell, ClassChecker=is_network, use_cycle_basis=True):
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
        use_cycle_basis (bool, optional):
            Whether to use the cycle basis method to restrict the possible reticulations

    Returns:
        dict(tuple(int): A dict with all valid C-root-edges as keys, and a tuple of
        reticulation nodes of a valid C-orientation with this root edge of `network`.
    """

    blobs_and_leaf_edges = get_blobs_and_leaf_edges(network)
    # Prepare the partially directed network that we will condense into T_CN
    partially_oriented_network = network.to_directed()
    # Empty list to store all orientations for all blobs
    blob_orientations = []
    # For each blob, compute the orientations, and orient the leaf-edges of the blobs in
    # partially_oriented_network only in the directions that are allowed by the orientations of the blobs.
    for blob, blob_leaf_edges in blobs_and_leaf_edges:
        if len(blob) <= 2:
            # For trivial biconnected components, add a trivial list to the blob orientations, 
            # so that indices still match between blob_orientations, and blobs
            blob_orientations += [[]]
            continue
        # Find all rootings of the blob
        blob_rootings = c_orientation_fpt_reticulation_number(
            blob, ell, ClassChecker, use_cycle_basis=use_cycle_basis
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

    # Create T_CN by condensing partially_oriented_network
    T_CN = nx.condensation(partially_oriented_network)
    # Find the root of T_CN if it exists; if it does not, the network is not C-orientable
    root_component = is_tree(T_CN)
    if not type(root_component) == int:
        return False
    root_component_nodes = T_CN.nodes(data=True)[root_component]["members"]

    # Go through all edges to find all root-edges
    rootings = dict()
    for root_edge in network.edges:
        # Check if the edge is in the root_component
        if not (root_edge[0] in root_component_nodes and root_edge[1] in root_component_nodes):
            continue
        extended_rooting = construct_full_orientation(
            blobs_and_leaf_edges, blob_orientations, root_edge
        )
        if not extended_rooting:
            continue
        rootings[root_edge] = extended_rooting
    return rootings
