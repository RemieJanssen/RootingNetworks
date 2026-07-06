import itertools
import math 
import networkx as nx

from phyloroot.class_checkers import is_network
from phyloroot.constrained_orientation import constrained_orientation_binary


def combinations_restricted(cycle_basis, internal_nodes, v_masks):
    """Generates all combinations of reticulations and restricts to the cycle basis.

    Args:
        cycle_basis (list(list(int))): A cycle basis of the network
        internal_nodes (list(int)): A list of internal nodes in the network
        v_masks (dict(int, int)): A dictionary mapping each internal node to a bitmask
            indicating which cycles it belongs to
    Yields:
        tuple(int): A combination of reticulations, one from each cycle in the cycle basis
    """
    full_mask = (1 << len(cycle_basis)) - 1

    for combination in itertools.combinations(internal_nodes, len(cycle_basis)):
        mask = 0
        # Check if the combination covers all cycles in the cycle basis
        # note: this does not guarantee that the combination has a 1-1 mapping to the cycles, 
        # but it is a necessary condition for such a mapping.
        for v in combination:
            mask |= v_masks[v]
        if mask == full_mask:
            yield combination

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


def root_at_edge_cycle_basis_combinations(network, root_edge, cycle_basis, internal_nodes, v_masks, class_checker=is_network):
    """Subroutine for c_orientation_exponential_cycle_basis.
    Checks C-orientation for a given root-edge.

    Args:
        network (nx.Graph): The network that is to be oriented
        root_edge (tuple(int, int)): An edge of `network`
        cycle_basis (list(list(int))): A cycle basis of `network`
        internal_nodes (list(int)): A list of internal nodes in the network
        v_masks (dict(int, int)): A dictionary mapping each internal node to a bitmask
            indicating which cycles it belongs to
        class_checker (function: nx.DiGraph -> Bool, optional):
            A function that determines whether a network is in a certain class

    Returns:
        tuple(int) or bool: The tuple of reticulation nodes of a valid
           C-orientation with root edge `root_edge` of `network` if it exists,
           and False otherwise
    """
    for reticulations_set in combinations_restricted(cycle_basis, internal_nodes, v_masks):
        result = constrained_orientation_binary(network, root_edge, reticulations_set)
        if result and class_checker(result):
            return list(reticulations_set)
    return False

def root_at_edge_cycle_basis_product(network, root_edge, cycle_basis, class_checker=is_network):
    """Subroutine for c_orientation_exponential_cycle_basis.
    Checks C-orientation for a given root-edge.

    Args:
        network (nx.Graph): The network that is to be oriented
        root_edge (tuple(int, int)): An edge of `network`
        cycle_basis (list(list(int))): A cycle basis of `network`
        class_checker (function: nx.DiGraph -> Bool, optional):
            A function that determines whether a network is in a certain class

    Returns:
        tuple(int) or bool: The tuple of reticulation nodes of a valid
           C-orientation with root edge `root_edge` of `network` if it exists,
           and False otherwise
    """
    for reticulations in itertools.product(*cycle_basis):
        reticulations_set = set(reticulations)
        if len(reticulations_set) != len(cycle_basis):
            continue
        result = constrained_orientation_binary(network, root_edge, reticulations_set)
        if result and class_checker(result):
            return list(reticulations_set)
    return False



def c_orientation_exponential(network, class_checker=is_network):
    """Solves C-orientation for the given network and class
    This uses the exponential-time algorithm (Algorithm 2) from
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


def c_orientation_exponential_cycle_basis(network, class_checker=is_network):
    """Solves C-orientation for the given network and class
    This uses the exponential-time algorithm (Algorithm 2) from
    the paper "Orienting Undirected Phylogenetic Networks".

    Uses Urata et al.'s cycle basis method to restrict the possible reticulations to the nodes in the cycle basis.
    If the product of the cycle lengths is smaller than the number of combinations of reticulations,
    we iterate over the product of the cycle lengths instead of all combinations of reticulations.
    Otherwise, we iterate over all combinations of reticulations as in the original algorithm, 
    and we make sure at least one node from each cycle is selected.
    Note that this does not guarantee a 1-1 mapping between cycles and reticulations, but it is a necessary condition for such a mapping.

    Args:
        network (nx.Graph): The network that is to be oriented
        class_checker (function: nx.DiGraph -> Bool, optional):
            A function that determines whether a network is in a certain class

    Returns:
        dict(tuple(int): A dict with all valid C-root-edges as keys, and a tuple of
        reticulation nodes of a valid C-orientation with this root edge of `network`.
    """
    rootings = dict()
    cycle_basis = sorted(nx.minimum_cycle_basis(network), key=len)
    # Get all internal nodes (nodes with degree > 1)
    internal_nodes = [v for v in network.nodes if network.degree(v) > 1]
    if math.comb(len(internal_nodes), len(cycle_basis)) < math.prod(len(c) for c in cycle_basis):
        # if the number of combinations of reticulations is smaller than the product of the cycle lengths, 
        # we iterate over all combinations of reticulations as in the original algorithm, and we make sure at least one node from each cycle is selected.
        
        # Create a mask for each internal node indicating which cycles it belongs to
        v_masks = dict()
        for v in internal_nodes:
            v_mask = 0
            for i, cycle in enumerate(cycle_basis):
                if v in cycle:
                    v_mask |= 1 << i
            v_masks[v] = v_mask

        for e in network.edges:
            result = root_at_edge_cycle_basis_combinations(network, e, cycle_basis, internal_nodes, v_masks, class_checker)
            if result:
                rootings[e] = result
    else:
        # if the product of the cycle lengths is smaller than the number of combinations of reticulations, 
        # we iterate over the product of the cycle lengths instead of all combinations of reticulations.
        for e in network.edges:
            result = root_at_edge_cycle_basis_product(network, e, cycle_basis, class_checker)
            if result:
                rootings[e] = result
    return rootings


