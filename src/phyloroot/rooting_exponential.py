import itertools

import networkx as nx

from phyloroot.class_checkers import is_network
from phyloroot.constrained_orientation import constrained_orientation_binary


def _combinations_from_product(cycle_basis, chosen=None, i=0):
    """Generates all combinations of reticulations from the cycle basis.

    Args:
        cycle_basis (list(list(int))): A cycle basis of the network
    Yields:
        tuple(int): A combination of reticulations, one from each cycle in the cycle basis    
    """
    chosen = chosen or set()

    if i == len(cycle_basis):
        yield chosen
        return

    for v in cycle_basis[i]:
        if v in chosen:
            continue 
        chosen.add(v)
        yield from _combinations_from_product(cycle_basis, chosen, i + 1)
        chosen.remove(v)


def combinations_from_product(cycle_basis):
    """Generates all combinations of reticulations from the cycle basis.

    Args:
        cycle_basis (list(list(int))): A cycle basis of the network
    Yields:
        tuple(int): A combination of reticulations, one from each cycle in the cycle basis    
    """
    for comb in _combinations_from_product(cycle_basis, None, 0):
        yield comb


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


def root_at_edge_cycle_basis(network, root_edge, cycle_basis, class_checker=is_network):
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
    for reticulations_set in combinations_from_product(cycle_basis):
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
    for e in network.edges:
        result = root_at_edge_cycle_basis(network, e, cycle_basis, class_checker)
        if result:
            rootings[e] = result
    return rootings


