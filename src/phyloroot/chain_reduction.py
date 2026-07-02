import networkx as nx
from copy import deepcopy


def get_generator_with_degree_two_nodes(network):
    """Returns the generator of the undirected network including degree-2 nodes.
    Assumes the network has no pendant subtrees"""
    generator = deepcopy(network)
    todoNodes = list(generator.nodes)
    while todoNodes:
        current = todoNodes.pop()
        if generator.degree(current) == 1:
            generator.remove_node(current)
    return generator


def get_side(generator, generator_node, side_node):
    """Finds the side of the generator with degree-2 nodes that starts with a given generator node and an adjacent node"""
    previous = generator_node
    current = side_node
    side = (
        previous,
        current,
    )
    while generator.degree(current) == 2:
        for nb in generator.neighbors(current):
            if nb != previous:
                next = nb
        previous, current = current, next
        side += (current,)
    return side


def get_all_sides(network):
    """Returns all sides of the undirected network
    Assumes the network no pendant subtrees"""
    generator = get_generator_with_degree_two_nodes(network)
    todoNodes = list(generator.nodes)
    sides = set()
    # Do something separate for cycles: in those cases, the sides do not end in generator nodes
    number_of_degree_3_nodes = 0
    while todoNodes:
        current = todoNodes.pop()
        if generator.degree(current) == 3:
            number_of_degree_3_nodes += 1
            # Now look for sides around the generator node
            for nb in network.neighbors(current):
                side = get_side(generator, current, nb)
                if not tuple(reversed(side)) in sides:
                    sides.add(side)
    if number_of_degree_3_nodes == 0:
        cycle = nx.cycle_basis(generator)[0]
        side = tuple([cycle[-1]] + cycle + [cycle[0]])
        sides.add(side)
    return sides


def reduce_chains(network, ell):
    """Reduces chains in the undirected network on all its sides to length at most `ell`
    When reducing, it removes internal nodes of the chain.
    Assumes the network has no pendant subtrees

    Returns the reduced network, and a dictionary of sides which maps the reduced side to the original side.
    """
    sides = get_all_sides(network)
    sidesDict = dict()
    reducedNetwork = nx.Graph()
    for side in sides:
        if len(side) > ell + 2:
            # Keep the first l-1 and the last leaf on the chain
            newSide = side[:ell] + side[-2:]
        else:
            newSide = side[:]
        sidesDict[newSide] = side
        # Now add the new side to the new network
        previousNode = newSide[0]
        for node in newSide[1:]:
            reducedNetwork.add_edge(previousNode, node)
            previousNode = node
    # Add the leaves to the reduced network
    for node in list(reducedNetwork.nodes):
        if reducedNetwork.degree(node) == 2:
            for nb in network.neighbors(node):
                if network.degree(nb) == 1:
                    leaf = nb
            reducedNetwork.add_edge(node, leaf)
    return reducedNetwork, sidesDict
