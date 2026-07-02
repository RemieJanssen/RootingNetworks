import networkx as nx

def constrained_orientation_binary(network,root_edge,reticulations):
    """Solves the constrained orientation problem for binary networks.
    Uses the Orientation Algorithm (Algorithm 1) from the paper
    Orienting Undirected Phylogenetic Networks.

    Args:
        network (nx.Graph): The network that is to be oriented.
        root_edge (int): The (ID of the) designated root edge in `network`.
        reticulations (_type_): The (IDs of the) designated reticulation nodes of `network`.

    Returns:
        nx.DiGraph or Bool: The oriented network if it exists, and False otherwise
    """
    number_of_edges = len(network.edges)
    network.remove_edges_from([root_edge])
    if len(network.nodes)+len(reticulations)!=number_of_edges+1:
        network.add_edges_from([root_edge])
        return False
    diNetwork=nx.DiGraph()
    diNetwork.add_nodes_from(network.nodes)
    diNetwork.add_node(-1)
    diNetwork.add_edges_from([(-1,root_edge[0]),(-1,root_edge[1])])
    readyNodes = (network.nodes - reticulations) & set(root_edge)
    numberOrientedEdges = 2
    while len(readyNodes)>0: #subdividing the root edge gives one additional edge
        orientingNode = readyNodes.pop()
        orientingNodeInDeg = diNetwork.in_degree(orientingNode)
        if orientingNodeInDeg>2 or (orientingNodeInDeg==2 and orientingNode not in reticulations):
            network.add_edges_from([root_edge])
            return False
        children = set(network.neighbors(orientingNode))-set(diNetwork.predecessors(orientingNode))
        newArcs = []
        for child in children:
            newArcs +=[(orientingNode,child)]
            if diNetwork.in_degree(child)==1 or child not in reticulations:
                readyNodes.add(child)
            numberOrientedEdges+=1
        diNetwork.add_edges_from(newArcs)

    network.add_edges_from([root_edge])
    if numberOrientedEdges < number_of_edges+1:
        return False
    return diNetwork
