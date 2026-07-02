def is_tree_child(network):
    """A function that checks whether a network is in the class of tree-child networks

    returns True if the given network is tree-child and False otherwise
    """
    for node in network.nodes:
        node_is_tree_child = False
        # A leaf is a tree-child node
        if network.out_degree(node) == 0:
            node_is_tree_child = True
        for child in network.successors(node):
            # The node is tree-child if at least one of its children is a leaf or a tree-node
            if network.out_degree(child) in [0, 2]:
                node_is_tree_child = True
        if not node_is_tree_child:
            return False
    return True
