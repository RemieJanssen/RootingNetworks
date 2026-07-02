def is_tree_based(network):
    """A function that returns True iff the given network is tree-based."""
    for node in network.nodes:
        if network.out_degree(node) == 1:
            # a network is tree-based if it has no W-fences.
            if is_endpoint_of_w_fence(network, node):
                return False
    return True


def is_endpoint_of_w_fence(network, reticulation):
    """Checks if the network has a W-fence with the given node at one of the endpoints of the fence"""
    for child in network.successors(reticulation):
        currentNode = child
    previousNode = reticulation
    currentlyAtTop = False
    done = False
    while not done:
        currentOutDegree = network.out_degree(currentNode)
        if currentOutDegree == 0:
            return False
        if currentOutDegree == 1:
            if currentlyAtTop:
                return True
            for node in network.predecessors(currentNode):
                if node != previousNode:
                    nextNode = node
        if currentOutDegree == 2:
            if not currentlyAtTop:
                return False
            for node in network.successors(currentNode):
                if node != previousNode:
                    nextNode = node
        previousNode, currentNode = currentNode, nextNode
        currentlyAtTop = not currentlyAtTop
    return False, "Error"
