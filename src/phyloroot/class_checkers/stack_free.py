def is_stack_free(network):
    """A function that returns True iff the given network is stack-free."""
    for node in network.nodes:
        sf_node = True
        if network.out_degree(node) == 1:
            for child in network.successors(node):
                outdegChild = network.out_degree(child)
                if outdegChild == 1:
                    sf_node = False
        if not sf_node:
            return False
    return True
