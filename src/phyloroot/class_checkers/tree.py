import networkx as nx

def is_tree(digraph):
    """Checks if a directed graph is a rooted tree. If so, it returns the root node"""
    if len(digraph.edges)!=len(digraph.nodes)-1 or not nx.is_weakly_connected(digraph):
        return False
    root_found = False
    root = False
    for node in digraph:
        if digraph.in_degree(node)==0:
            if root_found:
                return False
            root_found = True
            root = node
        if digraph.in_degree(node)>1:
            return False
    return root
