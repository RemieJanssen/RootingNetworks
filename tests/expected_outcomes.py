from phyloroot.rooting import is_network, is_tree_child, is_stack_free, is_orchard, is_tree_based

outcomes = {
    "Example1-Components.txt": {
        "tree-child": {"abbreviation": "TC", "class_checker": is_tree_child, "length": 3, "root_edges": [(1,2), (1,3), (1,4), (3,4), (3,5), (5,6), (4,8), (5,7), (7,8), (8,9)]},
        "stack-free": {"abbreviation": "SF", "class_checker": is_stack_free, "length": 3, "non_root_edges": []},
        "orchard": {"abbreviation": "O", "class_checker": is_orchard, "length": 3, "non_root_edges": [(7,10), (10,11), (10,12), (11,12), (11,13), (12,14)]},
        "TB": {"abbreviation": "TB", "class_checker": is_tree_based, "length": 2, "non_root_edges": []},
        "all": {"abbreviation": None, "class_checker": is_network, "length": 2, "non_root_edges": []},
    },
    "Example2-PetersenOneLeaf.txt": {
        "tree-child": {"abbreviation": "TC", "class_checker": is_tree_child, "length": 3, "root_edges": []},
        "stack-free": {"abbreviation": "SF", "class_checker": is_stack_free, "length": 3, "root_edges": []},
        "orchard": {"abbreviation": "O", "class_checker": is_orchard, "length": 3, "root_edges": []},
        "TB": {"abbreviation": "TB", "class_checker": is_tree_based, "length": 2, "non_root_edges": [(100,101), (1,100), (0,100), (0,4), (0,5), (1,2), (1,6)]},
        "all": {"abbreviation": None, "class_checker": is_network, "length": 2, "non_root_edges": [(100,101)]},
    },
    # "Example2.1-PetersenTwoLeaves.txt": {},
    "Example3-StackFree.txt": {
        "tree-child": {"abbreviation": "TC", "class_checker": is_tree_child, "length": 3, "root_edges": []},
        "stack-free": {"abbreviation": "SF", "class_checker": is_stack_free, "length": 3, "root_edges": []},
        # "orchard": {"abbreviation": "O", "class_checker": ClassOrchard, "length": 3, "root_edges": []},
        # "TB": {"abbreviation": "TB", "class_checker": ClassTreeBased, "length": 2, "root_edges": []},
        "all": {"abbreviation": None, "class_checker": is_network, "length": 2, "non_root_edges": []},
    },
    # "Example4-BigExampleFromPaper.txt": {},
    # "Example4.1-BigExampleFromPaperFewerLeaves.txt": {},
    # "Example4.2-BigExampleFromPaperReduced.txt": {},
    "Example5-Subtrees.txt": {
        "tree-child": {"abbreviation": "TC", "class_checker": is_tree_child, "length": 3, "root_edges": [(0,6), (0,7), (6,7), (6,8), (7,14)]},
        "stack-free": {"abbreviation": "SF", "class_checker": is_stack_free, "length": 3, "non_root_edges": []},
        "orchard": {"abbreviation": "O", "class_checker": is_orchard, "length": 3, "non_root_edges": [(8,9), (9,10), (9,11), (10,12), (10,13), (8,14), (14,15)]},
        "TB": {"abbreviation": "TB", "class_checker": is_tree_based, "length": 2, "non_root_edges": []},
        "all": {"abbreviation": None, "class_checker": is_network, "length": 2, "non_root_edges": []},
    }
}

