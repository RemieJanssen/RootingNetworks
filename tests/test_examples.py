import pytest

from expected_outcomes import outcomes
from phyloroot.main_cli import get_class_checker_and_chain_length, read_network_file
from phyloroot.rooting import LevelStuff

tests = [
    {
        "testname": f"{filename} {class_name}",
        "network_path": f"./examples/{filename}",
        "class_abbr": class_results["abbreviation"],
        "root_edges": class_results.get("root_edges"),
        "non_root_edges": class_results.get("non_root_edges")
    }
    for filename, file_tests in outcomes.items()
    for class_name, class_results in file_tests.items()
]

@pytest.mark.parametrize(
    "network_path, class_abbr, root_edges, non_root_edges",
    [
        (test["network_path"], test["class_abbr"], test["root_edges"], test["non_root_edges"])
        for test in tests
    ],
    ids=[test["testname"] for test in tests],
)
def test_level_class(network_path, class_abbr, root_edges, non_root_edges):
    network = read_network_file(network_path)
    ClassChecker, length = get_class_checker_and_chain_length(class_abbr)
    orientations = LevelStuff(network, length, ClassChecker)
    if not orientations:
        orientations = []
    found_root_edges = {tuple(sorted(edge)) for edge in orientations}
    if root_edges is not None:
        root_edges = {tuple(sorted(edge)) for edge in root_edges}
    else:
        non_root_edges = {tuple(sorted(edge)) for edge in non_root_edges}
        root_edges = {tuple(sorted(edge)) for edge in network.edges}.difference(non_root_edges)
    assert found_root_edges == root_edges

