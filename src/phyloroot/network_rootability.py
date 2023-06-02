from rooting import *
import networkx as nx
import os
import sys
import csv
import argparse


def read_network_file(filename):
    # Parse the network in <filename>
    edges = []
    with open(filename, "rt") as f:
        reader = csv.reader(f, delimiter=",", quotechar="|")
        for row in reader:
            if len(row) == 2:
                edges.append((int(row[0]), int(row[1])))

    # Build the unrooted network to be oriented
    network = nx.Graph()
    network.add_edges_from(edges)
    return network


def get_class_checker_and_chain_length(nwClass):
    """Set the class parameters corresponding to nwClass abbreviation"""
    ClassChecker = ClassAllNetworks
    length = 1
    if nwClass == "TC":
        ClassChecker = ClassTreeChild
        length = 3
    elif nwClass == "SF":
        ClassChecker = ClassStackFree
        length = 3
    elif nwClass == "O":
        ClassChecker = ClassOrchard
        length = 3
    elif nwClass == "TB":
        ClassChecker = ClassTreeBased
        length = 2
    return ClassChecker, length


def set_output(network, orientations, simple):
    output = ""
    if orientations:
        rootableEdges = orientations.keys()
        # print("rootable at: ", rootableEdges)
        # print("not at: ", set(network.edges)-set(rootableEdges))
        if not simple:
            output += "The root edges, reticulations and orientations are:"
        else:
            output += "The orientations of the network consist of the following root edges and reticulations:\r\n"
        for rootEdge, reticulations in orientations.items():
            # compute the actual orientation, instead of only the retculations.
            dinetwork = OrientationAlgorithmBinary(network, rootEdge, reticulations)
            if not simple:
                output += (
                    "\r\n\r\nroot edge:\r\n   "
                    + str(rootEdge)
                    + "\r\nreticulations:\r\n   "
                    + str(reticulations)
                    + "\r\norientation:\r\n   "
                    + str(dinetwork.edges)
                )
            else:
                output += "\r\n" + str(rootEdge) + " " + str(reticulations)
        if not simple:
            output += "\r\n\r\nThe network cannot be rooted at the edges:"
            for non_root in set(network.edges) - set(rootableEdges):
                output += "\r\n" + str(non_root)
    else:
        output = "There is no orientation of the given network in the desired class."
    return output


def cmd_parser():
    parser = argparse.ArgumentParser(
        description="finds the orientations of an undirected phylogenetic network that belong to a given class of directed networks."
    )
    parser.add_argument(
        "-f",
        "--file",
        help="input file with an undirected phylogenetic network as a list of edges. One edge per line, vertices are positive integers separated by a comma",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output file",
    )
    parser.add_argument(
        "-c",
        "--classname",
        help="A class of networks, choose from: TC (tree-child) TB (tree-based), SF (stack-free), and O (Orchard). If not set, this will return all orientations.",
    )
    parser.add_argument(
        "-s",
        "--simple",
        help="Use this option if you want output only consisting of the root-edge and the reticulation nodes for each orientation.",
        action="store_true",
    )
    return parser.parse_args()


def main():
    cmd_args = cmd_parser()
    network = read_network_file(cmd_args.file)
    ClassChecker, length = get_class_checker_and_chain_length(cmd_args.classname)
    orientations = LevelStuff(network, length, ClassChecker)
    output = set_output(network, orientations, cmd_args.simple)

    # write output to file or to the shell
    if cmd_args.output:
        with open(cmd_args.output, "w+") as f:
            f.write("The input network was:\r\n")
            f.write("  " + str(network.edges) + "\r\n")
            f.write("The class is:\r\n")
            f.write("  " + cmd_args.classname + "\r\n")
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
