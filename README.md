# Phyloroot

Phyloroot is a python module with tools for orienting undirected phylogenetic networks as oriented networks in several well-known classes.

## Install

Install as pypi package phyloroot:
```
    pip install phyloroot
```

## Usage

You can load the package methods with `import phyloroot` in python, or use the command line interface.
For example, to get all orchard orientations of the example network in `./examples/Example1-Components.txt`, execute the following command:
```
phyloroot -f ./examples/Example1-Components.txt -c O
```
For more information, use the help function:
```
phyloroot -h
```