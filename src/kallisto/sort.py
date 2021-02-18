# src/kallisto/sort.py


from collections import defaultdict

import click

from kallisto.data import chemical_symbols
import kallisto.reader.strucreader as ksr
from kallisto.units import Bohr


class Graph:
    """Define a molecular graph."""

    def __init__(self, inp: str, out: click.File):

        self.inp = inp
        self.graph = defaultdict(list)  # type: ignore
        self.molecule = ksr.constructMolecule(geometry=self.inp, out=out)
        self.nat = self.molecule.get_number_of_atoms()
        self.at = self.molecule.get_atomic_numbers()
        self.coordinates = self.molecule.get_positions()

    # function for adding an edge to graph
    def addEdge(self, u, v):
        self.graph[u].append(v)

    # method to print a breadth first search
    def BFS(self, s):

        # mark all vertices as not visited
        visited = [False] * (len(self.graph))

        # create a queue for BFS
        queue = []

        # mark source node as visited and enqueue it
        queue.append(s)
        visited[s] = True

        print("{:5}".format(self.nat), end=" \n")
        print("Created with kallisto")
        while queue:

            # dequeue a vertex from queue
            # and print it
            s = queue.pop(0)
            print(
                "{:3} {:9.4f} {:9.4f} {:9.4f}".format(
                    chemical_symbols[self.at[s]],
                    self.coordinates[s][0] * Bohr,
                    self.coordinates[s][1] * Bohr,
                    self.coordinates[s][2] * Bohr,
                ),
                end=" \n",
            )

            # get adjacent vertices of dequeued vertex s
            # If an adjacent has not been visited, then mark
            # it as visited and enqueue it
            for i in self.graph[s]:
                if visited[i] is False:
                    queue.append(i)
                    visited[i] = True
