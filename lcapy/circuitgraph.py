"""
This module provides a class to represent circuits as graphs.
This is primarily for loop analysis but is also used for nodal analysis.

Copyright 2019--2020 Michael Hayes, UCECE

"""

from matplotlib.pyplot import subplots, savefig
import networkx as nx


# MultiGraph handles parallel edges.

class CircuitGraph(nx.MultiGraph):

    def __init__(self, cct):

        super(CircuitGraph, self).__init__()
        self.cct = cct

        self.add_nodes_from(cct.node_list)

        node_map = cct.node_map

        for name in cct.branch_list:
            elt = cct.elements[name]
            if len(elt.nodenames) < 2:
                continue
            self.add_edge(node_map[elt.nodenames[0]],
                          node_map[elt.nodenames[1]], name=name)

        self.node_map = node_map

    def connected(self, node):
        """List of components connected to specified node."""

        for node1, edges in self.node_edges(node).items():
            for key, edge in edges.items():
                name = edge['name']
                elt = self.cct.elements[name]
                yield elt
            
    def all_loops(self):
        
        DG = nx.MultiDiGraph(self)
        cycles = list(nx.simple_cycles(DG))

        loops = []
        for cycle in cycles:
            if len(cycle) > 2:
                cycle = sorted(cycle)
                if cycle not in loops:
                    loops.append(cycle)        
        return loops

    def chordless_loops(self):

        loops = self.all_loops()
        sets = [set(loop) for loop in loops]

        rejects = []
        for i in range(len(sets)):
            for j in range(i + 1, len(sets)):
                if sets[i].union(sets[j]) == sets[i]:
                    rejects.append(i)
                elif sets[i].union(sets[j]) == sets[j]:
                    rejects.append(j)

        cloops = []
        for i, loop in enumerate(loops):
            if i not in rejects:
                cloops.append(loop)

        return cloops

    def loops(self):
        return self.chordless_loops()
    
    def draw(self, filename=None):
        """Use matplotlib to draw circuit graph."""

        fig, ax = subplots(1)

        G = self
        pos = nx.spring_layout(G)
        
        labels = dict(zip(G.nodes(), G.nodes()))
        nx.draw_networkx(G, pos, ax, labels=labels)
        
        edge_labels = dict([((u, v), d['name'])
                            for u, v, d in G.edges(data=True)])
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

        if filename is not None:
            savefig(filename, bbox_inches='tight')
        

    def node_edges(self, node):

        return self[node]
        
# neighbors(node) gives neighbouring nodes
