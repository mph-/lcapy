"""
This module provides a class to represent circuits as graphs.
This is primarily for loop analysis but is also used for nodal analysis.

Copyright 2019 Michael Hayes, UCECE

"""

from matplotlib.pyplot import subplots
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
            if len(elt.nodes) < 2:
                continue
            self.add_edge(node_map[elt.nodes[0]],
                          node_map[elt.nodes[1]], name=name)

        self.node_map = node_map

    def connected_elts(self, node):

        for node1, edges in self.node_edges(node).items():
            for key, edge in edges.items():
                name = edge['name']
                elt = self.cct.elements[name]
                yield elt
            
    def loops(self):
        
        DG = nx.MultiDiGraph(self)
        cycles = list(nx.simple_cycles(DG))

        loops = []
        for cycle in cycles:
            if len(cycle) > 2:
                cycle = sorted(cycle)
                if cycle not in loops:
                    loops.append(cycle)        
        return loops

    def draw(self):
        """Use matplotlib to draw circuit graph."""

        fig, ax = subplots(1)

        G = self
        pos = nx.spring_layout(G)
        
        labels = dict(zip(G.nodes(), G.nodes()))
        nx.draw_networkx(G, pos, ax, labels=labels)
        
        edge_labels = dict([((u, v), d['name'])
                            for u, v, d in G.edges(data=True)])
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

    def node_edges(self, node):

        return self[node]
        
# neighbors(node) gives neighbouring nodes
