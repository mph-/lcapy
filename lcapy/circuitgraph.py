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
                if sets[i].issubset(sets[j]):
                    rejects.append(j)
                elif sets[j].issubset(sets[i]):
                    rejects.append(i)

        cloops = []
        for i, loop in enumerate(loops):
            if i not in rejects:
                cloops.append(loop)

        return cloops

    def loops(self):
        if hasattr(self, '_loops'):
            return self._loops
        self._loops = self.chordless_loops()
        return self._loops
    
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

    @property
    def is_planar(self):
        """Return True for a planar network."""
        
        return nx.check_planarity(self)[0]

    def node_edges(self, node):

        return self[node]

    def component(self, node1, node2):
        
        return self.cct.elements[self.get_edge_data(node1, node2)[0]['name']]

    def loop_indices_for_cpt(self, elt):
        """Return list of tuples.  The first element of the tuple
        is the loop index the cpt belongs to; the second element indicates
        the cpt direction."""
        
        loops = self.loops()
        cloops = []

        # Map node names to equipotential node names.
        nodenames = [self.cct.node_map[nodename] for nodename in elt.nodenames]
        
        for n, loop in enumerate(loops):

            loop1 = loop.copy()
            loop1.append(loop1[0])

            for m in range(len(loop1) - 1):
                if (nodenames[0] == loop1[m] and
                    nodenames[1] == loop1[m + 1]):
                    cloops.append((n, True))
                    break
                elif (nodenames[1] == loop1[m] and
                      nodenames[0] == loop1[m + 1]):
                    cloops.append((n, False))
                    break            

        return cloops
    
# neighbors(node) gives neighbouring nodes
