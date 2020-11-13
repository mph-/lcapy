"""
This module provides a class to represent circuits as graphs.
This is primarily for loop analysis but is also used for nodal analysis.

Copyright 2019--2020 Michael Hayes, UCECE

"""

from matplotlib.pyplot import subplots, savefig
import networkx as nx


# V1 1 0 {u(t)}; down
# R1 1 2; right=2
# L1 2 3; down=2
# W1 0 3; right
# W 1 5; up
# W 2 6; up
# C1 5 6; right=2

class CircuitGraph(nx.Graph):

    def __init__(self, cct):

        super(CircuitGraph, self).__init__()
        self.cct = cct
        self.dummy = 0
        # Dummy nodes are used to avoid parallel edges.
        self.dummy_nodes = {}

        self.add_nodes_from(cct.node_list)

        node_map = cct.node_map

        for name in cct.branch_list:
            elt = cct.elements[name]
            if len(elt.nodenames) < 2:
                continue

            nodename1 = node_map[elt.nodenames[0]]
            nodename2 = node_map[elt.nodenames[1]]
            
            if self.has_edge(nodename1, nodename2):
                # Add dummy node in graph to avoid parallel edges.
                dummynode = '*%d' % self.dummy
                dummycpt = 'W%d' % self.dummy                
                self.add_edge(nodename1, dummynode, name=name)
                self.add_edge(dummynode, nodename2, name=dummycpt)
                self.dummy_nodes[dummynode] = nodename2
                self.dummy += 1
            else:
                self.add_edge(nodename1, nodename2, name=name)

        self.node_map = node_map

    def connected(self, node):
        """List of components connected to specified node."""

        for node1, edges in self.node_edges(node).items():
            if node1.startswith('*'):
                for elt in self.connected(node1):
                    yield elt
                continue
            
            for key, edge in edges.items():
                if not edge.startswith('W'):
                    elt = self.cct.elements[edge]
                    yield elt
            
    def all_loops(self):

        # This adds forward and backward edges.
        DG = nx.DiGraph(self)
        cycles = list(nx.simple_cycles(DG))

        loops = []
        for cycle in cycles:
            if len(cycle) <= 2:
                continue
            cycle = sorted(cycle)
            if cycle not in loops:
                loops.append(cycle)        
        return loops

    def chordless_loops(self):

        loops = self.all_loops()
        sets = [set(loop) for loop in loops]

        DG = nx.DiGraph(self)
        
        rejects = []
        for i in range(len(sets)):

            # Reject loops with chords.
            loop = loops[i]
            if len(loop) == 2:
                continue
            
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

        name = self.get_edge_data(node1, node2)['name']
        if name.startswith('W'):
            return None
        return self.cct.elements[name]

    def loops_for_cpt(self, elt):
        """Return list of tuples; one for each loop.  The first element of the
        tuple is the loop the cpt belongs to or an empty list; the
        second element indicates the cpt direction compared to the
        loop direction."""
        
        loops = self.loops()
        cloops = []

        # Map node names to equipotential node names.
        nodenames = [self.cct.node_map[nodename] for nodename in elt.nodenames]
        
        for n, loop in enumerate(loops):

            loop1 = loop.copy()
            loop1.append(loop1[0])

            def find(loop1, nodename1, nodename2):
                for m in range(len(loop1) - 1):
                    if (nodename1 == loop1[m] and
                        nodename2 == loop1[m + 1]):
                        return True

            if find(loop1, nodenames[0], nodenames[1]):
                cloops.append((loop, False))
            elif find(loop1, nodenames[1], nodenames[0]):
                cloops.append((loop, True))
            else:
                cloops.append(([], None))

        return cloops
    
# neighbors(node) gives neighbouring nodes
