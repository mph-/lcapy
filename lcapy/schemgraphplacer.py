"""
This module performs schematic component placement using a graph algorithm.

Copyright 2021 Michael Hayes, UCECE
"""

from .schemgraph import Graph
from .schemmisc import Pos

# Components are positioned using two graphs; one graph for
# the x direction and the other for the y direction. 
#
# There is a naming confusion.  We have network nodes (electrical
# nodes) and nodes in the graph used for component placement.
# Let's call the latter gnodes.  The names of these gnodes are a
# tuple of the common network nodes.
#
# x and y component positioning are performed independently.  Let's
# consider the x or horizontal positioning.  There are three stages:
#   1. Component nodes that share a y position are linked; this can
#      occur, for example, for a vertically oriented component.
#      This helps to reduce the size of the graph.
#   2. The x positions of the components are used to determine the
#      graph edges.
#   3. The longest path through the graph is found and the x positions
#      of the nodes are assigned based on the distance along the
#      longest path.

class SchemGraphPlacer(object):

    def __init__(self, elements, nodes):

        self.elements = elements
        self.nodes = nodes


    def _xlink(self, elt, graph):

        if not elt.place or elt.free:
            return
        
        xvals = elt.xvals
        nodes = elt.nodes
        for m1, n1 in enumerate(nodes):
            for m2, n2 in enumerate(nodes[m1 + 1:], m1 + 1):
                if xvals[m2] == xvals[m1]:
                    graph.link(n1.name, n2.name)

    def _ylink(self, elt, graph):

        if not elt.place or elt.free:
            return        

        yvals = elt.yvals
        nodes = elt.nodes        
        for m1, n1 in enumerate(nodes):
            for m2, n2 in enumerate(nodes[m1 + 1:], m1 + 1):
                if yvals[m2] == yvals[m1]:
                    graph.link(n1.name, n2.name)
        
    def _make_graphs(self, debug=None):

        # The x and y positions of a component node are determined
        # independently.  The principle is that each component has a
        # minimum size (usually 1 but changeable with the size option)
        # but its wires can be stretched.

        # When solving the x position, first nodes that must be
        # vertically aligned (with the up or down option) are combined
        # into a set.  Then the left and right options are used to
        # form a graph.  This graph is traversed to find the longest
        # path and in the process each node gets assigned the longest
        # distance from the root of the graph.  To centre components,
        # a reverse graph is created and the distances are averaged.

        self.xgraph = Graph('horizontal', self.nodes, debug)
        self.ygraph = Graph('vertical', self.nodes, debug)

        # Use components in orthogonal directions as constraints.  The
        # nodes of orthogonal components get combined into a
        # common node.
        for m, elt in enumerate(self.elements.values()):

            if elt.offset != 0:
                raise ValueError('offset field should be removed')
            if elt.directive or elt.ignore:
                continue
            
            self._xlink(elt, self.xgraph)
            self._ylink(elt, self.ygraph)

        # Now form forward and reverse directed graph using components
        # in the desired directions.
        # Note, this must be done after the linking step.
        for m, elt in enumerate(self.elements.values()):
            if elt.directive or elt.ignore:
                continue            
            elt.xplace(self.xgraph)
            elt.yplace(self.ygraph)
            
    def positions_calculate(self, node_spacing):

        self._make_graphs()

        xpos, width = self.xgraph.analyse()
        ypos, height = self.ygraph.analyse()

        scale = node_spacing
        for n, node in self.nodes.items():
            node.pos = Pos(xpos[n] * scale, ypos[n] * scale)

        return width, height
