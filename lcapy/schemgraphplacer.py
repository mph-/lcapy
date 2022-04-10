"""
This module performs schematic component placement using a graph algorithm.

Copyright 2021 Michael Hayes, UCECE
"""

from .schemgraph import Graph
from .schemmisc import Pos
from .schemplacerbase import SchemPlacerBase

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


class SchemGraphPlacer(SchemPlacerBase):

    def __init__(self, elements, nodes, debug=0):

        self.elements = elements
        self.nodes = nodes
        self.debug = debug

        self.xgraph = Graph('horizontal', nodes, debug)
        self.ygraph = Graph('vertical', nodes, debug)
