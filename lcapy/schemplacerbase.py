"""This module provides a base class for schematic placement algorithms.

Copyright 2021 Michael Hayes, UCECE
"""

from .schemmisc import Pos
from numpy import argsort


class SchemPlacerBase(object):

    def _xlink(self, elt):

        xvals = elt.xvals
        nodes = elt.nodes
        for m1, n1 in enumerate(nodes):
            for m2, n2 in enumerate(nodes[m1 + 1:], m1 + 1):
                if xvals[m2] == xvals[m1]:
                    self.xgraph.link(n1.name, n2.name)

    def _ylink(self, elt):

        yvals = elt.yvals
        nodes = elt.nodes
        for m1, n1 in enumerate(nodes):
            for m2, n2 in enumerate(nodes[m1 + 1:], m1 + 1):
                if yvals[m2] == yvals[m1]:
                    self.ygraph.link(n1.name, n2.name)

    def _place(self, elt, graph, vals):

        if elt.free:
            return

        if elt.offset != 0:
            print('TODO: offset %s by %f' % (elt, elt.offset))

        size = elt.size
        nodes = elt.nodes
        idx = argsort(vals)[::-1]
        for i in range(len(idx) - 1):
            m1 = idx[i]
            m2 = idx[i + 1]
            n1 = nodes[m1]
            n2 = nodes[m2]
            value = (vals[m2] - vals[m1]) * size
            graph.add(elt, n1.name, n2.name, value, elt.stretch)

    def _xplace(self, elt):

        self._place(elt, self.xgraph, elt.xvals)

    def _yplace(self, elt):

        self._place(elt, self.ygraph, elt.yvals)

    def _make_graphs(self):

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

        # Use components in orthogonal directions as constraints.  The
        # nodes of orthogonal components get combined into a
        # common node.
        for m, elt in enumerate(self.elements.values()):

            if elt.offset != 0:
                raise ValueError('offset field should be removed')
            if elt.directive or elt.ignore or not elt.place or elt.free:
                continue

            self._xlink(elt)
            self._ylink(elt)

        # Now form forward and reverse directed graph using components
        # in the desired directions.
        # Note, this must be done after the linking step.
        for m, elt in enumerate(self.elements.values()):
            if elt.directive or elt.ignore or elt.free or not elt.place:
                continue
            self._xplace(elt)
            self._yplace(elt)

    def solve(self, node_spacing):

        self._make_graphs()

        xpos, width = self.xgraph.solve()
        ypos, height = self.ygraph.solve()

        scale = node_spacing
        for n, node in self.nodes.items():
            node.pos = Pos(xpos[n] * scale, ypos[n] * scale)

        return width, height
