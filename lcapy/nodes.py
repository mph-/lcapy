"""This module handles the nodes for a circuit.

Copyright 2023 Michael Hayes, UCECE
"""

from .attrdict import AttrDict


class Nodes(AttrDict):

    def _remove(self, node_name):
        """Remove node from dict of nodes.  Note, it is best to use Node.remove."""

        node = self[node_name]
        if node.count != 0:
            raise ValueError('Cannot remove node; it is needed for ' +
                             ', '.join(node.components))

        self.pop(node_name)

    def debug(self):

        s = ''
        for node in self:
            s += node.debug()
        return s

    def new_name(self):

        num = 1
        while True:
            name = str(num)
            if not name in self:
                return name
            num += 1

    def closest(self, x, y):

        for node in self:
            x1, y1 = node.pos
            rsq = (x1 - x)**2 + (y1 - y)**2
            if rsq < 0.1:
                return node
        return None

    def by_position(self, position):

        x, y = position

        for node in self:
            if abs(node.x - x) < 1e-5 and abs(node.y - y) < 1e-5:
                return node
        return None
