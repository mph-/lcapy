"""This module handles the nodes for a circuit.

Copyright 2023 Michael Hayes, UCECE
"""

from .attrdict import AttrDict
from .schemmisc import Pos


class Nodes(AttrDict):

    def _remove(self, node_name):
        """Remove node from dict of nodes.  Note, it is best to use Node.remove;
        this will not remove the Node until all uses are removed."""

        node = self[node_name]
        if node.count != 0:
            raise ValueError('Cannot remove node; it is needed for ' +
                             ', '.join(node.components))

        self.pop(node_name)

    def debug(self):

        s = ''
        for node in self.values():
            s += node.debug()
        return s

    def by_position(self, position):

        x, y = position

        for node in self.values():
            if abs(node.x - x) < 1e-5 and abs(node.y - y) < 1e-5:
                return node
        return None


def parse_nodes(nodesstr):

    from .utils import split_parens

    node_positions = {}

    # Ignore {}
    nodesstr = nodesstr[1:-1]
    entries = split_parens(nodesstr, ',')
    for entry in entries:
        parts = entry.split('@')
        nodename = parts[0].strip()
        # Ignore ()
        values = parts[1][1:-1]
        parts = values.split(',')
        x = float(parts[0])
        y = float(parts[1])
        pos = Pos(x, y)

        node_positions[nodename] = pos

    return node_positions
