"""This module handles the nodes for a circuit.

Copyright 2023 Michael Hayes, UCECE
"""

from .attrdict import AttrDict
from .schemmisc import Pos
from .node import Node


class Nodes(AttrDict):
    """An instance of this stores all the defined Nodes for a circuit as a
    dictionary.  Elements in the direction can also be accessed as
    attributes.  Each Node has a list of associated components.

    """

    def add(self, node_name, cpt, cct):
        """Add a new node with associated component and circuit."""

        if node_name in self:
            node = self[node_name]
        else:
            node = Node(cct, node_name)
            self[node_name] = node
        node.append(cpt)
        return node

    def remove(self, node_name, cpt):
        """Detach node from component and delete the node
        if there are no more uses."""

        if node_name not in self:
            raise ValueError('Cannot remove unknown node', node_name)
        node = self[node_name]
        node.remove(cpt)

    def rename(self, old_node_name, new_node_name, cpts=None):
        """Rename old node name to new node name if the old node is connected
        to a component in the list cpts or if cpts is None.  Returns
        the new node."""

        if old_node_name not in self:
            raise ValueError('Cannot rename unknown node', old_node_name)
        node = self[old_node_name]
        return node.rename(new_node_name, cpts)

    def _delete(self, node_name):
        """Delete node from dict of nodes.  Note, it is best to use remove
        since this will not delete the Node until all uses are
        removed."""

        node = self[node_name]
        if len(node.connected) != 0:
            raise ValueError('Cannot remove node ' + str(node_name) +
                             '; it is needed for ' +
                             ', '.join([str(name) for name in node.connected]))

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
        if entry == '':
            continue
        parts = entry.split('@')
        node_name = parts[0].strip()
        # Ignore ()
        values = parts[1][1:-1]
        parts = values.split(',')
        x = float(parts[0])
        y = float(parts[1])
        pos = Pos(x, y)

        node_positions[node_name] = pos

    return node_positions
