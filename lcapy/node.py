"""This module provides support for nodes.

Copyright 2020--2023 Michael Hayes, UCECE

"""

from .immittancemixin import ImmittanceMixin


class DummyNode:

    def __init__(self, name):
        self.name = name


class Node(ImmittanceMixin):

    def __init__(self, cct, name):

        self.cct = cct
        self._name = name
        self.pos = None
        self.port = False
        parts = name.split('_')
        self.rootname = parts[0] if name[0] != '_' else name
        self.primary = len(parts) == 1
        # List of elements connected to this node.
        self._connected = []
        self._count = 0

    def __repr__(self):
        return "Node('%s')" % self.name

    def __str__(self):

        if self.pos is not None:
            # Note, need to have float type otherwise 0 becomes empty string.
            x = str(round(self.x, 2)).rstrip('0').rstrip('.')
            y = str(round(self.y, 2)).rstrip('0').rstrip('.')

            return '%s@(%s, %s)' % (self.name, x, y)
        else:
            return self.name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        # Use new node name for nodes
        self.cct.nodes[name] = self.cct.nodes.pop(self._name)
        self._name = name

    @property
    def x(self):

        if self.pos is None:
            raise ValueError('Node position undefined')
        return self.pos.x

    @property
    def y(self):

        if self.pos is None:
            raise ValueError('Node position undefined')
        return self.pos.y

    @property
    def V(self):
        """Node voltage with respect to ground."""

        return self.cct.get_Vd(self.name, '0')

    @property
    def v(self):
        """Node time-domain voltage with respect to ground."""

        return self.cct.get_vd(self.name, '0')

    @property
    def admittance(self):
        """Driving-point admittance between node and ground."""

        return self.cct.admittance(self.name, '0')

    @property
    def impedance(self):
        """Driving-point impedance between node and ground."""

        return self.cct.impedance(self.name, '0')

    def append(self, cpt):

        if cpt.type in ('P', ):
            self.port = True
        if cpt.type not in ('A', 'O', 'P'):
            self._count += 1

        self._connected.append(cpt)

    def remove(self, cpt):

        if len(self._connected) == 0:
            raise RuntimeError(
                'Removing node %s with no connections' % str(self))

        for c in self._connected:
            if c.name == cpt.name:
                self._connected.remove(c)
                break
        if cpt.type not in ('A', 'O', 'P'):
            self._count -= 1

        if self.count == 0:
            self.cct.nodes._remove(self.name)

    @property
    def count(self):
        """Number of elements electrically connected to the node.  This
        includes wires but not open-circuits, ports, and annotations.

        """

        return self._count

    def oneport(self, node=0):
        """Create oneport object with respect to specified node
        (default ground)."""

        return self.cct.oneport(self.name, node)

    def thevenin(self, node=0):
        """Create Thevenin oneport object with respect to specified
        node (default ground)."""

        return self.cct.thevenin(self.name, node)

    def norton(self, node=0):
        """Create Norton oneport object with respect to specified node
        (default ground)."""

        return self.cct.norton(self.name, node)

    @property
    def connected(self):
        """Return list of components connected to the node."""

        return self._connected

    def is_connected(self, cpt):
        """Return True if cpt is connected to the node."""

        if isinstance(cpt, str):
            for cpt1 in self.connected:
                if cpt1.name == cpt:
                    return True
            return False

        return cpt in self.connected

    def is_wired_to(self, node):
        """Return True if the nodes are wired together."""

        node = str(node)
        for nodes in self.cct.equipotential_nodes.values():
            if self.name in nodes and node in nodes:
                return True
        return False

    def wired_to(self):
        """Return list of names of nodes that are wired to this node."""

        for nodes in self.cct.equipotential_nodes.values():
            if self.name in nodes:
                nodes.remove(self.name)
                return nodes
        return []

    @property
    def is_port(self):
        """Return True if node is a port"""

        return self._port

    @property
    def is_dangling(self):
        """Return True if node has a single electrical connection"""

        return self._count <= 1

    @property
    def is_ground(self):
        """Return True if node is a ground"""

        return self.name.startswith('0')

    def debug(self):

        s = str(self) + ', count=%s' % self.count

        names = [cpt.name for cpt in self.connected]

        s += ', cpts=[%s]' % ', '.join(names) + '\n'

        return s
