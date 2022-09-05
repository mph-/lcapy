"""This module provides support for nodes.

Copyright 2020--2022 Michael Hayes, UCECE

"""

from .immittancemixin import ImmittanceMixin


class Node(ImmittanceMixin):

    def __init__(self, cct, name):

        self.cct = cct
        self.name = name
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
        if cpt.type not in ('A', 'O'):
            self._count += 1

        self._connected.append(cpt)

    @property
    def count(self):
        """Number of elements (including wires but not open-circuits and
        annotations) connected to the node"""

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
        """Return True if node has a single connection"""

        return self._count <= 1

    @property
    def is_ground(self):
        """Return True if node is a ground"""

        return self.name.startswith('0')
