"""This module provides support for nodes.

Copyright 2020--2024 Michael Hayes, UCECE

"""

from .immittancemixin import ImmittanceMixin
from copy import copy


class DummyNode:

    def __init__(self, name):
        self.name = name


class Node(ImmittanceMixin):
    """Each node has an instance of this.   """

    def __init__(self, cct, name):

        self.cct = cct
        self._name = name
        self.pos = None
        parts = name.split('_')
        self.rootname = parts[0] if name[0] != '_' else name
        self.primary = len(parts) == 1
        # List of components connected to this node.
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
        """Rename node name."""

        self.rename(name)

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

    @property
    def connected_nodes(self):
        """Return list of nodes that are connected to the node
        by a component."""

        nodes = []
        for cpt in self._connected:
            for node in cpt.nodes:
                if node is not self and node not in nodes:
                    nodes.append(node)
        return nodes

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
        """Return True if node is a port."""

        for cpt in self.connected:
            if cpt.type == 'P':
                return True
        return False

    @property
    def port(self):
        """Alias for is_port."""

        return self.is_port

    @property
    def is_dangling(self):
        """Return True if node has a single electrical connection"""

        return self._count <= 1

    @property
    def is_ground(self):
        """Return True if node is a ground."""

        return self.name.startswith('0')

    def append(self, cpt):
        """Attach component to the node."""

        if cpt.type not in ('A', 'O'):
            self._count += 1

        self._connected.append(cpt)

    def clone(self, name):
        """Create a copy of the node but with a new name."""

        node = Node(self.cct, name)
        node.pos = copy(self.pos)
        return node

    def remove(self, cpt):
        """Detach component from node.  If no other components
        are attached, the node is deleted."""

        if len(self._connected) == 0:
            raise RuntimeError(
                'Removing node %s with no connections' % str(self))

        for c in self._connected:
            if c.name == cpt.name:
                self._connected.remove(c)
                break
        if cpt.type not in ('A', 'O'):
            self._count -= 1

        if self.count == 0:
            self.cct.nodes._delete(self.name)

    def _rename1(self, name):
        # Case 1.  The name is new and is applied to all components

        if name in self.cct.nodes:
            raise ValueError('Node %s is not new' % name)

        self.cct.nodes[name] = self.cct.nodes.pop(self._name)
        self._name = name
        return self

    def _rename2(self, name):
        # Case 2.  The name is not new and is applied to all components

        if name not in self.cct.nodes:
            raise ValueError('Node %s is new' % name)

        self._connected.extend(self.cct.nodes[name]._connected)
        self._count += self.cct.nodes[name]._count

        self.cct.nodes[name] = self.cct.nodes.pop(self._name)
        self._name = name
        return self

    def _rename3(self, name, cpts):
        # Case 3.  The name is new and is applied to some of the components

        if name in self.cct.nodes:
            raise ValueError('Node %s is not new' % name)

        new_node = self.clone(name)
        self.cct.nodes[name] = new_node

        for cpt in cpts:
            for m, node in enumerate(cpt.nodes):
                if node.name == self.name:
                    if cpt.type not in ('A', 'O'):
                        self._count -= 1
                    new_node.append(cpt)
                    cpt.nodes[m] = new_node
        return new_node

    def _rename4(self, name, cpts):
        # Case 4.  The name is not new and is applied to some of the components

        if name not in self.cct.nodes:
            raise ValueError('Node %s is new' % name)

        new_node = self.cct.nodes[name]
        self.cct.nodes[name] = new_node

        for cpt in cpts:
            for m, node in enumerate(cpt.nodes):
                if node.name == self.name:
                    if cpt.type not in ('A', 'O'):
                        self._count -= 1
                    new_node.append(cpt)
                    cpt.nodes[m] = new_node
        return new_node

    def rename(self, name, cpts=None):
        """Rename node name if connected to a component in the list cpts or if
        cpts is None.  Return new node."""

        # Case 1. The name is new and is applied to all components
        # Case 2. The name is not new and is applied to all components
        # Case 3. The name is new and is applied to some of the components
        # Case 4. The name is not new and is applied to some of the components
        #
        # all and some applies to the components that are connected
        # to this node

        if cpts is None:
            if name not in self.cct.nodes:
                # Case 1.
                return self._rename1(name)
            else:
                # Case 2.
                return self._rename2(name)

        # Convert names to cpts
        if not isinstance(cpts, (list, tuple)):
            cpts = [cpts]

        tcpts = []
        for cpt in cpts:
            if isinstance(cpt, str):
                tcpts.append(self.cct[cpt])
            else:
                tcpts.append(cpt)
        cpts = tcpts

        other = []
        for cpt in self.connected:
            if cpt not in cpts:
                other.append(cpt)

        if other == []:
            # Case 1
            if name not in self.cct.nodes:
                return self._rename1(name)

            # Case 2
            return self._rename2(name)

        # Case 3
        if name not in self.cct.nodes:
            return self._rename3(name, cpts)

        # Case 4
        return self._rename4(name, cpts)

    @property
    def count(self):
        """Number of components electrically connected to the node.  This
        includes wires but not annotations.

        """

        return self._count

    def debug(self):

        s = str(self) + ', count=%s' % self.count

        names = [cpt.name for cpt in self.connected]

        s += ', cpts=[%s]' % ', '.join(names) + '\n'

        return s
