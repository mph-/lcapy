"""This module provides node support for schematic drawing.

Copyright 2022 Michael Hayes, UCECE

"""

from .latex import latex_format_node_label


class Node:

    def __init__(self, name):

        self.name = name
        self._port = False
        self._count = 0
        parts = name.split('_')
        # Primary 1, 2, a, a_3.  Not primary _1, _2, _a, _a_3, 0_3, a_b_c
        self.primary = (name[0] != '_' and len(parts) <= 2) and not (
            name[0].isdigit() and len(parts) != 1)

        # See MX1.sch for example.
        if '._' in name:
            self.primary = False

        self.elt_list = []
        self.pos = 'unknown'
        # Sanitised name
        self.s = name.replace('.', '@')
        self.label = latex_format_node_label(self.name)
        self.pin = False
        self.pinlabel = ''
        self.namepos = None
        self.pinname = ''
        self.pinpos = None
        self.clock = False
        self.auxiliary = None
        # Reference to node that defines the pin
        self.ref = None
        # Node has an implicit connection, say to ground.
        self.implicit = False
        # Symbol to use if node is an autoground
        self._symbol = None
        self.split_count = 0

    @property
    def basename(self):
        fields = self.name.split('.')
        return fields[-1]

    @property
    def cptname(self):
        fields = self.name.split('.')
        if len(fields) < 2:
            return None
        return fields[-2]

    def __repr__(self):
        return '%s @ (%s)' % (self.name, self.pos)

    def append(self, elt):
        """Add new element to the node"""

        if elt.type == 'P':
            self._port = True

        self.elt_list.append(elt)
        if elt.type not in ('A', 'O'):
            self._count += 1

    @property
    def count(self):
        """Number of elements (including wires but not open-circuits and
        annotations) connected to the node"""

        return self._count

    def belongs(self, cpt_name):
        return self.cptname == cpt_name

    def visible(self, draw_nodes):
        """Return True if node drawn.
        `draw_nodes' can be `all', 'none', 'connections', 'primary', None,
        True, or False."""

        if self.auxiliary:
            return False

        if draw_nodes in ('all', True):
            return True

        if self.pin:
            return False

        if self.implicit:
            return False

        if self._port:
            return True

        if draw_nodes in ('none', None, False):
            return False

        # Implied port
        if self._count == 1:
            return True

        if draw_nodes in ('connections', 'connected'):
            return self._count > 2

        if draw_nodes == 'primary':
            return self.primary

        raise ValueError('Unknown argument %s for draw_nodes' % draw_nodes)

    @property
    def is_port(self):
        """Return True if node is a port"""

        return self._port

    @property
    def is_dangling(self):
        """Return True if node has a single connection"""

        return self._count == 1

    @property
    def is_ground(self):
        """Return True if node is a ground"""

        return self.name.startswith('0')

    def show_label(self, label_nodes):

        if self.label == '':
            return False

        name = self.basename

        # pins is for backward compatibility
        if label_nodes in ('none', 'pins', 'false', False):
            return False
        elif label_nodes in ('all', 'true', True):
            return True
        elif label_nodes == 'alpha':
            return self.primary and name[0].isalpha()
        elif label_nodes == 'primary':
            return self.primary

        # handle label_nodes = '{1, 2}' etc.
        if label_nodes[0] == '{' and label_nodes[-1] == '}':
            label_nodes = label_nodes[1:-1]
        labels = [foo.strip() for foo in label_nodes.split(',')]
        return self.label in labels

    def debug(self):
        print(' %s @ (%s), count=%d, pin=%s' % (self.name, self.pos,
                                                self._count, self.pin))

    def split(self, elt):

        self.split_count += 1
        if self.split_count == 1:
            return self

        name = self.name + '_split%d' % self.split_count

        new_node = Node(name)
        self._count -= 1
        new_node._count = 1
        new_node._port = self._port
        new_node.label = self.label
        new_node._port = elt.type == 'P'

        return new_node

    @property
    def symbol(self):
        """Return node symbol"""

        return self._symbol

    @symbol.setter
    def symbol(self, value):
        """Set node symbol."""

        if self._symbol != None and self._symbol != value:
            raise ValueError('Conflicting symbols %s and %s' %
                             (self._symbol, value))
        self._symbol = value
