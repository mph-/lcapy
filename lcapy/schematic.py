"""
This module performs schematic drawing using circuitikz from a netlist.

>>> from lcapy import Schematic
>>> sch = Schematic()
>>> sch.add('P1 1 0.1; down')
>>> sch.add('R1 3 1; right')
>>> sch.add('L1 2 3; right')
>>> sch.add('C1 3 0; down')
>>> sch.add('P2 2 0.2; down')
>>> sch.add('W 0 0.1; right')
>>> sch.add('W 0.2 0.2; right')
>>> sch.draw()

Copyright 2014, 2015 Michael Hayes, UCECE
"""

from __future__ import print_function
import numpy as np
import re
from lcapy.core import Expr
from os import system, path, remove


__all__ = ('Schematic', )


# Regular expression alternate matches stop with first match so need
# to have longer names first.
cpt_types = ['R', 'C', 'L', 'Z', 'Y', 'V', 'I', 'W', 'P', 'E', 
             'D', 'J', 'M', 'Q', 'TF', 'TP', 'K']
cpt_types.sort(lambda x, y: cmp(len(y), len(x)))

cpt_type_pattern = re.compile(r"(%s)([\w']*)" % '|'.join(cpt_types))

label_pattern = re.compile(r"([\w']*)(_{[\w]})?")


import math


class EngFormat(object):

    def __init__(self, value, unit=''):

        self.value = value
        self.unit = unit

    def latex(self, trim=True, hundreds=False):

        prefixes = ('f', 'p', 'n', '$\mu$', 'm', '', 'k', 'M', 'G', 'T')

        sfmax = 3

        value = self.value
        m = math.log10(abs(value))

        if hundreds:
            # Generate 100 m
            n = int(math.floor(m / 3))
            k = int(math.floor(m)) - n * 3
        else:
            # Generate 0.1
            n = int(round(m / 3))
            k = int(round(m)) - n * 3

        dp = sfmax - k

        idx = n + 5
        if idx < 0:
            idx = 0
            return '%e\,' % value + self.unit
        elif idx >= len(prefixes):
            idx = len(prefixes) - 1
            return '%e\,' % value + self.unit

        fmt = '%%.%df' % dp

        n = idx - 5
        value = value * 10**(-3 * n)

        string = fmt % value

        if trim:
            # Remove trailing zeroes after decimal point
            string = string.rstrip('0').rstrip('.')

        return string + '\,' + prefixes[idx] + r'\mbox{' + self.unit + r'}'


class Opts(dict):

    def _parse(self, string):

        for part in string.split(','):
            part = part.strip()
            if part == '':
                continue

            if part in ('up', 'down', 'left', 'right'):
                self['dir'] = part
                continue

            fields = part.split('=')
            key = fields[0].strip()
            arg = fields[1].strip() if len(fields) > 1 else ''
            self[key] = arg

    def __init__(self, arg):

        if isinstance(arg, str):
            self._parse(arg)
            return

        for key, val in arg.iteritems():
            self[key] = val

    def format(self):

        return ', '.join(['%s=%s' % (key, val)
                          for key, val in self.iteritems()])


class Cnodes(object):

    """Common nodes"""

    def __init__(self, nodes):

        self.sets = {}
        for node in nodes:
            self.sets[node] = (node, )

    def link(self, n1, n2):
        """Make nodes n1 and n2 share common node"""

        set1 = self.sets[n1]
        set2 = self.sets[n2]
        newset = set1 + set2

        for n in self.sets[n1]:
            self.sets[n] = newset
        for n in self.sets[n2]:
            self.sets[n] = newset

    def _analyse(self):

        # Add dymmy cnode at start
        unique = ['dummy'] + list(set(self.sets.values()))
        node_map = {}
        for node, nodes in self.sets.iteritems():
            node_map[node] = unique.index(nodes)

        self._node_map = node_map
        self._nodes = unique

    @property
    def node_map(self):
        """Return mapping of node number to common node number"""

        if not hasattr(self, '_node_map'):
            self._analyse()

        return self._node_map

    def map(self, nodes):

        if not isinstance(nodes, (tuple, list)):
            nodes = list[nodes]

        return [self.node_map[node] for node in nodes]

    @property
    def nodes(self):
        """Return mapping of common node number to tuple of shared nodes"""

        if not hasattr(self, '_nodes'):
            self._analyse()

        return self._nodes

    @property
    def size(self):
        """Return number of common nodes"""

        return len(self.nodes)


class Graph(dict):

    def __init__(self, size, name):

        self.name = name
        for m in range(size):
            self[m] = []

    def add(self, n1, n2, size):
        self[n1].append((n2, size))


    def longest_path(self):
        """Find longest path through DAG.  all_nodes is an iterable for all
        the nodes in the graph, from_nodes is a directory indexed by node
        that stores a tuple of tuples.  The first tuple element is the
        parent node and the second element is the minimium size of the
        component connecting the nodes.
        """

        all_nodes = self.keys()
        from_nodes = self

        memo = {}

        def get_longest(to_node):

            if to_node in memo:
                return memo[to_node]

            best = 0
            for from_node, size in from_nodes[to_node]:
                best = max(best, get_longest(from_node) + size)

            memo[to_node] = best

            return best

        try:
            length, node = max([(get_longest(to_node), to_node)
                                for to_node in all_nodes])
        except RuntimeError:
            raise RuntimeError(
                ('The schematic graph %s is dodgy, probably a component'
                 ' is connected to the wrong node\n%s') 
                % (self.name, from_nodes))

        return length, node, memo


class Graphs(object):

    def __init__(self, size, name):

        self.fwd = Graph(size, name + ' forward')
        self.rev = Graph(size, name + ' reverse')

    def add(self, n1, n2, size):
        self.fwd.add(n1, n2, size)
        self.rev.add(n2, n1, size)

    @property
    def nodes(self):
        return self.fwd.keys()


class Node(object):

    def __init__(self, name):

        self.name = name
        self._port = False
        self._count = 0
        self._cpt_count = 0
        parts = name.split('_')
        self.rootname = parts[0] if name[0] != '_' else name
        self.primary = len(parts) == 1
        self.list = []

    def append(self, elt):
        """Add new element to the node"""

        if elt.cpt_type == 'P':
            self._port = True

        self.list.append(elt)
        self._count += 1
        if elt.cpt_type != 'W':
            self._cpt_count += 1

    @property
    def count(self):
        """Number of elements (including wires) connected to the node"""

        return self._count

    @property
    def cpt_count(self):
        """Number of elements (not including wires) connected to the node"""

        return self._cpt_count

    @property
    def visible(self):
        """Return true if node drawn"""

        # Should show node when cpt_count >= 2
        # Could add blobs when count > 2

        return self.name.find('_') == -1

    @property
    def port(self):
        """Return true if node is a port"""

        return self._port or self.count == 1


class Pos(object):

    def __init__(self, x, y):

        self.x = x
        self.y = y

    def __mul__(self, scale):

        return Pos(self.x * scale, self.y * scale)

    def __str__(self):

        xstr = ('%.2f' % self.x).rstrip('0').rstrip('.')
        ystr = ('%.2f' % self.y).rstrip('0').rstrip('.')

        return "%s,%s" % (xstr, ystr)

    def __repr__(self):

        return 'Pos(%s)' % self

    @property
    def xy(self):

        return np.array((self.x, self.y))


class NetElement(object):

    cpt_type_counter = 0

    def __init__(self, name, n1, n2, *args, **opts):

        match = cpt_type_pattern.match(name)

        if not match:
            raise ValueError('Unknown schematic component %s' % name)

        # Circuitikz does not like a . in a name
        if n1.find('.') != -1:
            raise ValueError('Cannot have . in node name %s' % n1)
        if n2.find('.') != -1:
            raise ValueError('Cannot have . in node name %s' % n2)

        cpt_type = match.groups()[0]
        id = match.groups()[1]

        if id is None:
            NetElement.cpt_type_counter += 1
            id = '#%d' % NetElement.cpt_type_counter
            name = cpt_type + id

        # Component arguments
        self.args = args

        if args != ():
            if cpt_type in ('V', 'I') and args[0] in (
                    'ac', 'dc', 'step', 'acstep', 'impulse', 's'):
                cpt_type = cpt_type + args[0]
                args = args[1:]
            elif cpt_type == 'E' and args[0] == 'opamp':
                cpt_type = 'opamp'
                args = args[1:]

        # Tuple of nodes
        self.nodes = (n1, n2)
        # Identifier for component, e.g., 'R1'
        self.name = name
        # Type of component, e.g., 'V'
        self.cpt_type = cpt_type

        if cpt_type in ('J', 'M', 'Q'):
            if len(args) < 1:
                raise ValueError(
                    'Component type %s requires 3 nodes' % cpt_type)
            self.nodes += (args[0], )
            args = args[1:]
            
            if cpt_type == 'Q':
                self.sub_type = 'npn'
                if len(args) > 0:
                    if args[0] not in ('npn', 'pnp'):
                        raise ValueError('Bad argument %s for BJT' % args[0])
                    self.sub_type = args[0]
                args = args[1:]
            elif cpt_type == 'J':
                self.sub_type = 'njf'
                if len(args) > 0:
                    if args[0] not in ('npj', 'pjf'):
                        raise ValueError('Bad argument %s for JFET' % args[0])
                    self.sub_type = args[0]
                args = args[1:]
            elif cpt_type == 'M':
                self.sub_type = 'nmos'
                if len(args) > 0:
                    if args[0] not in ('nmos', 'pmos'):
                        raise ValueError('Bad argument %s for MOSFET' % args[0])
                    self.sub_type = args[0]
                args = args[1:]

        elif cpt_type in ('E', 'F', 'G', 'H', 'TF', 'TP', 'opamp'):
            if len(args) < 2:
                raise ValueError(
                    'Component type %s requires 4 nodes' % cpt_type)
            self.nodes += (args[0], args[1])
            args = args[2:]

        elif cpt_type == 'TP' and len(args) != 5:
            raise ValueError('TP component requires 5 args')

        elif cpt_type == 'D':
            self.sub_type = ''
            if len(args) > 0:
                if args[0] not in ('led', 'zener', 'photo', 'tunnel', 
                                   'schottky'):
                    raise ValueError('Bad argument %s for diode' % args[0])
                self.sub_type = args[0]
                args = args[1:]

        # There are two possible labels for a component:
        # 1. Component identifier, e.g., R1
        # 2. Component value, expression, or symbol
        identifier_label = name
        value_label = None

        if cpt_type in ('P', 'W') or identifier_label.find('#') != -1:
            identifier_label = None

        if 'dir' not in opts:
            opts['dir'] = None
        if 'size' not in opts:
            opts['size'] = 1

        if opts['dir'] is None:
            opts['dir'] = 'down' if cpt_type in ('P', ) else 'right'
        if len(args) > 0:

            # TODO, extend for mechanical and acoustical components.
            units_map = {'V': 'V', 'I': 'A', 'R': '$\Omega$',
                         'C': 'F', 'L': 'H'}

            expr = args[0]
            if cpt_type in ('Vimpulse', 'Iimpulse'):
                expr = '(%s) * DiracDelta(t)' % expr
                value_label = Expr(expr).latex()
            elif cpt_type in ('Vstep', 'Istep'):
                expr = '(%s) * Heaviside(t)' % expr
                value_label = Expr(expr).latex()
            elif cpt_type in ('Vs', 'Is'):
                value_label = Expr(expr).latex()
            elif cpt_type == 'TF':
                value_label = '1:%s' % args[0]
            elif cpt_type not in ('TP',):
                try:
                    value = float(args[0])
                    if cpt_type[0] in units_map:
                        value_label = EngFormat(
                            value, units_map[cpt_type[0]]).latex()
                    else:
                        value_label = Expr(expr).latex()

                except ValueError:
                    value_label = Expr(expr).latex()

        # Currently, we only annnotated the component with the value,
        # expression, or symbol.  If this is not specified, it
        # defaults to the component identifier.  Note, some objects
        # we do not want to label, such as wires and ports.

        tex_label = value_label
        if tex_label is None:
            if identifier_label is not None:
                tex_label = Expr(identifier_label).latex()

        # Default label to use when drawing with LaTeX
        if tex_label is None:
            self.tex_label = ''
        else:
            self.tex_label = '$%s$' % tex_label

        self.identifier_label = identifier_label
        self.value_label = value_label
        # Drawing hints
        self.opts = Opts(opts)

    def __repr__(self):

        str = ', '.join(arg.__str__()
                        for arg in [self.name] + list(self.nodes))
        return 'NetElement(%s)' % str

    def __str__(self):

        return ' '.join(['%s' % arg for arg in (self.name, ) + self.nodes])


class Schematic(object):

    def __init__(self, filename=None, **kwargs):

        self.elements = {}
        self.nodes = {}
        # Shared nodes (with same voltage)
        self.snodes = {}
        self.hints = False

        if filename is not None:
            self.netfile_add(filename)

    def __getitem__(self, name):
        """Return component by name"""

        return self.elements[name]

    def netfile_add(self, filename):
        """Add the nets from file with specified filename"""

        file = open(filename, 'r')

        lines = file.readlines()

        for line in lines:
            # Skip comments
            if line[0] in ('#', '%'):
                continue
            line = line.strip()
            if line != '':
                self.add(line)

    def netlist(self):
        """Return the current netlist"""

        return '\n'.join([elt.__str__() for elt in self.elements.values()])

    def _invalidate(self):

        for attr in ('_xnodes', '_ynodes', '_coords'):
            if hasattr(self, attr):
                delattr(self, attr)

    def _node_add(self, node, elt):

        if node not in self.nodes:
            self.nodes[node] = Node(node)
        self.nodes[node].append(elt)

        vnode = self.nodes[node].rootname

        if vnode not in self.snodes:
            self.snodes[vnode] = []

        if node not in self.snodes[vnode]:
            self.snodes[vnode].append(node)

    def _elt_add(self, elt):

        self._invalidate()

        if elt.name in self.elements:
            print('Overriding component %s' % elt.name)
            # Need to search lists and update component.

        self.elements[elt.name] = elt

        # Ignore nodes for mutual inductance.
        if elt.cpt_type == 'K':
            return

        for node in elt.nodes:
            self._node_add(node, elt)

    def add(self, string):
        """The general form is: 'Name Np Nm symbol'
        where Np is the positive nose and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        fields = string.split(';')
        string = fields[1].strip() if len(fields) > 1 else ''

        if string != '':
            self.hints = True

        opts = Opts(string)

        parts = re.split(r'[\s]+', fields[0].strip())
        elt = NetElement(*parts, **opts)

        self._elt_add(elt)


# Transformer
#   n4      n2
#
#   n3      n1
#
# For horiz. node layout want to make (n3, n4) and (n1, n2) in same cnodes
# For vert. node layout want to make (n2, n4) and (n1, n3) in same cnodes

    def _make_graphs(self, dirs):

        cnodes = Cnodes(self.nodes)

        # Use components in orthogonal directions as constraints.  The
        # nodes of orthogonal components get combined into a
        # common node.
        for m, elt in enumerate(self.elements.values()):

            if elt.cpt_type in ('TF', 'TP'):
                if dirs[0] == 'right':
                    # Ensure nodes have same x value
                    cnodes.link(*elt.nodes[0:2])
                    cnodes.link(*elt.nodes[2:4])
                else:
                    cnodes.link(*elt.nodes[0:4:2])
                    cnodes.link(*elt.nodes[1:4:2])
                continue
            elif elt.cpt_type == 'opamp':
                if dirs[0] == 'right':
                    # Ensure input nodes have same x value
                    cnodes.link(*elt.nodes[2:4])
                continue
            elif elt.cpt_type == 'K':

                # Should check that these inductors exist.
                L1 = self.elements[elt.nodes[0]]
                L2 = self.elements[elt.nodes[1]]

                # TODO, generalise
                if L1.opts['dir'] != 'down' or L2.opts['dir'] != 'down':
                    raise ValueError(
                        'Can only handle vertical mutual inductors')
                nodes = L2.nodes + L1.nodes
                n1, n2, n3, n4 = nodes

                # Provide horizontal constraints (the inductors
                # provide the vertical constraints).
                if dirs[0] != 'right':
                    cnodes.link(n3, n1)
                    cnodes.link(n4, n2)
                continue
            elif elt.cpt_type in ('J', 'M', 'Q'):
                n1, n2, n3 = elt.nodes                

                if dirs[0] == 'right':
                    # Horizontal constraint for C and E.
                    cnodes.link(n3, n1)
                continue

            if elt.opts['dir'] in dirs:
                continue

            cnodes.link(*elt.nodes[0:2])

        # Now form forward and reverse directed graphs using components
        # in the desired directions.
        graphs = Graphs(cnodes.size, 
                        'vertical' if dirs[0] == 'up' else 'horizontal')

        for m, elt in enumerate(self.elements.values()):

            size = float(elt.opts['size'])

            if elt.cpt_type in ('TF', 'TP'):
                # m1, m2 output nodes; m3, m4 input nodes
                m1, m2, m3, m4 = cnodes.map(elt.nodes)

                scale = {'TF': 0.5, 'TP': 2}

                if dirs[0] == 'right':
                    graphs.add(m3, m1, scale[elt.cpt_type] * size)
                    graphs.add(m4, m2, scale[elt.cpt_type] * size)
                else:
                    graphs.add(m2, m1, size)
                    graphs.add(m4, m3, size)
                continue
            elif elt.cpt_type == 'opamp':
                # m1, m2 output nodes; m3, m4 input nodes
                # m2 assumed ground...
                m1, m2, m3, m4 = cnodes.map(elt.nodes)

                if dirs[0] == 'right':
                    graphs.add(m3, m1, size * 2 + 0.4)
                    graphs.add(m4, m1, size * 2 + 0.4)
                else:
                    if 'mirror' in elt.opts:
                        graphs.add(m3, m1, 0.5 * size)
                        graphs.add(m1, m4, 0.5 * size)
                    else:
                        graphs.add(m4, m1, 0.5 * size)
                        graphs.add(m1, m3, 0.5 * size)
                continue
            elif elt.cpt_type == 'K':

                L1 = self.elements[elt.nodes[0]]
                L2 = self.elements[elt.nodes[1]]

                nodes = L2.nodes + L1.nodes
                m1, m2, m3, m4 = cnodes.map(nodes)

                scale = 0.8
                if dirs[0] == 'right':
                    graphs.add(m3, m1, scale * size)
                    graphs.add(m4, m2, scale * size)
                continue
            elif elt.cpt_type == 'J':
                # D, G, S
                m1, m2, m3 = cnodes.map(elt.nodes)

                yscale = 1.5
                xscale = 1.0
                if dirs[0] == 'right':
                    if elt.opts['dir'] == 'left':
                        m1, m2 = m2, m1
                    graphs.add(m2, m1, xscale * size)
                else:
                    if elt.sub_type == 'pjf':
                        graphs.add(m3, m2, yscale * size * 0.68)
                        graphs.add(m2, m1, yscale * size * 0.32)
                    else:
                        graphs.add(m3, m2, yscale * size * 0.32)
                        graphs.add(m2, m1, yscale * size * 0.68)
                continue
            elif elt.cpt_type in ('M', 'Q'):
                # C, B, E or D, G, S
                m1, m2, m3 = cnodes.map(elt.nodes)

                yscale = 1.5
                xscale = 0.8
                if dirs[0] == 'right':
                    if elt.opts['dir'] == 'left':
                        m1, m2 = m2, m1
                    graphs.add(m2, m1, xscale * size)
                else:
                    graphs.add(m3, m2, yscale * size * 0.5)
                    graphs.add(m2, m1, yscale * size * 0.5)
                continue

            if elt.opts['dir'] not in dirs:
                continue

            # Handle bipoles here.
            m1, m2 = cnodes.map(elt.nodes)

            if elt.opts['dir'] == dirs[0]:
                graphs.add(m1, m2, size)
            elif elt.opts['dir'] == dirs[1]:
                graphs.add(m2, m1, size)

        # Chain all potential start nodes to node 0.
        orphans = []
        rorphans = []
        for m in range(1, cnodes.size):
            if graphs.fwd[m] == []:
                orphans.append((m, 0))
            if graphs.rev[m] == []:
                rorphans.append((m, 0))
        graphs.fwd[0] = rorphans
        graphs.rev[0] = orphans

        if False:
            print(graphs.fwd)
            print(graphs.rev)
            print(cnodes.node_map)
            import pdb
            pdb.set_trace()

        # Find longest path through the graphs.
        length, node, memo = graphs.fwd.longest_path()
        length, node, memor = graphs.rev.longest_path()

        pos = {}
        posr = {}
        posa = {}
        for cnode in graphs.fwd.keys():
            if cnode == 0:
                continue

            for node in cnodes.nodes[cnode]:
                pos[node] = length - memo[cnode]
                posr[node] = memor[cnode]
                posa[node] = 0.5 * (pos[node] + posr[node])

        if False:
            print(pos)
            print(posr)
        return posa, cnodes.nodes

    def _positions_calculate(self):

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

        xpos, self._xnodes = self._make_graphs(('right', 'left'))
        ypos, self._ynodes = self._make_graphs(('up', 'down'))

        coords = {}
        for node in xpos.keys():
            coords[node] = Pos(xpos[node], ypos[node])

        self._coords = coords

    @property
    def xnodes(self):

        if not hasattr(self, '_xnodes'):
            self._positions_calculate()
        return self._xnodes

    @property
    def ynodes(self):

        if not hasattr(self, '_ynodes'):
            self._positions_calculate()
        return self._ynodes

    @property
    def coords(self):
        """Directory of position tuples indexed by node name"""

        if not hasattr(self, '_coords'):
            self._positions_calculate()
        return self._coords

    def _make_wires1(self, snode_list):

        num_wires = len(snode_list) - 1
        if num_wires == 0:
            return []

        wires = []

        # TODO: remove overdrawn wires...
        for n in range(num_wires):
            n1 = snode_list[n]
            n2 = snode_list[n + 1]

            wires.append(NetElement('W_', n1, n2))

        return wires

    def _make_wires(self):
        """Create implict wires between common nodes."""

        wires = []

        snode_dir = self.snodes

        for m, snode_list in enumerate(snode_dir.values()):
            wires.extend(self._make_wires1(snode_list))

        return wires

    def _node_str(self, n1, n2, draw_nodes=True):

        node1, node2 = self.nodes[n1], self.nodes[n2]

        if node1.port:
            node_str = 'o'
        else:
            node_str = '*' if draw_nodes and node1.visible else ''

        node_str += '-'

        if node2.port:
            node_str += 'o'
        else:
            node_str += '*' if draw_nodes and node2.visible else ''

        if node_str == '-':
            node_str = ''

        return node_str

    def _tikz_draw_opamp(self, elt, draw_labels):

        if elt.opts['dir'] != 'right':
            raise ValueError('Cannot draw opamp %s in direction %s'
                             % (elt.name, elt.opts['dir']))

        n1, n2, n3, n4 = elt.nodes

        p1, p2, p3, p4 = [self.coords[n] for n in elt.nodes]

        centre = Pos(0.5 * (p3.x + p1.x), p1.y)

        labelstr = elt.tex_label if draw_labels else ''
        argstr = '' if 'mirror' in elt.opts else 'yscale=-1'

        s = r'  \draw (%s) node[op amp, %s, scale=%.1f] (opamp) {};' % (
            centre, argstr, self.scale * 2)
        # Draw label separately to avoid being scaled by 2.
        s += r'  \draw (%s) node [] {%s};' % (centre, labelstr)
        return s

    def _tikz_draw_TF1(self, elt, nodes, draw_labels, link=False):

        p1, p2, p3, p4 = [self.coords[n] for n in nodes]

        xoffset = 0.06
        yoffset = 0.40

        primary_dot = Pos(p3.x - xoffset, 0.5 * (p3.y + p4.y) + yoffset)
        secondary_dot = Pos(p1.x + xoffset, 0.5 * (p1.y + p2.y) + yoffset)

        centre = Pos(0.5 * (p3.x + p1.x), 0.5 * (p2.y + p1.y))
        labelpos = Pos(centre.x, primary_dot.y)

        labelstr = elt.tex_label if draw_labels else ''

        s = r'  \draw (%s) node[circ] {};''\n' % primary_dot
        s += r'  \draw (%s) node[circ] {};''\n' % secondary_dot
        s += r'  \draw (%s) node[minimum width=%.1f] {%s};''\n' % (
            labelpos, 0.5, labelstr)

        if link:
            width = p1.x - p3.x
            arcpos = Pos((p1.x + p3.x) / 2, secondary_dot.y - width / 2 + 0.2)

            s += r'  \draw [<->] ([shift=(45:%.2f)]%s) arc(45:135:%.2f);' % (
                width / 2, arcpos, width / 2)
            s += '\n'

        return s

    def _tikz_draw_TF(self, elt, draw_labels):

        if elt.opts['dir'] != 'right':
            raise ValueError('Cannot draw transformer %s in direction %s'
                             % (elt.name, elt.opts['dir']))

        n1, n2, n3, n4 = elt.nodes

        s = r'  \draw (%s) to [inductor] (%s);''\n' % (n3, n4)
        s += r'  \draw (%s) to [inductor] (%s);''\n' % (n1, n2)
        s += self._tikz_draw_TF1(elt, elt.nodes, draw_labels)
        return s

    def _tikz_draw_TP(self, elt, draw_labels):

        if elt.opts['dir'] != 'right':
            raise ValueError('Cannot draw twoport network %s in direction %s'
                             % (elt.name, elt.opts['dir']))

        p1, p2, p3, p4 = [self.coords[n] for n in elt.nodes]
        width = p2.x - p4.x
        # height = p1.y - p2.y
        extra = 0.25
        p1.y += extra
        p2.y -= extra
        p3.y += extra
        p4.y -= extra
        centre = Pos(0.5 * (p3.x + p1.x), 0.5 * (p2.y + p1.y))
        top = Pos(centre.x, p1.y + 0.15)

        labelstr = elt.tex_label if draw_labels else ''
        titlestr = "%s-parameter two-port" % elt.args[2]

        s = r'  \draw (%s) -- (%s) -- (%s) -- (%s) -- (%s);''\n' % (
            p4, p3, p1, p2, p4)
        s += r'  \draw (%s) node[minimum width=%.1f] {%s};''\n' % (
            centre, width, titlestr)
        s += r'  \draw (%s) node[minimum width=%.1f] {%s};''\n' % (
            top, width, labelstr)
        return s

    def _tikz_draw_K(self, elt, draw_labels):

        if elt.opts['dir'] != 'right':
            raise ValueError('Cannot draw mutual coupling %s in direction %s'
                             % (elt.name, elt.opts['dir']))

        L1 = self.elements[elt.nodes[0]]
        L2 = self.elements[elt.nodes[1]]

        nodes = L2.nodes + L1.nodes

        s = self._tikz_draw_TF1(elt, nodes, draw_labels, link=True)
        return s

    def _tikz_draw_Q(self, elt, draw_labels):

        # For common base, will need to support up and down.
        if elt.opts['dir'] not in ('left', 'right'):
            raise ValueError('Cannot draw transistor %s in direction %s'
                             '; try left or right'
                             % (elt.name, elt.opts['dir']))

        n1, n2, n3 = elt.nodes

        p1, p2, p3 = [self.coords[n] for n in elt.nodes]

        centre = Pos(p3.x, 0.5 * (p1.y + p3.y))

        labelstr = elt.tex_label if draw_labels else ''
        sub_type = elt.sub_type.replace('jf', 'jfet')
        argstr = '' if elt.opts['dir'] == 'right' else 'xscale=-1'
        if 'mirror' in elt.opts:
            argstr += ', yscale=-1'

        s = r'  \draw (%s) node[%s, %s, scale=%.1f] () {};' % (
            centre, sub_type, argstr, self.scale * 2)
        s += r'  \draw (%s) node [] {%s};' % (centre, labelstr)
        return s

    def _tikz_draw_cpt(self, elt, draw_labels, draw_nodes):

        # Mapping of component names to circuitikz names.
        cpt_type_map = {'R': 'R', 'C': 'C', 'L': 'L',
                        'Vac': 'sV', 'Vdc': 'V', 'Iac': 'sI', 'Idc': 'I',
                        'Vacstep': 'sV', 'Vstep': 'V',
                        'Iacstep': 'sI', 'Istep': 'I',
                        'Vimpulse': 'V', 'Iimpulse': 'I',
                        'Vs': 'V', 'Is': 'I',
                        'V': 'V', 'I': 'I', 'v': 'V', 'i': 'I',
                        'P': 'open', 'W': 'short',
                        'TF': 'transformer', 'D' : 'D',
                        'Z': 'Z', 'Y': 'Y'}

        cpt_type = cpt_type_map[elt.cpt_type]
        if cpt_type == 'R' and 'variable' in elt.opts:
            cpt_type = 'vR'

        n1, n2 = elt.nodes[0:2]

        modifier = ''
        if cpt_type in ('V', 'Vdc', 'I', 'Idc'):

            # circuitikz expects the positive node first, except for
            # voltage and current sources!  So swap the nodes
            # otherwise they are drawn the wrong way around.
            n1, n2 = n2, n1

            if elt.opts['dir'] in ('down', 'left'):
                # Draw label on RHS for vertical cpt and below
                # for horizontal cpt.
                modifier = '_'
        else:
            if elt.opts['dir'] in ('up', 'right'):
                # Draw label on RHS for vertical cpt and below
                # for horizontal cpt.
                modifier = '_'

        # Current, voltage, label options.
        # It might be better to allow any options and prune out
        # dir and size.
        opts_str = ''
        for opt in ('i', 'i_', 'i^', 'i_>', 'i_<', 'i^>', 'i^<',
                    'i>_', 'i<_', 'i>^', 'i<^',
                    'v', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<',
                    'l', 'l^', 'l_'):
            if opt in elt.opts:
                opts_str += '%s=$%s$, ' % (opt, elt.opts[opt])

        node_str = self._node_str(n1, n2, draw_nodes)

        label_str = ''
        keys = elt.opts.keys()
        if draw_labels and not ('l' in keys or 'l_' in keys or 'l^' in keys):
            if cpt_type not in ('open', 'short'):
                label_str = ', l%s=%s' % (modifier, elt.tex_label)

        if cpt_type in ('Y', 'Z'):
            cpt_type = 'european resistor'
        elif cpt_type == 'D':
            diode_type_map = {'' : 'D', 'led' : 'leD', 'tunnel' : 'tD',
                              'zener' : 'zD', 'photo' : 'pD',
                              'schottky' : 'sD'}
            cpt_type = diode_type_map[elt.sub_type]

        s = r'  \draw (%s) to [%s%s, %s%s] (%s);''\n' % (
            n1, cpt_type, label_str, opts_str, node_str, n2)
        return s

    def _tikz_draw(self, draw_labels=True, draw_nodes=True,
                   label_nodes=True, args=None):

        # Preamble
        if args is None:
            args = ''

        opts = r'scale=%.2f,/tikz/circuitikz/bipoles/length=%.1fcm,%s' % (
            self.node_spacing, self.cpt_size, args)
        s = r'\begin{tikzpicture}[%s]''\n' % opts

        # Write coordinates
        for coord in self.coords.keys():
            s += r'  \coordinate (%s) at (%s);''\n' % (
                coord, self.coords[coord])

        draw = {'TF': self._tikz_draw_TF,
                'TP': self._tikz_draw_TP,
                'K': self._tikz_draw_K,
                'J': self._tikz_draw_Q,
                'M': self._tikz_draw_Q,
                'Q': self._tikz_draw_Q,
                'opamp': self._tikz_draw_opamp}

        # Draw components
        for m, elt in enumerate(self.elements.values()):

            if elt.cpt_type in draw:
                s += draw[elt.cpt_type](elt, draw_labels)
            else:
                s += self._tikz_draw_cpt(elt, draw_labels, draw_nodes)

        wires = self._make_wires()

        if False:
            # Draw implict wires
            for wire in wires:
                n1, n2 = wire.nodes

                node_str = self._node_str(n1, n2, draw_nodes)
                s += r'  \draw (%s) to [short, %s] (%s);''\n' % (
                    n1, node_str, n2)

        # Label primary nodes
        if label_nodes:
            for m, node in enumerate(self.nodes.values()):
                if not node.primary:
                    continue
                s += r'  \draw {[anchor=south east] (%s) node {%s}};''\n' % (
                    node.name, node.name)

        s += r'\end{tikzpicture}''\n'

        return s

    def _tmpfilename(self, suffix=''):

        from tempfile import NamedTemporaryFile

        filename = NamedTemporaryFile(suffix=suffix, delete=False).name
        return filename

    def tikz_draw(self, filename, args, **kwargs):

        root, ext = path.splitext(filename)

        content = self._tikz_draw(args=args, **kwargs)

        if ext == '.pytex':
            open(filename, 'w').write(content)
            return

        template = ('\\documentclass[a4paper]{standalone}\n'
                    '\\usepackage[americanvoltages]{circuitikz}\n'
                    '\\begin{document}\n%s\\end{document}')
        content = template % content

        texfilename = filename.replace(ext, '.tex')
        open(texfilename, 'w').write(content)

        if ext == '.tex':
            return

        dirname = path.dirname(texfilename)
        baseroot = path.basename(root)
        chdir = '' if dirname == '' else 'cd %s; ' % dirname

        system('%spdflatex -interaction batchmode %s.tex' % (chdir, baseroot))
        system('rm -f %s.aux %s.log %s.tex' % (root, root, root))

        pdf_filename = root + '.pdf'
        if not path.exists(pdf_filename):
            raise RuntimeError('Could not generate %s with pdflatex' % 
                               pdf_filename)

        if ext == '.pdf':
            return

        if ext == '.svg':
            svg_filename = root + '.svg'
            system('pdf2svg %s %s' % (pdf_filename, svg_filename))
            if not path.exists(svg_filename):
                raise RuntimeError('Could not generate %s with pdf2svg' % 
                                   svg_filename)
            remove(pdf_filename)
            return

        if ext == '.png':
            png_filename = root + '.png'
            system('convert -density 200 %s %s' %
                   (pdf_filename, png_filename))
            if not path.exists(png_filename):
                raise RuntimeError('Could not generate %s with convert' % 
                                   png_filename)
            remove(pdf_filename)
            return

        raise ValueError('Cannot create file of type %s' % ext)

    def draw(self, filename=None, args=None, stretch=1, scale=1, **kwargs):

        self.node_spacing = 2 * stretch * scale
        self.cpt_size = 1.5 * scale
        self.scale = scale

        def in_ipynb():
            try:
                ip = get_ipython()
                cfg = ip.config
                if cfg['IPKernelApp']['parent_appname'] == 'ipython-notebook':
                    return True
                else:
                    return False
            except (NameError, KeyError):
                return False

        if not self.hints:
            raise RuntimeWarning('No schematic drawing hints provided!')

        png = 'png' in kwargs and kwargs.pop('png')
        svg = 'svg' in kwargs and kwargs.pop('svg')

        if not png and not svg:
            png = True

        if in_ipynb() and filename is None:

            if png:
                from IPython.display import Image

                pngfilename = self._tmpfilename('.png')
                self.tikz_draw(pngfilename, args=args, **kwargs)

                # Display image.
                result = Image(pngfilename)
                return result

            if svg:
                from IPython.display import SVG

                svgfilename = self._tmpfilename('.svg')
                self.tikz_draw(svgfilename, args=args, **kwargs)

                # Display image.  There is a problem displaying
                # multiple SVG files since the later ones inherit
                # the namespace of the first ones.
                result = SVG(svgfilename)
                return result

        display = False
        if filename is None:
            filename = self._tmpfilename('.png')
            display = True

        self.tikz_draw(filename=filename, args=args, **kwargs)
        
        if display:
            # TODO display as SVG...

            from matplotlib.pyplot import figure
            from matplotlib.image import imread

            img = imread(filename)

            fig = figure()
            ax = fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('equal')
            ax.axis('off')
        return None


def test():

    sch = Schematic()

    sch.add('P1 1 0.1')
    sch.add('R1 1 3; right')
    sch.add('L1 3 2; right')
    sch.add('C1 3 0; down')
    sch.add('P2 2 0.2')
    sch.add('W 0.1 0; right')
    sch.add('W 0 0.2; right')

    sch.draw()
    return sch
