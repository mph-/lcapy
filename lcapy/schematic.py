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
from lcapy.latex import latex_str
from lcapy.core import Expr
from os import system, path, remove, mkdir, chdir


__all__ = ('Schematic', )

# Perhaps AM for ammeter, VM for voltmeter, SW for switch, VR for
# variable resistor?  Currently, a variable resistor is supported with
# the variable option.


# Regular expression alternate matches stop with first match so need
# to have longer names first.
cpt_types = ['R', 'C', 'L', 'Z', 'Y', 'V', 'I', 'W', 'O', 'P', 'E', 
             'D', 'J', 'M', 'Q', 'SW', 'TF', 'TP', 'K']
cpt_types.sort(lambda x, y: cmp(len(y), len(x)))

cpt_type_pattern = re.compile(r"(%s)([\w']*)" % '|'.join(cpt_types))

label_pattern = re.compile(r"([\w']*)(_{[\w]})?")


import math


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
            arg = '='.join(fields[1:]).strip() if len(fields) > 1 else ''
            self[key] = arg

    def __init__(self, arg=None):

        if arg is None:
            return

        if isinstance(arg, str):
            self._parse(arg)
            return

        for key, val in arg.iteritems():
            self[key] = val

    @property
    def size(self):

        size = self.get('size', 1)
        return float(size)

    def format(self):

        return ', '.join(['%s=%s' % (key, val)
                          for key, val in self.iteritems()])

    def copy(self):
        
        return self.__class__(super(Opts, self).copy())

    def strip(self, *args):

        stripped = Opts()
        for opt in args:
            if opt in self:
                stripped[opt] = self.pop(opt)        
        return stripped

    def strip_voltage_labels(self):

        return self.strip('v', 'vr', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<')

    def strip_current_labels(self):

        return self.strip('i', 'ir', 'i_', 'i^', 'i_>', 'i_<', 'i^>', 'i^<',
                          'i>_', 'i<_', 'i>^', 'i<^')

    def strip_labels(self):

        return self.strip('l', 'l^', 'l_')

    def strip_all_labels(self):

        self.strip_voltage_labels()
        self.strip_current_labels()
        self.strip_labels()


class EngFormat(object):

    def __init__(self, value, unit=''):

        self.value = value
        self.unit = unit

    def latex(self, trim=True, hundreds=False):

        prefixes = ('f', 'p', 'n', '$\mu$', 'm', '', 'k', 'M', 'G', 'T')

        sfmax = 3

        value = self.value
        m = math.log10(abs(value))

        if m < -1 or m >= 3.0:
            if hundreds:
                # Generate 100 m
                n = int(math.floor(m / 3))
                k = int(math.floor(m)) - n * 3
            else:
                # Generate 0.1
                n = int(round(m / 3))
                k = int(round(m)) - n * 3
        else:
            n = 0
            k = m - 1

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

        return string + '\,' + r'\mbox{' + prefixes[idx] + self.unit + r'}'


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

        # Add dummy cnode at start
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
            nodes = list([nodes])

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

        if self[0] == []:
            raise ValueError("Cannot find start node for graph '%s'. "
                             "Probably a component has an incorrect direction."
                             % self.name)

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
                ("The schematic graph '%s' is dodgy, probably a component"
                 " is connected to the wrong node\n%s") 
                % (self.name, from_nodes))

        return length, node, memo


class Graphs(object):

    def __init__(self, size, name):

        self.fwd = Graph(size, 'forward ' + name)
        self.rev = Graph(size, 'reverse ' + name)
        self.size = size

    def add(self, n1, n2, size):
        self.fwd.add(n1, n2, size)
        self.rev.add(n2, n1, size)

    @property
    def nodes(self):
        return self.fwd.keys()

    def add_start_nodes(self):

        # Chain all potential start nodes to node 0.
        orphans = []
        rorphans = []
        for m in range(1, self.size):
            if self.fwd[m] == []:
                orphans.append((m, 0))
            if self.rev[m] == []:
                rorphans.append((m, 0))
        self.fwd[0] = rorphans
        self.rev[0] = orphans


class Node(object):

    def __init__(self, name):

        self.name = name
        self._port = False
        self._count = 0
        parts = name.split('_')
        self.rootname = parts[0] if name[0] != '_' else name
        self.primary = len(parts) == 1
        self.list = []

    def append(self, elt):
        """Add new element to the node"""

        if elt.cpt_type == 'P':
            self._port = True

        self.list.append(elt)
        if elt.cpt_type not in ('O', ):
            self._count += 1

    @property
    def count(self):
        """Number of elements (including wires but not open-circuits)
        connected to the node"""

        return self._count

    def visible(self, draw_nodes):
        """Return true if node drawn"""

        if self.port:
            return True

        if draw_nodes in ('none', None, False):
            return False
        
        if draw_nodes == 'all':
            return True

        if draw_nodes == 'connections':
            return self.count > 2

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

        def tex_name(name, subscript=''):

            if len(name) > 1:
                name = r'\mathrm{%s}' % name
            if len(subscript) > 1:
                subscript = r'\mathrm{%s}' % subscript
            if len(subscript) == 0:
                return name
        
            return '%s_{%s}' % (name, subscript)

        match = cpt_type_pattern.match(name)

        if not match:
            raise ValueError('Unknown schematic component %s' % name)

        # Circuitikz does not like a . in a name
        if n1.find('.') != -1:
            raise ValueError('Cannot have . in node name %s' % n1)
        if n2.find('.') != -1:
            raise ValueError('Cannot have . in node name %s' % n2)

        cpt_type = match.groups()[0]
        cpt_id = match.groups()[1]

        if cpt_id is None:
            NetElement.cpt_type_counter += 1
            name = cpt_type + '#%d' % NetElement.cpt_type_counter

        # There are two possible labels for a component:
        # 1. Component identifier, e.g., R1
        # 2. Component value, expression, or symbol
        id_label = tex_name(cpt_type, cpt_id)
        value_label = None

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
            elif cpt_type == 'SW' and args[0] in ('nc', 'no', 'push'):
                cpt_type = cpt_type + args[0]
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

        if cpt_type in ('O', 'P', 'W') or id_label.find('#') != -1:
            id_label = None

        if 'dir' not in opts:
            opts['dir'] = None
        if 'size' not in opts:
            opts['size'] = 1

        if opts['dir'] is None:
            opts['dir'] = 'down' if cpt_type in ('O', 'P') else 'right'

        if len(args) > 0:

            # TODO, extend for mechanical and acoustical components.
            units_map = {'V': 'V', 'I': 'A', 'R': '$\Omega$',
                         'C': 'F', 'L': 'H'}

            expr = args[0]
            if cpt_type in ('Vimpulse', 'Iimpulse'):
                expr = '(%s) * DiracDelta(t)' % expr
                value_label = Expr(expr, cache=False).latex()
            elif cpt_type in ('Vstep', 'Istep'):
                expr = '(%s) * Heaviside(t)' % expr
                value_label = Expr(expr, cache=False).latex()
            elif cpt_type in ('Vs', 'Is'):
                value_label = Expr(expr, cache=False).latex()
            elif cpt_type == 'TF':
                value_label = '1:%s' % args[0]
            elif cpt_type not in ('TP',):
                try:
                    value = float(args[0])
                    if cpt_type[0] in units_map:
                        value_label = EngFormat(
                            value, units_map[cpt_type[0]]).latex()
                    else:
                        value_label = Expr(expr, cache=False).latex()

                except ValueError:
                    value_label = Expr(expr, cache=False).latex()

        # Currently, we only annnotated the component with the value,
        # expression, or symbol.  If this is not specified, it
        # defaults to the component identifier.  Note, some objects
        # we do not want to label, such as wires and ports.

        self.id_label = '' if id_label is None else latex_str(id_label)
        self.value_label = '' if value_label is None else latex_str(value_label)
        self.default_label = self.id_label if self.value_label == '' else self.value_label

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
        try:
            return self.elements[name]
        except KeyError:
            raise AttributeError('Unknown component %s' % name)

    def netfile_add(self, filename):
        """Add the nets from file with specified filename"""

        file = open(filename, 'r')

        lines = file.readlines()

        for line in lines:
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

        nodes = elt.nodes
        # The controlling nodes are not drawn.
        if elt.cpt_type in ('E', 'F', 'G', 'H'):
            nodes = nodes[0:2]

        for node in nodes:
            self._node_add(node, elt)

    def add(self, string):
        """The general form is: 'Name Np Nm symbol'
        where Np is the positive nose and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        # Ignore comments
        string = string.strip()
        if string == '' or string[0] in ('#', '%'):
            return

        if '\n' in string:
            lines = string.split('\n')
            for line in lines:
                self.add(line)
            return

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
                xscale = 1.0 if elt.cpt_type == 'M' else 0.85
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

            # Handle bipoles here.  Treat VCCS etc as bipole.
            m1, m2 = cnodes.map(elt.nodes[0:2])

            if elt.opts['dir'] == dirs[0]:
                graphs.add(m1, m2, size)
            elif elt.opts['dir'] == dirs[1]:
                graphs.add(m2, m1, size)

        graphs.add_start_nodes()

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
        return posa, cnodes.nodes, length

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

        xpos, self._xnodes, self.width = self._make_graphs(('right', 'left'))
        ypos, self._ynodes, self.height = self._make_graphs(('up', 'down'))

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

        node_str = ''
        if node1.visible(draw_nodes):
            node_str = 'o' if node1.port else '*'

        node_str += '-'

        if node2.visible(draw_nodes):
            node_str += 'o' if node2.port else '*'

        if node_str == '-':
            node_str = ''
        
        return node_str


    def _tikz_draw_node(self, n, draw_nodes=True):
        
        s = ''
        if not draw_nodes:
            return s

        node = self.nodes[n]
        if not node.visible(draw_nodes):
            return s

        pos = self.coords[n]
        if node.port:
            s = r'  \draw (%s) node[ocirc] {};''\n' % pos
        else:
            s = r'  \draw (%s) node[circ] {};''\n' % pos

        return s

    def _tikz_draw_nodes(self, elt, draw_nodes=True):

        s = ''
        for n in elt.nodes:
            s += self._tikz_draw_node(n, draw_nodes)
        return s

    def _tikz_draw_opamp(self, elt, label_values, draw_nodes):

        if elt.opts['dir'] != 'right':
            raise ValueError('Cannot draw opamp %s in direction %s'
                             % (elt.name, elt.opts['dir']))

        n1, n2, n3, n4 = elt.nodes

        p1, p2, p3, p4 = [self.coords[n] for n in elt.nodes]

        centre = Pos(0.5 * (p3.x + p1.x), p1.y)

        label_str = '$%s$' % elt.default_label if label_values else ''
        args_str = '' if 'mirror' in elt.opts else 'yscale=-1'
        for key, val in elt.opts.iteritems():
            if key in ('color', ):
                args_str += '%s=%s, ' % (key, val)                

        s = r'  \draw (%s) node[op amp, %s, scale=%.1f] (opamp) {};' % (
            centre, args_str, self.scale * 2)
        # Draw label separately to avoid being scaled by 2.
        s += r'  \draw (%s) node [] {%s};' % (centre, label_str)
        
        s += self._tikz_draw_nodes(elt, draw_nodes)
        return s

    def _tikz_draw_TF1(self, elt, nodes, label_values, link=False):

        p1, p2, p3, p4 = [self.coords[n] for n in nodes]

        xoffset = 0.06
        yoffset = 0.40

        primary_dot = Pos(p3.x - xoffset, 0.5 * (p3.y + p4.y) + yoffset)
        secondary_dot = Pos(p1.x + xoffset, 0.5 * (p1.y + p2.y) + yoffset)

        centre = Pos(0.5 * (p3.x + p1.x), 0.5 * (p2.y + p1.y))
        labelpos = Pos(centre.x, primary_dot.y)

        label_str = '$%s$' % elt.default_label if label_values else ''

        s = r'  \draw (%s) node[circ] {};''\n' % primary_dot
        s += r'  \draw (%s) node[circ] {};''\n' % secondary_dot
        s += r'  \draw (%s) node[minimum width=%.1f] {%s};''\n' % (
            labelpos, 0.5, label_str)

        if link:
            width = p1.x - p3.x
            arcpos = Pos((p1.x + p3.x) / 2, secondary_dot.y - width / 2 + 0.2)

            s += r'  \draw [<->] ([shift=(45:%.2f)]%s) arc(45:135:%.2f);' % (
                width / 2, arcpos, width / 2)
            s += '\n'

        return s

    def _tikz_draw_TF(self, elt, label_values, draw_nodes):

        if elt.opts['dir'] != 'right':
            raise ValueError('Cannot draw transformer %s in direction %s'
                             % (elt.name, elt.opts['dir']))

        n1, n2, n3, n4 = elt.nodes

        s = r'  \draw (%s) to [inductor] (%s);''\n' % (n3, n4)
        s += r'  \draw (%s) to [inductor] (%s);''\n' % (n1, n2)
        s += self._tikz_draw_TF1(elt, elt.nodes, label_values)

        s += self._tikz_draw_nodes(elt, draw_nodes)
        return s

    def _tikz_draw_TP(self, elt, label_values, draw_nodes):

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

        label_str = '$%s$' % elt.default_label if label_values else ''
        titlestr = "%s-parameter two-port" % elt.args[2]

        s = r'  \draw (%s) -- (%s) -- (%s) -- (%s) -- (%s);''\n' % (
            p4, p3, p1, p2, p4)
        s += r'  \draw (%s) node[minimum width=%.1f] {%s};''\n' % (
            centre, width, titlestr)
        s += r'  \draw (%s) node[minimum width=%.1f] {%s};''\n' % (
            top, width, label_str)

        s += self._tikz_draw_nodes(elt, draw_nodes)
        return s

    def _tikz_draw_K(self, elt, label_values, draw_nodes):

        if elt.opts['dir'] != 'right':
            raise ValueError('Cannot draw mutual coupling %s in direction %s'
                             % (elt.name, elt.opts['dir']))

        L1 = self.elements[elt.nodes[0]]
        L2 = self.elements[elt.nodes[1]]

        nodes = L2.nodes + L1.nodes

        s = self._tikz_draw_TF1(elt, nodes, label_values, link=True)
        return s

    def _tikz_draw_Q(self, elt, label_values, draw_nodes):

        # For common base, will need to support up and down.
        if elt.opts['dir'] not in ('left', 'right'):
            raise ValueError('Cannot draw transistor %s in direction %s'
                             '; try left or right'
                             % (elt.name, elt.opts['dir']))

        n1, n2, n3 = elt.nodes

        p1, p2, p3 = [self.coords[n] for n in elt.nodes]

        centre = Pos(p3.x, 0.5 * (p1.y + p3.y))

        label_str = '$%s$' % elt.default_label if label_values else ''
        sub_type = elt.sub_type.replace('jf', 'jfet')
        args_str = '' if elt.opts['dir'] == 'right' else 'xscale=-1'
        if 'mirror' in elt.opts:
            args_str += ', yscale=-1'
        for key, val in elt.opts.iteritems():
            if key in ('color', ):
                args_str += '%s=%s, ' % (key, val)                

        s = r'  \draw (%s) node[%s, %s, scale=%.1f] (T) {};''\n' % (
            centre, sub_type, args_str, self.scale * 2)
        s += r'  \draw (%s) node [] {%s};''\n'% (centre, label_str)

        if sub_type in ('pnp', 'pmos', 'pjfet'):
            n1, n3 = n3, n1

        # Add additional wires.
        if elt.cpt_type in ('J', 'M'):
            s += r'  \draw (T.D) -- (%s) (T.G) -- (%s) (T.S) -- (%s);''\n' % (n1, n2, n3)
        else:
            s += r'  \draw (T.C) -- (%s) (T.B) -- (%s) (T.E) -- (%s);''\n' % (n1, n2, n3)

        s += self._tikz_draw_nodes(elt, draw_nodes)
        return s

    def _tikz_draw_cpt(self, elt, label_values, draw_nodes, label_ids):

        n1, n2 = elt.nodes[0:2]

        if elt.cpt_type == 'W' and ('implicit' in elt.opts
                                    or 'ground' in elt.opts
                                    or 'sground' in elt.opts):
                                    
            # Draw implict wires, i.e., connections to ground, etc.

            kind = ''
            if 'implicit' in elt.opts:
                kind = 'sground'
            if 'ground' in elt.opts:
                kind = 'ground'
            if 'sground' in elt.opts:
                kind = 'sground'

            args = [kind]
            if elt.opts['dir'] == 'up':
                args.append('yscale=-1')
            if elt.opts['dir'] == 'left':
                args.append('xscale=-1')
            args_str = ','.join(args)

            offset = 0.0
            anchor = 'south west'
            p = self.coords[n2]
            if elt.opts['dir'] == 'down':
                offset = 0.25
                anchor = 'north west'
            elif elt.opts['dir'] == 'up':
                offset = -0.25
            pos = Pos(p.x, p.y + offset)

            s = r'  \draw (%s) to [short] (%s);''\n' % (n1, pos)
            s += r'  \draw (%s) node[%s] {};''\n' % (pos, args_str)

            if 'l' in elt.opts:
                label_str = '${%s}$' % latex_str(elt.opts['l'])
                s += r'  \draw {[anchor=%s] (%s) node {%s}};''\n' % (anchor, n2, label_str)
            return s

        # Mapping of component names to circuitikz names.
        cpt_type_map = {'R': 'R', 'C': 'C', 'L': 'L',
                        'Vac': 'sV', 'Vdc': 'V', 'Iac': 'sI', 'Idc': 'I',
                        'Vacstep': 'sV', 'Vstep': 'V',
                        'Iacstep': 'sI', 'Istep': 'I',
                        'Vimpulse': 'V', 'Iimpulse': 'I',
                        'Vs': 'V', 'Is': 'I',
                        'V': 'V', 'I': 'I', 'v': 'V', 'i': 'I',
                        'O' : 'open', 'P': 'open', 'W': 'short',
                        'TF': 'transformer', 'D' : 'D',
                        'Z': 'Z', 'Y': 'Y',
                        'E' : 'american controlled voltage source',
                        'SWnc' : 'opening switch', 'SWno' : 'closing switch',
                        'SW' : 'closing switch', 'SWpush' : 'push button'}

        cpt_type = cpt_type_map[elt.cpt_type]
        if cpt_type == 'R' and 'variable' in elt.opts:
            cpt_type = 'vR'

        id_pos = '_'
        voltage_pos = '^'
        if elt.cpt_type in ('V', 'Vdc', 'Vstep', 'Vac', 'Vacstep', 'Vimpulse', 'v',
                            'I', 'Idc', 'Istep', 'Iac', 'Iacstep', 'Iimpulse', 'i',
                            'E', 'Vs', 'Is'):

            # circuitikz expects the positive node first, except for
            # voltage and current sources!  So swap the nodes
            # otherwise they are drawn the wrong way around.
            n1, n2 = n2, n1

            if elt.opts['dir'] in ('down', 'right'):
                # Draw label on LHS for vertical cpt and below
                # for horizontal cpt.
                id_pos = '^'
                voltage_pos = '_'
        else:
            if elt.opts['dir'] in ('up', 'left'):
                # Draw label on LHS for vertical cpt and below
                # for horizontal cpt.
                id_pos = '^'
                voltage_pos = '_'

        # Add modifier to place voltage label on other side
        # from component identifier label.
        if 'v' in elt.opts:
            elt.opts['v' + voltage_pos] = elt.opts.pop('v')

        # Reversed voltage.
        if 'vr' in elt.opts:
            elt.opts['v' + voltage_pos + '>'] = elt.opts.pop('vr')

        current_pos = id_pos
        # Add modifier to place current label on other side
        # from voltage marks.
        if 'i' in elt.opts:
            elt.opts['i' + current_pos] = elt.opts.pop('i')

        # Reversed current.
        if 'ir' in elt.opts:
            elt.opts['i' + current_pos + '<'] = elt.opts.pop('ir')

        # Current, voltage, label options.
        # It might be better to allow any options and prune out
        # dir and size.
        voltage_str = ''
        current_str = ''
        label_str = ''
        args_str = ''
        for key, val in elt.opts.iteritems():
            if key in ('i', 'i_', 'i^', 'i_>', 'i_<', 'i^>', 'i^<',
                       'i>_', 'i<_', 'i>^', 'i<^'):
                current_str += '%s=${%s}$, ' % (key, latex_str(val))
            elif key in ('v', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<'):
                voltage_str += '%s=${%s}$, ' % (key, latex_str(val))
            elif key in ('l', 'l^', 'l_'):
                label_str += '%s=${%s}$, ' % (key, latex_str(val))
            elif key in ('color', ):
                args_str += '%s=%s, ' % (key, val)                

        node_str = self._node_str(n1, n2, draw_nodes)

        args_str += voltage_str + current_str

        # Generate default label unless specified.
        if label_str == '':
            if cpt_type not in ('open', 'short'):
                
                label_str = ', l%s=${%s}$' % (id_pos, elt.default_label)
                
                if label_ids and elt.value_label != '':
                    label_str += r', l%s={$%s$\\${%s}$}' % (id_pos, elt.id_label, elt.value_label)
        else:
            label_str = ', ' + label_str

        if not label_values:
            label_str = ''

        if not label_values and label_ids:
            label_str = ', l%s=${%s}$' % (id_pos, elt.id_label)

        if cpt_type in ('Y', 'Z'):
            cpt_type = 'european resistor'
        elif cpt_type == 'D':
            diode_type_map = {'' : 'D', 'led' : 'leD', 'tunnel' : 'tD',
                              'zener' : 'zD', 'photo' : 'pD',
                              'schottky' : 'sD'}
            cpt_type = diode_type_map[elt.sub_type]

        s = r'  \draw (%s) to [align=right, %s%s, %s%s] (%s);''\n' % (
            n1, cpt_type, label_str, args_str, node_str, n2)
        return s

    def _tikz_draw(self, label_values=True, draw_nodes=True, label_ids=True,
                   label_nodes='primary', args=None):

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
                s += draw[elt.cpt_type](elt, label_values, draw_nodes)
            else:
                s += self._tikz_draw_cpt(elt, label_values, draw_nodes, label_ids)

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
                if label_nodes == 'primary' and not node.primary:
                    continue
                s += r'  \draw {[anchor=south east] (%s) node {%s}};''\n' % (
                    node.name, node.name.replace('_', r'\_'))

        s += r'\end{tikzpicture}''\n'

        return s

    def _tmpfilename(self, suffix=''):

        from tempfile import NamedTemporaryFile
        
        tempdir = 'temp'
        if not path.exists(tempdir):
            mkdir(tempdir)
        filename = NamedTemporaryFile(suffix=suffix, dir=tempdir, 
                                      delete=False).name
        return filename

    def _convert_pdf_svg(self, pdf_filename, svg_filename):

        system('pdf2svg %s %s' % (pdf_filename, svg_filename))
        if not path.exists(svg_filename):
            raise RuntimeError('Could not generate %s with pdf2svg' % 
                               svg_filename)

    def _convert_pdf_png(self, pdf_filename, png_filename, oversample=1):

        system('convert -density %d %s %s' %
               (oversample * 100, pdf_filename, png_filename))
        if path.exists(png_filename):
            return

        # Windows has a program called convert, try im-convert
        # for image magick convert.
        system('im-convert -density %d %s %s' %
               (oversample * 100, pdf_filename, png_filename))
        if path.exists(png_filename):
            return

        raise RuntimeError('Could not generate %s with convert' % 
                           png_filename)

    def tikz_draw(self, filename, args, **kwargs):

        root, ext = path.splitext(filename)

        debug = kwargs.pop('debug', False)
        oversample = float(kwargs.pop('oversample', 2))

        content = self._tikz_draw(args=args, **kwargs)

        if debug:
            print('width = %d, height = %d, oversample = %d'
                  % (self.width, self.height, oversample))

        if ext == '.pytex':
            open(filename, 'w').write(content)
            return

        template = ('\\documentclass[a4paper]{standalone}\n'
                    '\\usepackage[americanvoltages,americancurrents]{circuitikz}\n'
                    '\\begin{document}\n%s\\end{document}')
        content = template % content

        texfilename = filename.replace(ext, '.tex')
        open(texfilename, 'w').write(content)

        if ext == '.tex':
            return

        dirname = path.dirname(texfilename)
        baseroot = path.basename(root)
        if dirname != '':
            chdir(path.abspath(dirname))

        system('pdflatex -interaction batchmode %s.tex' % baseroot)
        if not debug:
            try:
                remove(baseroot + '.aux')
                remove(baseroot + '.log')
                remove(baseroot + '.tex')
            except:
                pass

        pdf_filename = root + '.pdf'
        if not path.exists(pdf_filename):
            raise RuntimeError('Could not generate %s with pdflatex' % 
                               pdf_filename)

        if ext == '.pdf':
            return

        if ext == '.svg':
            self._convert_pdf_svg(pdf_filename, root + '.svg')
            if not debug:
                remove(pdf_filename)
            return

        if ext == '.png':
            self._convert_pdf_png(pdf_filename, root + '.png', oversample)
            if not debug:
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

                kernapp = cfg['IPKernelApp']

                # Check if processing ipynb file.
                if 'connection_file' in kernapp:
                    return True
                elif kernapp and kernapp['parent_appname'] == 'ipython-notebook':
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
                from IPython.display import Image, display_png

                pngfilename = self._tmpfilename('.png')
                self.tikz_draw(pngfilename, args=args, **kwargs)

                # Create and display PNG image object.
                # There are two problems:
                # 1. The image metadata (width, height) is ignored
                #    when the ipynb file is loaded.
                # 2. The image metadata (width, height) is not stored
                #    when the ipynb file is written non-interactively.
                display_png(Image(filename=pngfilename,
                                  width=self.width * 100, 
                                  height=self.height * 100))
                return

            if svg:
                from IPython.display import SVG, display_svg

                svgfilename = self._tmpfilename('.svg')
                self.tikz_draw(svgfilename, args=args, **kwargs)

                # Create and display SVG image object.
                # Note, there is a problem displaying multiple SVG
                # files since the later ones inherit the namespace of
                # the first ones.
                display_svg(SVG(filename=pngfilename, 
                                width=self.width * 100, height=self.height * 100))
                return

        display = False
        if filename is None:
            filename = self._tmpfilename('.png')
            display = True

        self.tikz_draw(filename=filename, args=args, **kwargs)
        
        if display:
            # TODO display as SVG so have scaled fonts...

            from matplotlib.pyplot import figure
            from matplotlib.image import imread

            img = imread(filename)

            fig = figure()
            ax = fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('equal')
            ax.axis('off')

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
