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

Copyright 2014 Michael Hayes, UCECE
"""

from __future__ import print_function
import numpy as np
import re
import sympy as sym
from lcapy.core import Expr

__all__ = ('Schematic', )


# Mapping of component names to circuitikz names.   The keys define
# the allowable component names.
cpt_type_map = {'R' : 'R', 'C' : 'C', 'L' : 'L', 'G' : 'G',
                'Vac' : 'sV', 'Vdc' : 'V', 'Iac' : 'sI', 'Idc' : 'I', 
                'Vacstep' : 'sV', 'Vstep' : 'V', 'Iacstep' : 'sI', 'Istep' : 'I', 
                'Vimpulse' : 'V', 'Iimpulse' : 'I',
                'Vs' : 'V', 'Is' : 'I',
                'V' : 'V', 'I' : 'I', 'v' : 'V', 'i' : 'I',
                'P' : 'open', 'W' : 'short', 
                'TF' : 'transformer', 
                'Z' : 'Z', 'Y' : 'Y'}


# Regular expression alternate matches stop with first match so need
# to have longer names first.
cpt_types = ['R', 'C', 'L', 'G', 'Z', 'Y', 'V', 'I', 'W', 'P', 'E', 'TF']
cpt_types.sort(lambda x, y: cmp(len(y), len(x)))

cpt_type_pattern = re.compile(r'(%s)(\w)?' % '|'.join(cpt_types))


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
            string = string.rstrip('0').rstrip('.')
            
        return string + '\,' + prefixes[idx] + self.unit


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

        return ', '.join(['%s=%s' % (key, val) for key, val in self.iteritems()])


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
        map = {}
        for node, nodes in self.sets.iteritems():
            map[node] = unique.index(nodes)
        
        self._map = map
        self._nodes = unique


    @property
    def map(self):
        """Return mapping of node number to common node number"""

        if not hasattr(self, '_map'):
            self._analyse()

        return self._map


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


def longest_path(all_nodes, from_nodes):
    """Find longest path through DAG.  all_nodes is an iterable for all
    the nodes in the graph, from_nodes is a directory indexed by node
    that stores a tuple of tuples.  The first tuple element is the
    parent node and the second element is the minimium size of the
    component connecting the nodes.
    """

    memo = {}

    def get_longest(to_node):

        if to_node in memo:
            return memo[to_node]

        best = 0
        for from_node, size in from_nodes[to_node]:
            best = max(best, get_longest(from_node) + size)

        memo[to_node] = best

        return best

    length, node = max([(get_longest(to_node), to_node) for to_node in all_nodes])
    return length, node, memo


class Node(object):

    def __init__(self, name):

        self.name = name
        self.pos = None
        self.port = False
        parts = name.split('_')
        self.rootname = parts[0] if name[0] != '_' else name
        self.primary = len(parts) == 1
        self.list = []

    
    def append(self, elt):

        if elt.cpt_type in ('P', ):
            self.port = True

        self.list.append(elt)


class NetElement(object):

    cpt_type_counter = 0

    def __init__(self, name, node1, node2, *args, **opts):

        match = cpt_type_pattern.match(name)

        if not match:
            raise ValueError('Unknown schematic component %s' % name)

        # Circuitikz does not like a . in a name
        if node1.find('.') != -1:
            raise ValueError('Cannot have . in node name %s' % node1)
        if node2.find('.') != -1:
            raise ValueError('Cannot have . in node name %s' % node2)

        cpt_type = match.groups()[0]
        id = match.groups()[1]

        if id is None:
            NetElement.cpt_type_counter += 1
            id = '#%d' % NetElement.cpt_type_counter
            name = cpt_type + id

        cpt_type_orig = cpt_type
        if args != ():
            if cpt_type in ('V', 'I') and args[0] in ('ac', 'dc', 'step', 'acstep', 'impulse', 's'):
                cpt_type = cpt_type + args[0]
                args = args[1:]

        # Tuple of nodes
        self.nodes = (node1, node2)
        # Identifier for component, e.g., 'R1'
        self.name = name
        # Type of component, e.g., 'V'
        self.cpt_type = cpt_type
        # Component arguments
        self.args = args

        if cpt_type in ('E', 'TF'):
            if len(args) < 3:
                raise ValueError('Component type %s requires 3 args' % cpt_type)
            self.nodes += (args[0], args[1])
            args = args[2:]

        symbol = None
        self.symbol = symbol

        autolabel = symbol
        if autolabel is None:
            autolabel = cpt_type_orig + '_{' + id + '}'

        if cpt_type in ('P', 'W') or autolabel.find('#') != -1:
            autolabel = ''

        if not opts.has_key('dir'):
            opts['dir'] = None
        if not opts.has_key('size'):
            opts['size'] = 1

        if opts['dir'] is None:
            opts['dir'] = 'down' if cpt_type in ('P', ) else 'right'

        if len(args) > 0:

            # TODO, extend for mechanical and acoustical components.
            units_map = {'V' : 'V', 'I' : 'A', 'R' : '$\Omega$',
                         'G' : 'S', 'C' : 'F', 'L' : 'H'}

            expr = args[0]
            if cpt_type in ('Vimpulse', 'Iimpulse'):
                expr = '(%s) * DiracDelta(t)' % expr
                autolabel = Expr(expr).latex()
            elif cpt_type in ('Vstep', 'Istep'):
                expr = '(%s) * Heaviside(t)' % expr
                autolabel = Expr(expr).latex()
            elif cpt_type in ('Vs', 'Is'):
                autolabel = Expr(expr).latex()
            else:
                try:
                    value = float(args[0])
                    if cpt_type[0] in units_map:
                        autolabel = EngFormat(value, units_map[cpt_type[0]]).latex()

                except ValueError:
                    autolabel = Expr(expr).latex()

        # Default label to use when drawing
        self.autolabel = autolabel
        # Drawing hints
        self.opts = Opts(opts)


    def __repr__(self):

        str = ', '.join(arg.__str__() for arg in [self.name] + list(self.nodes))
        return 'NetElement(%s)' % str


    def __str__(self):

        return ' '.join(['%s' % arg for arg in (self.name, ) + self.nodes])


class Schematic(object):

    def __init__(self, filename=None):

        self.elements = {}
        self.nodes = {}
        # Shared nodes (with same voltage)
        self.snodes = {}
        self.scale = 2
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

        if not self.nodes.has_key(node):
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

            if elt.cpt_type == 'TF' and dirs[0] == 'right':
                cnodes.link(*elt.nodes[0:2])
                cnodes.link(*elt.nodes[2:4])
                continue

            if elt.opts['dir'] in dirs:
                continue

            cnodes.link(*elt.nodes[0:2])

        # Now form forward and reverse directed graphs using components
        # in the desired directions.
        graph = {}
        rgraph = {}
        for m in range(cnodes.size):
            graph[m] = []
            rgraph[m] = []

        for m, elt in enumerate(self.elements.values()):
            if elt.opts['dir'] not in dirs:
                continue

            m1 = cnodes.map[elt.nodes[0]]
            m2 = cnodes.map[elt.nodes[1]]

            size = float(elt.opts['size'])

            if elt.opts['dir'] == dirs[0]:
                graph[m1].append((m2, size))
                rgraph[m2].append((m1, size))
            elif elt.opts['dir'] == dirs[1]:
                graph[m2].append((m1, size))
                rgraph[m1].append((m2, size))

        # Chain all potential start nodes to node 0.
        orphans = []
        rorphans = []
        for m in range(1, cnodes.size):
            if graph[m] == []:
                orphans.append((m, 0))
            if rgraph[m] == []:
                rorphans.append((m, 0))
        graph[0] = rorphans
        rgraph[0] = orphans

        if False:
            print(graph)
            print(rgraph)
            print(cnodes.map)
            import pdb
            pdb.set_trace()


        # Find longest path through the graphs.
        length, node, memo = longest_path(graph.keys(), graph)
        length, node, memor = longest_path(graph.keys(), rgraph)

        pos = {}
        posr = {}
        posa = {}
        for cnode in graph.keys():
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
            coords[node] = np.array((xpos[node], ypos[node]))

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

        if self.nodes[n1].port:
            node_str = 'o'
        else:
            node_str = '*' if draw_nodes and n1.find('_') == - 1 else ''
            
        node_str += '-'

        if self.nodes[n2].port:
            node_str += 'o'
        else:
            node_str += '*' if draw_nodes and n2.find('_') == - 1 else ''

        if node_str == '-':
            node_str = ''
        
        return node_str


    def tikz_draw(self, draw_labels=True, draw_nodes=True, label_nodes=True,
                  filename=None, args=None):

        if filename != None and filename != '':
            outfile = open(filename, 'w')
        else:
            import sys
            outfile = sys.stdout

        # Preamble
        if args is None: args = ''
        print(r'\begin{tikzpicture}[%s]' % args, file=outfile)

        # Write coordinates
        for coord in self.coords.keys():
            print(r'    \coordinate (%s) at (%.1f, %.1f);' % (coord, self.coords[coord][0] * self.scale, self.coords[coord][1] * self.scale), file=outfile)


        # Draw components
        for m, elt in enumerate(self.elements.values()):

            n1, n2 = elt.nodes[0:2]

            cpt_type = cpt_type_map[elt.cpt_type]

            # circuittikz expects the positive node first, except for 
            # voltage and current sources!   So swap the nodes otherwise
            # they are drawn the wrong way around.
            modifier = ''
            if elt.opts['dir'] == 'down' and cpt_type in ('V', 'Vdc', 'I', 'Idc'):
                n1, n2 = n2, n1
                # Draw label on RHS for vertical cpt.
                modifier = '_'
      
            # If have a left drawn cpt, then switch nodes so that
            # label defaults to top but then have to switch current
            # and voltage directions.
            if elt.opts['dir'] == 'left':
                n1, n2 = n2, n1
                if elt.opts.has_key('i'):
                    elt.opts['i<^'] = elt.opts.pop('i')
                if elt.opts.has_key('v'):
                    elt.opts['v_>'] = elt.opts.pop('v')

            # Current, voltage, label options.
            # It might be better to allow any options and prune out
            # dir and size.
            opts_str = ''
            for opt in ('i', 'i_', 'i^', 'i_>', 'i_<', 'i^>', 'i^<', 
                        'i>_', 'i<_', 'i>^', 'i<^', 
                        'v', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<', 'l', 'l^', 'l_'):
                if elt.opts.has_key(opt):
                    opts_str += '%s=$%s$, ' % (opt, elt.opts[opt])

            node_str = self._node_str(n1, n2, draw_nodes)
               
            label_str =''
            if draw_labels and not ('l' in elt.opts.keys() or 'l_' in elt.opts.keys() or 'l^' in elt.opts.keys()):
                if cpt_type not in ('open', 'short'):
                    label_parts = elt.autolabel.split('\\,')
                    if len(label_parts) > 1:
                        label_str = ', l%s=$%s\mbox{%s}$' % (modifier, label_parts[0], label_parts[1])
                    else:
                        label_str = ', l%s=$%s$' % (modifier, label_parts[0])

            if cpt_type in ('Y', 'Z'):
                cpt_type = 'european resistor'

            print(r'    \draw (%s) to [%s%s, %s%s] (%s);' % (n1, cpt_type, label_str, opts_str, node_str, n2), file=outfile)

        wires = self._make_wires()

        # Draw implict wires
        for wire in wires:
            n1, n2 = wire.nodes

            node_str = self._node_str(n1, n2, draw_nodes)
            print(r'    \draw (%s) to [short, %s] (%s);' % (n1, node_str, n2),
                  file=outfile)
    
        # Label primary nodes
        if label_nodes:
            for m, node in enumerate(self.nodes.values()):
                if not node.primary:
                    continue
                print(r'    \draw {[anchor=south east] (%s) node {%s}};' % (node.name, node.name), file=outfile)

        print(r'\end{tikzpicture}', file=outfile)


    def schemdraw_draw(self, draw_labels=True, draw_nodes=True, 
                       label_nodes=True, filename=None, args=None):

        from SchemDraw import Drawing
        import SchemDraw.elements as e

        cpt_type_map2 = {'R' : e.RES, 'C' : e.CAP, 'L' : e.INDUCTOR2, 
                         'Vac' : e.SOURCE_SIN, 'Vdc' : e.SOURCE_V,
                         'Iac' : e.SOURCE_SIN, 'Idc' : e.SOURCE_I, 
                         'V' : e.SOURCE_V, 'I' : e.SOURCE_I, 
                         'Vs' : e.SOURCE_V, 'Is' : e.SOURCE_I, 
                         'v' : e.SOURCE_V, 'i' : e.SOURCE_I,
                         'P' : e.GAP_LABEL, 'port' : e.GAP_LABEL,
                         'W' : e.LINE, 'wire' : e.LINE,
                         'Y' : e.RBOX, 'Z' : e.RBOX}        

        # Preamble
        if args is None: args = ''
        
        drw = Drawing()

        # Update element positions
        self.coords

        # Draw components
        for m, elt in enumerate(self.elements.values()):

            # FIXME: transformer needs to be drawn as a pair of inductors
            cpt_type = cpt_type_map2[elt.cpt_type]

            n1, n2 = elt.nodes[0:2]

            pos1 = self.coords[n1]
            pos2 = self.coords[n2]

            if draw_labels:
                drw.add(cpt_type, xy=pos2 * self.scale, 
                        to=pos1 * self.scale, 
                        label=elt.autolabel.replace('\\,', ' '))
            else:
                drw.add(cpt_type, xy=pos2 * self.scale,
                        to=pos1 * self.scale)

        if draw_nodes:
            for m, node in enumerate(self.nodes.values()):
                label_str = node.name if draw_labels and node.primary else ''
                if node.port:
                    drw.add(e.DOT_OPEN, xy=self.coords[node.name] * self.scale,
                            label=label_str)
                elif node.primary:
                    drw.add(e.DOT, xy=self.coords[node.name] * self.scale, 
                            label=label_str)

        drw.draw()
        if filename is not None:
            drw.save(filename)


    def draw(self, draw_labels=True, draw_nodes=True, label_nodes=True,
             filename=None, args=None, scale=1, tex=False):

        self.scale = scale * 2

        if not self.hints:
            raise RuntimeWarning('No schematic drawing hints provided!')

        if tex or (filename is not None and filename.endswith('.tex')):
            self.tikz_draw(draw_labels=draw_labels, draw_nodes=draw_nodes,
                           label_nodes=label_nodes, filename=filename,
                           args=args)            
        else:
            self.schemdraw_draw(draw_labels=draw_labels, draw_nodes=draw_nodes, 
                                label_nodes=label_nodes, filename=filename)


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
