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

Copyright 2014, 2015, 2016 Michael Hayes, UCECE
"""

# Components are positioned using two graphs; one graph for
# the x direction and the other for the y direction. 
#
# There is naming confusion.  We have network nodes (electrical
# nodes) and nodes in the graph used for component placement.
# Let's call the latter gnodes.  The names of these gnodes are a
# tuple of the common network nodes.
#
# x and y component positioning are performed independently.  Let's
# consider the x or horizontal positioning.  There are three stages:
#   1. Component nodes that share a y position are linked; this can
#      occur, for example, for a vertically oriented component.
#      This helps to reduce the size of the graph.
#   2. The x positions of the components are used to determine the
#      graph edges.
#   3. The longest path through the graph is found and the x positions
#      of the nodes are assigned based on the distance along the
#      longest path.


from __future__ import print_function
import numpy as np
import re
from lcapy.latex import latex_str, format_label
from lcapy.core import Expr
import lcapy.grammar as grammar
from lcapy.parser import Parser
import lcapy.schemcpts as cpts
from lcapy.schemmisc import Pos, Opts
from os import system, path, remove, mkdir, chdir, getcwd
import math

__all__ = ('Schematic', )

parser = Parser(cpts, grammar)

def tmpfilename(suffix=''):

    from tempfile import gettempdir, NamedTemporaryFile
    
    # Searches using TMPDIR, TEMP, TMP environment variables
    tempdir = gettempdir()
    
    filename = NamedTemporaryFile(suffix=suffix, dir=tempdir, 
                                  delete=False).name
    return filename


def display_matplotlib(filename):
        
    from matplotlib.pyplot import figure
    from matplotlib.image import imread
    
    img = imread(filename)
    
    fig = figure()
    ax = fig.add_subplot(111)
    ax.imshow(img)
    ax.axis('equal')
    ax.axis('off')


class SchematicOpts(Opts):

    def __init__(self):

        super (SchematicOpts, self).__init__(
            {'draw_nodes': 'primary',
             'label_values': True,
             'label_ids': True,
             'label_nodes': 'primary',
             'scale' : 1.0,
             'cpt_size' : 1.5,
             'node_spacing' : 2.0,
             'append' : '',
             'help_lines' : 0.0,
             'style' : 'american'})


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


class Cnodes(dict):
    """Common nodes"""

    def __init__(self, nodes):

        super (Cnodes, self).__init__()
        for node in nodes:
            self[node] = (node, )

    def link(self, n1, n2):
        """Make nodes n1 and n2 share common node"""

        set1 = self[n1]
        set2 = self[n2]
        # Convert to set to remove duplicates.
        newset = tuple(set(set1 + set2))

        for n in self[n1]:
            self[n] = newset
        for n in self[n2]:
            self[n] = newset


class Gedge(object):

    def __init__(self, cpt, from_gnode, to_gnode, size, stretch=False):

        self.cpt = cpt
        self.from_gnode = from_gnode
        self.to_gnode = to_gnode
        self.size = size
        self.stretch = stretch
        self.stretched = False

    def __repr__(self):

        return '%s --%s%s--> %s' % (
            self.from_gnode.fmt_name, self.size, '*' if self.stretch else '',
            self.to_gnode.fmt_name)


class Gnode(object):

    def __init__(self, name):

        self.name = name
        self.prev = None
        self.next = None
        self.dist = 0
        self.pos = None
        self.fedges = []
        self.redges = []

    def add_fedge(self, edge):

        for edge1 in self.fedges:
            if (edge1.cpt == edge.cpt and 
                edge1.to_gnode == edge.to_gnode and
                edge1.from_gnode == edge.from_gnode):
                return

        self.fedges.append(edge)

    def add_redge(self, edge):

        for edge1 in self.redges:
            if (edge1.cpt == edge.cpt and 
                edge1.to_gnode == edge.to_gnode and
                edge1.from_gnode == edge.from_gnode):
                return

        self.redges.append(edge)

    @property
    def fmt_name(self):

        if isinstance(self.name, tuple):
            return ', '.join(self.name)
        return self.name

    def __repr__(self):
        
        return self.fmt_name


class Graph(dict):

    def __init__(self, name, nodes):

        self.cnodes = Cnodes(nodes)
        self.name = name

    def __repr__(self):

        s = ''
        for gnode in self.values():
            if gnode.fedges == []:
                s += '@%s %s\n' % (gnode.dist, gnode.fmt_name)
            else:
                for edge in gnode.fedges:
                    s += '@%s %s\n' % (gnode.dist, edge)
        return s

    def __getitem__(self, key):

        # If key is an integer, convert to a string.
        if isinstance(key, int):
            key = '%d' % key

        # Allow indexing by a node name rather than by cnode tuple.
        if key in self.cnodes:
            key = self.cnodes[key]

        return super(Graph, self).__getitem__(key)

    def link(self, n1, n2):
        self.cnodes.link(n1, n2)

    def add(self, cpt, n1, n2, size, stretch):

        if size == 0:
            return

        if size < 0:
            n1, n2 = n2, n1
            size = -size

        if n1 in self.cnodes:
            n1 = self.cnodes[n1]
        if n2 in self.cnodes:
            n2 = self.cnodes[n2]

        gnode1 = self.add_node(n1)
        gnode2 = self.add_node(n2)

        self.add_edges(cpt, gnode1, gnode2, size, stretch)

    def add_node(self, n):
        
        if n not in self:
            self[n] = Gnode(n)
        return self[n]

    def add_edges(self, cpt, gnode1, gnode2, size, stretch):

        gnode1.add_fedge(Gedge(cpt, gnode1, gnode2, size, stretch))
        gnode2.add_redge(Gedge(cpt, gnode2, gnode1, size, stretch))

    @property
    def nodes(self):
        return self.keys()

    @property
    def all_nodes(self):
        # Use set to remove duplicates.
        return list(set(self.cnodes.values()))

    def add_start_nodes(self):

        if 'start' in self:
            return

        gnodes = self.values()
        start = self.add_node('start')
        end = self.add_node('end')

        for gnode in gnodes:
            if gnode.redges == []:
                self.add_edges(None, start, gnode, 0, True)
            if gnode.fedges == []:
                self.add_edges(None, gnode, end, 0, True)

    def assign_fixed(self, gnode, unknown):

        for edge in gnode.fedges:
            if (not edge.stretch and edge.to_gnode.name not in unknown 
                and edge.to_gnode.name != 'end'):
                gnode.pos = edge.to_gnode.dist - edge.size
                return True
        for edge in gnode.redges:
            if (not edge.stretch and edge.to_gnode.name not in unknown
                and edge.to_gnode.name != 'start'):
                gnode.pos = edge.to_gnode.dist + edge.size
                return True
        return False

    def analyse(self, stage=None):

        self.add_start_nodes()

        for gnode in self.values():
            gnode.pos = None

        unknown = set(self.keys())
        unknown.discard('start')
        unknown.discard('end')

        if unknown == set():
            pos = {}
            for n, gnode in self.cnodes.items():
                pos[n] = 0
            return pos, 0

        # Find longest path through the graph.
        # self.longest_forward_path(self['start'])
        self.longest_path_to_known(self['start'])

        try:
            # Nodes on the longest path have known positions.
            gnode = self['end']
            while gnode != None and gnode.name != 'start':
                unknown.discard(gnode.name)
                gnode.path = True
                gnode.pos = gnode.dist
                gnode = gnode.prev.from_gnode
                self['start'].pos = 0
        except AttributeError:
            raise AttributeError("The %s schematic graph is dodgy, probably a "
                                 "component is attached to the wrong nodes:\n%s" % (self.name, self))            

        if stage == 1:
            return
            
        # Assign node positions to nodes with fixed edge lengths to
        # nodes with known positions.  Iterate until no more changes.
        # This stage is not needed but provides a minor optimisation.
        changes = True
        while changes and unknown != set():
            for n in unknown:
                gnode = self[n]
                changes = self.assign_fixed(gnode, unknown)
                if changes:
                    unknown.discard(n)
                    break

        if stage == 2:
            return

        # For each node, find the longest path between nodes with
        # known positions.  The path between the nodes with known
        # positions is then traversed, measuring the total separation
        # and the number of stretchable components.  (Note, there may
        # be multiple paths of the same length; it does not matter
        # which one is chosen.) The stretch is then calculated from
        # the distance of the path, the total separation, and the
        # number of stretches.  If only one known position is found,
        # the node is dangling and the stretch is zero.
        
        # TODO, check for multiple paths with a conflict, for example,
        # a stretchy path of a longser distance than a fixed path.

        # This can be optimised by processing all the nodes on
        # a path between nodes with known positions at the same time.
        for n in unknown:
            gnode = self[n]

            to_stretches = 0
            to_dist = 0
            self.longest_path_to_known(gnode)
            to_gnode = gnode
            while to_gnode.next is not None:
                if to_gnode.next.stretch:
                    to_stretches += 1
                to_dist += to_gnode.next.size
                to_gnode = to_gnode.next.to_gnode

            from_stretches = 0
            from_dist = 0
            self.longest_path_to_known(gnode, False)
            from_gnode = gnode
            while from_gnode.next is not None:
                if from_gnode.next.stretch:
                    from_stretches += 1
                from_dist += from_gnode.next.size
                from_gnode = from_gnode.next.to_gnode

            if from_gnode.name == 'start' and to_gnode.name == 'end':
                # There is a chance that we have naively picked an unlucky gnode.
                # We should defer this gnode until the next pass
                raise ValueError('Disconnected %s schematic graph for node %s\n%s' % (self.name, gnode, self))

            if from_gnode.name == 'start':
                # Have dangling node, so no stretch needed.
                gnode.pos = to_gnode.pos - to_dist
                continue

            if to_gnode.name == 'end':
                # Have dangling node, so no stretch needed.
                gnode.pos = from_gnode.pos + from_dist
                continue

            separation = to_gnode.pos - from_gnode.pos
            extent = to_dist + from_dist
            if extent - separation > 1e-6:
                raise ValueError('Inconsistent %s schematic graph, component will not fit:  separation %s between %s and %s, need %s.\n%s' % (self.name, separation, from_gnode, to_gnode, extent, self))

            if from_stretches == 0:
                gnode.pos = from_gnode.pos + from_dist
            elif to_stretches == 0:
                gnode.pos = to_gnode.pos - to_dist
            else:
                stretch = (separation - extent) / (to_stretches + from_stretches)
                gnode.pos = from_gnode.pos + from_dist + stretch * from_stretches

        self.check_positions()
  
        try:
            pos = {}
            for n, gnode in self.cnodes.items():
                pos[n] = self[gnode].pos

        except KeyError:
            # TODO determine which components are not connected.
            raise KeyError("The %s schematic graph is dodgy, probably a "
                           "component is unattached:\n%s" % (self.name, self))

        distance_max = self['end'].pos

        return pos, distance_max


    def longest_path_to_known(self, start, forward=True):
        """Find longest path through DAG to a node with a known dist."""

        for gnode in self.values():
            gnode.dist = -1
            gnode.prev = None
            gnode.next = None

        def traverse(gnode):

            if gnode.pos is not None:
                return True

            edges = gnode.fedges if forward else gnode.redges
            for edge in edges:
                next_gnode = edge.to_gnode
                dist = gnode.dist + edge.size
                if dist > next_gnode.dist:
                    next_gnode.dist = dist
                    next_gnode.prev = edge
                    gnode.next = edge
                    if traverse(next_gnode):
                        return True
            return False

        start.dist = 0
        try:
            traverse(start)
        except RuntimeError:
            raise RuntimeError(
                ("The %s schematic graph is dodgy, probably a component"
                 " is connected to the wrong node:\n%s") % (self.name, self))

    def check_positions(self):

        for gnode in self.values():
            for edge in gnode.fedges:
                dist = edge.to_gnode.pos - gnode.pos
                if edge.stretch:
                    if dist - edge.size < -1e-6:
                        print('Distance conflict %s for %s between nodes %s and %s,'
                              ' due to incompatible sizes' % (dist,
                                  edge.cpt.name, gnode.fmt_name,
                                  edge.to_gnode.fmt_name))
                else:
                    if abs(dist - edge.size) > 1e-6:
                        print('Stretch conflict %s vs %s for %s between nodes %s and %s,'
                              ' due to incompatible sizes' % (dist, edge.size,
                                  edge.cpt.name, gnode.fmt_name,
                                  edge.to_gnode.fmt_name))

    def dot(self, filename=None, stage=None):

        if filename is None:
            filename = tmpfilename('.png')
            self.dot(filename=filename, stage=stage)
            display_matplotlib(filename)
            return

        base, ext = path.splitext(filename)
        if ext in ('.pdf', '.png'):
            dotfilename = filename + '.dot'
            self.dot(dotfilename, stage=stage)
            system('dot -T %s -o %s %s' % (ext[1:], filename, dotfilename))
            remove(dotfilename)            
            return

        if stage != 0:
            self.analyse(stage=stage)

        dotfile = open (filename, 'w')
        dotfile.write ('strict digraph {\n\tgraph [rankdir=LR];\n')

        for gnode in self.values():
            colour = 'red' if gnode.pos is not None else 'blue'
            if gnode.name in ('start', 'end'):
                colour = 'green'
            dotfile.write ('\t"%s"\t [style=filled, color=%s, xlabel="@%s"];\n' % (gnode.fmt_name, colour, gnode.pos))

        for gnode in self.values():
            for edge in gnode.fedges:
                colour = 'black' if edge.stretch else 'red'
                dotfile.write ('\t"%s" ->\t"%s" [ color=%s, label="%s%s" ];\n' % (
                    gnode.fmt_name, edge.to_gnode.fmt_name, colour, 
                    edge.size, '*' if edge.stretch else ''))

        dotfile.write ('}\n')
        dotfile.close ()


class Node(object):

    def __init__(self, name):

        self.name = name
        self._port = False
        self._count = 0
        parts = name.split('_')
        self.rootname = parts[0] if name[0] != '_' else name
        self.primary = len(parts) == 1
        self.list = []
        self.pos = 'unknown'
        self.pin = False
        
    def __repr__(self):
        return '%s @ (%s)' % (self.name, self.pos)

    def append(self, elt):
        """Add new element to the node"""

        if elt.type == 'P':
            self._port = True

        self.list.append(elt)
        if elt.type not in ('O', ):
            self._count += 1

    @property
    def count(self):
        """Number of elements (including wires but not open-circuits)
        connected to the node"""

        return self._count

    def visible(self, draw_nodes):
        """Return true if node drawn"""

        if draw_nodes == 'all':
            return True

        if '@' in self.name:
            return False

        if self.port:
            return True

        if draw_nodes in ('none', None, False):
            return False
        
        if draw_nodes == 'connections':
            return self.count > 2

        return '_' not in self.name

    @property
    def port(self):
        """Return true if node is a port"""

        return self._port or self.count == 1


class Schematic(object):

    def __init__(self, filename=None, **kwargs):

        self.anon = {}
        self.elements = {}
        self.nodes = {}
        # Shared nodes (with same voltage)
        self.snodes = {}
        self.hints = False

        if filename is not None:
            self.netfile_add(filename)

    def __repr__(self):
        
        return self.netlist()

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

        for attr in ('xgraph', 'ygraph'):
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


    def parse(self, string):
        """The general form is: 'Name Np Nm symbol'
        where Np is the positive nose and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        def tex_name(name, subscript=None):

            if subscript is None:
                subscript = ''

            if len(name) > 1:
                name = r'\mathrm{%s}' % name
            if len(subscript) > 1:
                subscript = r'\mathrm{%s}' % subscript
            if len(subscript) == 0:
                return name
        
            return '%s_{%s}' % (name, subscript)

        if '\n' in string:
            lines = string.split('\n')
            for line in lines:
                self.add(line.strip())
            return

        if string == '':
            return None            

        if string[0] == ';':
            #TODO
            #self.opts.add(string[1:])
            return None

        cpt = parser.parse(string, self)
        if cpt is None:
            return

        # There are two possible labels for a component:
        # 1. Component identifier, e.g., R1
        # 2. Component value, expression, or symbol
        id_label = tex_name(cpt.type, cpt.id)
        value_label = None

        if cpt.type in ('O', 'P', 'W') or id_label.find('#') != -1:
            id_label = None

        if cpt.args != ():

            # TODO, extend for mechanical and acoustical components.
            units_map = {'V': 'V', 'I': 'A', 'R': '$\Omega$',
                         'C': 'F', 'L': 'H'}

            expr = cpt.args[0]
            if cpt.classname in ('Vstep', 'Istep'):
                expr = '(%s) * Heaviside(t)' % expr
                value_label = Expr(expr, cache=False).latex()
            elif cpt.classname in ('Vs', 'Is'):
                value_label = Expr(expr, cache=False).latex()
            elif cpt.classname == 'TF':
                value_label = '1:%s' % expr
            elif cpt.classname not in ('TP',):
                try:
                    value = float(expr)
                    if cpt.type in units_map:
                        value_label = EngFormat(
                            value, units_map[cpt.type]).latex()
                    else:
                        value_label = Expr(expr, cache=False).latex()

                except ValueError:
                    value_label = Expr(expr, cache=False).latex()

        # Currently, we only annnotated the component with the value,
        # expression, or symbol.  If this is not specified, it
        # defaults to the component identifier.  Note, some objects
        # we do not want to label, such as wires and ports.

        cpt.id_label = '' if id_label is None else format_label(id_label)
        cpt.value_label = '' if value_label is None else format_label(value_label)
        cpt.default_label = cpt.id_label if cpt.value_label == '' else cpt.value_label

        return cpt

    def add(self, string):
        """The general form is: 'Name Np Nm symbol'
        where Np is the positive nose and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        cpt = self.parse(string)
        if cpt is None:
            return

        if cpt.opts_string != '':
            self.hints = True

        self._invalidate()

        if cpt.name in self.elements:
            print('Overriding component %s' % cpt.name)
            # Need to search lists and update component.

        self.elements[cpt.name] = cpt

        for node in cpt.vnodes:
            self._node_add(node, cpt)

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

        self.xgraph = Graph('horizontal', self.nodes)
        self.ygraph = Graph('vertical', self.nodes)

        # Use components in orthogonal directions as constraints.  The
        # nodes of orthogonal components get combined into a
        # common node.
        for m, elt in enumerate(self.elements.values()):
            elt.xlink(self.xgraph)
            elt.ylink(self.ygraph)

        # Now form forward and reverse directed graph using components
        # in the desired directions.
        for m, elt in enumerate(self.elements.values()):
            elt.xplace(self.xgraph)
            elt.yplace(self.ygraph)

        xpos, self.width = self.xgraph.analyse()
        ypos, self.height = self.ygraph.analyse()

        scale = self.node_spacing
        for n, node in self.nodes.items():
            node.pos = Pos(xpos[n] * scale, ypos[n] * scale)

    def _make_wires1(self, snode_list):

        num_wires = len(snode_list) - 1
        if num_wires == 0:
            return []

        wires = []

        # TODO: remove overdrawn wires...
        for n in range(num_wires):
            n1 = snode_list[n]
            n2 = snode_list[n + 1]
            
            wires.append(self.parse('W_ %s %s' % (n1, n2)))

        return wires

    def _make_wires(self):
        """Create implicit wires between common nodes."""

        wires = []

        snode_dir = self.snodes

        for m, snode_list in enumerate(snode_dir.values()):
            wires.extend(self._make_wires1(snode_list))

        return wires

    def _tikz_draw(self, style_args='', **kwargs):

        self._positions_calculate()

        # Note, scale does not scale the font size.
        opts = r'scale=%.2f,transform shape,/tikz/circuitikz/bipoles/length=%.2fcm,%s' % (
            self.scale, self.cpt_size, style_args)
        s = r'\begin{tikzpicture}[%s]''\n' % opts

        help = float(kwargs.pop('help_lines', 0))
        if help != 0:
            start = Pos(-0.5, -0.5) * self.node_spacing
            stop = Pos(self.width + 0.5, self.height + 0.5) * self.node_spacing

            s += r'\draw[help lines, blue] (%s) grid [xstep=%s, ystep=%s] (%s);''\n' % (
                start, help, help, stop)

        # Write coordinates
        for n, node in self.nodes.items():
            s += r'  \coordinate (%s) at (%s);''\n' % (
                n, node.pos)

        # Draw components
        for m, elt in enumerate(self.elements.values()):
            s += elt.draw(**kwargs)

        wires = self._make_wires()

        label_nodes = kwargs.get('label_nodes', 'primary')

        # Label primary nodes
        if label_nodes:
            for m, node in enumerate(self.nodes.values()):
                name = node.name
                name = name.split('@')[-1]

                if label_nodes == 'alpha':
                    if not name[0].isalpha():
                        continue
                elif label_nodes == 'primary':
                    if not node.primary:
                        continue
                anchors = {False: 'south east', 
                           'l' : 'west', 'r' : 'east', 
                           't' : 'north', 'b' : 'south'}
                anchor = anchors[node.pin]

                s += r'  \draw {[anchor=%s] (%s) node {%s}};''\n' % (
                    anchor, node.name, name.replace('_', r'\_'))

        s += '  ' + kwargs.pop('append', '')

        s += r'\end{tikzpicture}''\n'

        return s

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

    def tikz_draw(self, filename=None, **kwargs):

        root, ext = path.splitext(filename)

        debug = kwargs.pop('debug', False)
        oversample = float(kwargs.pop('oversample', 2))
        style = kwargs.pop('style', 'american')
        self.cpt_size = float(kwargs.pop('cpt_size', 1.2))
        self.node_spacing = float(kwargs.pop('node_spacing', 2.0))
        self.scale = float(kwargs.pop('scale', 1.0))

        if style == 'american':
            style_args = 'american currents, american voltages'
        elif style == 'british':
            style_args = 'american currents, european voltages'
        elif style == 'european':
            style_args = ('european currents, european voltages'
                          ', european inductors, european resistors')
        else:
            raise ValueError('Unknown style %s' % style)

        # For debugging when do not want to write to file
        nosave = kwargs.pop('nosave', False)

        content = self._tikz_draw(style_args=style_args, **kwargs)
        
        if nosave:
            return

        if debug:
            print('width = %d, height = %d, oversample = %d, cpt_size = %.2f, node_spacing = %.2f, scale = %.2f'
                  % (self.width, self.height, oversample, 
                     self.cpt_size, self.node_spacing, self.scale))
            print(self.nodes)
            # print(self.xgraph.cnodes)
            # print(self.ygraph.cnodes)

        if ext == '.pytex':
            open(filename, 'w').write(content)
            return

        template = ('\\documentclass[a4paper]{standalone}\n'
                    '\\usepackage{circuitikz}\n'
                    '\\usetikzlibrary{fit}\n'
                    '\\begin{document}\n%s\\end{document}')
        content = template % content

        texfilename = filename.replace(ext, '.tex')
        open(texfilename, 'w').write(content)

        if ext == '.tex':
            return

        dirname = path.dirname(texfilename)
        baseroot = path.basename(root)
        cwd = getcwd()
        if dirname != '':
            chdir(path.abspath(dirname))

        system('pdflatex -interaction batchmode %s.tex' % baseroot)

        if dirname != '':
            chdir(cwd)            

        if not debug:
            try:
                remove(root + '.aux')
                remove(root + '.log')
                remove(root + '.tex')
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

    def draw(self, filename=None, opts={}, **kwargs):
        """
        filename specifies the name of the file to produce.  If None,
        the schematic is displayed on the screen.

        kwargs include:
           label_ids: True to show component ids
           label_values: True to display component values
           draw_nodes: True to show nodes
           style: 'american', 'british', or 'european'
           scale: schematic scale factor, default 1.0
           node_spacing: spacing between component nodes, default 2.0
           cpt_size: size of a component, default 1.5
           oversample: oversampling factor for png or pdf files
           help_lines: distance between lines in grid, default 0.0 (disabled)
           debug: True to display debug information
        """

        for key, val in opts.items():
            if key not in kwargs or kwargs[key] is None:
                kwargs[key] = val

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

                pngfilename = tmpfilename('.png')
                self.tikz_draw(pngfilename, **kwargs)

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

                svgfilename = tmpfilename('.svg')
                self.tikz_draw(svgfilename, **kwargs)

                # Create and display SVG image object.
                # Note, there is a problem displaying multiple SVG
                # files since the later ones inherit the namespace of
                # the first ones.
                display_svg(SVG(filename=pngfilename, 
                                width=self.width * 100, height=self.height * 100))
                return

        if filename is None:
            filename = tmpfilename('.png')
            self.tikz_draw(filename=filename, **kwargs)
            display_matplotlib(filename)
            return
        
        self.tikz_draw(filename=filename, **kwargs)

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
