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

Copyright 2014--2019 Michael Hayes, UCECE
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
from .latex import latex_str, format_label
from .expr import Expr
from . import schemcpts
import sympy as sym
from .schemgraph import Graph
from .schemmisc import Pos, Opts
from .netfile import NetfileMixin
from .system import run_dot, run_latex, convert_pdf_png, convert_pdf_svg
from .system import tmpfilename, circuitikz_version, latex_cleanup
from os import path, remove
from collections import OrderedDict
import math

__all__ = ('Schematic', )


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

    def math_latex(self):
        """Make latex math-mode string."""

        return '$' + self.latex() + '$'
        
    def latex(self, trim=True, hundreds=False):
        """If hundreds True format like 100 pF rather than 0.1 nF"""

        prefixes = ('f', 'p', 'n', '$\mu$', 'm', '', 'k', 'M', 'G', 'T')

        sfmax = 3

        value = self.value
        if value == 0:
            return '0' + self.unit

        m = math.log10(abs(value))

        if m < -1 or m >= 3.0:
            if not hundreds:
                m += 1
            n = int(math.floor(m / 3))
            k = int(math.floor(m)) - n * 3
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


class Node(object):

    def __init__(self, name):

        self.name = name
        self._port = False
        self._count = 0
        parts = name.split('_')
        self.rootname = parts[0] if name[0] != '_' else name
        self.primary = len(parts) == 1
        self.elt_list = []
        self.pos = 'unknown'
        self.pinpos = None
        self.pin = False
        # Sanitised name
        self.s = name.replace('.', '@')
        self.label = name
        
    def __repr__(self):
        return '%s @ (%s)' % (self.name, self.pos)

    def append(self, elt):
        """Add new element to the node"""

        if elt.type == 'P':
            self._port = True

        self.elt_list.append(elt)
        if elt.type not in ('O', ):
            self._count += 1

    @property
    def count(self):
        """Number of elements (including wires but not open-circuits)
        connected to the node"""

        return self._count

    def visible(self, draw_nodes):
        """Return true if node drawn.
        `draw_nodes' can be `all', 'none', 'connections', 'primary', None,
        True, or False."""

        if draw_nodes in ('all', True):
            return True

        if self.pin:
            return False

        if self._port:
            return True

        if draw_nodes in ('none', None, False):
            return False

        if '_' in self.name:
            return False
        
        # Implied port
        if self.count == 1:
            return True

        if draw_nodes == 'connections':
            return self.count > 2

        if draw_nodes == 'primary':        
            return True
        
        raise ValueError('Unknown argument %s for draw_nodes' % draw_nodes)

    @property
    def port(self):
        """Return true if node is a port"""

        return self._port or self.count == 1


class Schematic(NetfileMixin):

    def __init__(self, filename=None, **kwargs):

        self.elements = OrderedDict()
        self.nodes = {}
        # Shared nodes (with same voltage)
        self.snodes = {}
        self.hints = False
        self._init_parser(schemcpts)
        self.cpt_size = 1.2
        self.node_spacing = 2.0
        self.scale = 1.0
        self.dummy_node = 0

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
        
        self._netfile_add(filename)

    def add(self, string):
        """Add a component to the netlist.
        The general form is: 'Name Np Nm args'
        where Np is the positive node and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        self._add(string)

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

        if elt.type in ('U', 'MX', 'S', 'SP', 'TP', 'TR'):
            self.nodes[node].pin = True

        if vnode not in self.snodes:
            self.snodes[vnode] = []

        if node not in self.snodes[vnode]:
            self.snodes[vnode].append(node)

    def _cpt_add(self, cpt):

        if cpt.offset != 0 and len(cpt.node_names) == 2:
            # Change component and add wires.

            n1 = cpt.node_names[0]
            n2 = cpt.node_names[1]
            self.dummy_node += 1                            
            on1 = n1 + '_o%d' % self.dummy_node
            self.dummy_node += 1                            
            on2 = n2 + '_o%d' % self.dummy_node

            angle = 90
            size = cpt.offset
            if size < 0:
                size = -size
                angle = -angle
            
            w1 = 'W %s %s; rotate=%s, size=%s' % (n1, on1, cpt.angle + angle, size)
            w2 = 'W %s %s; rotate=%s, size=%s' % (n2, on2, cpt.angle + angle, size)

            self.add(w1)
            self.add(w2)

            # Add open-circuit to ensure alignment.
            o = 'O %s %s; rotate=%s, size=%s' % (n1, n2, cpt.angle, cpt.size)
            self.add(o)
            

            # Rename nodes
            parts = cpt.net.split(' ')
            parts[1] = on1
            parts[2] = on2
            net = ' '.join(parts)
            cpt.opts.strip('offset')
            self.add('%s; %s' % (net, cpt.opts))
            return
        
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

        # There are two possible labels for a component:
        # 1. Component identifier, e.g., R1
        # 2. Component value, expression, or symbol
        id_label = tex_name(cpt.type, cpt.id)
        value_label = None

        if cpt.type in ('O', 'P', 'W') or id_label.find('#') != -1:
            id_label = None

        if cpt.type in ('S', 'U'):
            value_label = ''
            
        if cpt.args != ():

            # TODO, extend for mechanical and acoustical components.
            units_map = {'V': 'V', 'I': 'A', 'R': '$\Omega$',
                         'C': 'F', 'L': 'H'}

            expr = cpt.args[0]
            if cpt.classname in ('Vstep', 'Istep'):
                expr = '(%s) * Heaviside(t)' % expr
                value_label = Expr(expr, cache=False).math_latex()
            elif cpt.classname in ('Vs', 'Is'):
                value_label = Expr(expr, cache=False).math_latex()
            elif cpt.classname == 'TF':
                expr = sym.sympify(expr)
                if expr.is_Pow and expr.args[1] == -1:
                    value_label = '%s:1' % (1 / expr)
                else:
                    value_label = '1:%s' % expr
            elif cpt.classname not in ('TP',):
                try:
                    value = float(expr)
                    if cpt.type in units_map:
                        value_label = EngFormat(
                            value, units_map[cpt.type]).math_latex()
                    else:
                        value_label = Expr(expr, cache=False).math_latex()

                except ValueError:
                    # This catches non numeric arg.
                    value_label = Expr(expr, cache=False).math_latex()

        # Currently, we only annnotated the component with the value,
        # expression, or symbol.  If this is not specified, it
        # defaults to the component identifier.  Note, some objects
        # we do not want to label, such as wires and ports.
        cpt.id_label = '' if id_label is None else format_label(id_label)
        cpt.value_label = cpt.id_label if value_label is None else format_label(value_label)

        if cpt.opts_string != '':
            self.hints = True

        self._invalidate()

        if cpt.name in self.elements:
            print('Overriding component %s' % cpt.name)
            # Need to search lists and update component.

        self.elements[cpt.name] = cpt

        for node in cpt.required_node_names + cpt.extra_node_names:
            self._node_add(node, cpt)

    def check_nodes(self):

        def check_explicit_node(node):

            node_name = node.name            
            for elt in node.elt_list:
                if node_name in elt.explicit_node_names:
                    return True
            return False
        
        for node in self.nodes.values():
            if check_explicit_node(node):
                continue
            
            node_name = node.name
            if '.' not in node_name:
                # Something is wrong...
                raise ValueError('Unknown node %s.' % node_name)

            parts = node_name.split('.')
            cpt_name = '.'.join(parts[0:-1])
            if cpt_name not in self.elements:
                raise ValueError('Unknown component name %s for node %s'
                                 % (cpt_name, node_name))
            node_names = self.elements[cpt_name].node_names
            if node_name not in node_names:
                raise ValueError('Unknown node %s. Possible names: %s' %
                                 (node.name, ', '.join(node_names)))
            
    def make_graphs(self):

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

        self.check_nodes()
        
        self.xgraph = Graph('horizontal', self.nodes)
        self.ygraph = Graph('vertical', self.nodes)

        # Use components in orthogonal directions as constraints.  The
        # nodes of orthogonal components get combined into a
        # common node.
        for m, elt in enumerate(self.elements.values()):

            if elt.offset != 0:
                raise ValueError('offset field should be removed')
            
            elt.xlink(self.xgraph)
            elt.ylink(self.ygraph)

        # Now form forward and reverse directed graph using components
        # in the desired directions.
        # Note, this must be done after the linking step.
        for m, elt in enumerate(self.elements.values()):
            elt.xplace(self.xgraph)
            elt.yplace(self.ygraph)

    def _positions_calculate(self):

        self.make_graphs()

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
            
            wires.append(self._parse('W_ %s %s' % (n1, n2)))

        return wires

    def _make_wires(self):
        """Create implicit wires between common nodes."""

        wires = []

        snode_dir = self.snodes

        for m, snode_list in enumerate(snode_dir.values()):
            wires.extend(self._make_wires1(snode_list))

        return wires

    def _label_nodes(self, **kwargs):

        label_nodes = kwargs.get('label_nodes', 'primary')

        s = ''
        if not label_nodes:
            return s

        for m, node in enumerate(self.nodes.values()):

            name = node.name
            name = name.split('.')[-1]
            
            if label_nodes == 'pins':
                if not node.pin:
                    continue
            elif label_nodes == 'alpha':
                if not node.primary or not name[0].isalpha():
                    continue
            elif label_nodes == 'primary':
                if not node.primary:
                    continue
            anchors = {None: 'south east', 
                       'l' : 'west', 'r' : 'east', 
                       't' : 'north', 'b' : 'south'}
            anchor = anchors[node.pinpos]

            if node.pin and node.pinpos is None:
                continue

            s += r'  \draw[anchor=%s] (%s) node {%s};''\n' % (
                anchor, node.s, node.label.replace('_', r'\_'))
        return s

    def _tikz_draw(self, style_args='', **kwargs):

        self._positions_calculate()

        # Note, scale does not scale the font size.
        opts = r'scale=%.2f,transform shape,/tikz/circuitikz/bipoles/length=%.2fcm,%s' % (
            self.scale, self.cpt_size, style_args)
        s = r'\begin{tikzpicture}[%s]''\n' % opts

        help = float(kwargs.pop('help_lines', 0))
        color = kwargs.pop('color', 'blue')
        if help != 0:
            start = Pos(-0.5, -0.5) * self.node_spacing
            stop = Pos(self.width + 0.5, self.height + 0.5) * self.node_spacing

            s += r'\draw[help lines, %s] (%s) grid [xstep=%s, ystep=%s] (%s);''\n' % (
                color, start, help, help, stop)

        # Write coordinates
        for n in self.nodes.values():
            s += r'  \coordinate (%s) at (%s);''\n' % (n.s, n.pos)

        # Draw components
        for m, elt in enumerate(self.elements.values()):
            s += elt.draw(**kwargs)

        wires = self._make_wires()

        s += self._label_nodes(**kwargs)

        s += '  ' + kwargs.pop('append', '')

        s += r'\end{tikzpicture}''\n'

        return s

    def tikz_draw(self, filename, **kwargs):

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

        self.circuitikz_version = circuitikz_version()
        if self.circuitikz_version is None:
            raise RuntimeError('circuitikz is not installed')

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

        if ext in ('.pytex', '.schtex'):
            open(filename, 'w').write(content)
            return

        # Need amsmath for operatorname
        template = ('\\documentclass[a4paper]{standalone}\n'
                    '\\usepackage{amsmath}\n'
                    '\\usepackage{circuitikz}\n'
                    '\\usetikzlibrary{fit, shapes}\n'
                    '\\begin{document}\n%s\\end{document}')
        content = template % content

        tex_filename = filename.replace(ext, '.tex')
        open(tex_filename, 'w').write(content)

        if ext == '.tex':
            return

        pdf_filename = tex_filename.replace('.tex', '.pdf')
        run_latex(tex_filename)
        if not debug:
            latex_cleanup(tex_filename, pdf_filename)

        if not path.exists(pdf_filename):
            raise RuntimeError('Could not generate %s with pdflatex' % 
                               pdf_filename)

        if ext == '.pdf':
            return

        if ext == '.svg':
            convert_pdf_svg(pdf_filename, root + '.svg')
            if not debug:
                remove(pdf_filename)
            return

        if ext == '.png':
            convert_pdf_png(pdf_filename, root + '.png', oversample)
            if not debug:
                remove(pdf_filename)
            return

        raise RuntimeError('Cannot create file of type %s' % ext)

    def draw(self, filename=None, opts={}, **kwargs):
        """
        filename specifies the name of the file to produce.  If None,
        the schematic is displayed on the screen.

        Note, if using Jupyter, then need to first issue command %matplotlib inline

        kwargs include:
           label_ids: True to show component ids
           label_values: True to display component values
           draw_nodes: True to show all nodes, False to show no nodes, 
             'primary' to show primary nodes,
             'connections' to show nodes that connect more than two components,
             'all' to show all nodes
           label_nodes: True to label all nodes, False to label no nodes, 
             'primary' to label primary nodes,
             'alpha' to label nodes starting with a letter,
             'pins' to label nodes that are pins on a chip,
             'all' to label all nodes
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
                                width=self.width * 100,
                                height=self.height * 100))
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
