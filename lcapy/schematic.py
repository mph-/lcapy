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

Copyright 2014--2020 Michael Hayes, UCECE
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

# Strings starting with ;; are schematic options.  They are parsed in
# netfile.py and added to the opts attribute of the netlist.  They get
# passed to the draw method of the Schematic class.


from __future__ import print_function
from .latex import latex_format_label, latex_format_node_label
from .expr import Expr
from . import schemcpts
import sympy as sym
from .schemgraph import Graph
from .schemmisc import Pos
from .opts import Opts
from .netfile import NetfileMixin
from .system import run_latex, convert_pdf_png, convert_pdf_svg
from .system import tmpfilename, circuitikz_version, latex_cleanup
from os import path, remove
from collections import OrderedDict
import math

__all__ = ('Schematic', )


def display_matplotlib(filename, dpi=150, scale=2):
        
    from matplotlib.pyplot import figure, figaspect
    from matplotlib.image import imread
    
    img = imread(filename)
    h, w, d = img.shape
    width = scale * w / dpi
    height = scale * h / dpi
    fig = figure(figsize=(width, height))
    ax = fig.add_subplot(111)
    ax.imshow(img)
    ax.axis('off')
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)


def png_image_size(filename):

    import struct
    
    with open(filename, 'rb') as f:
        header = f.read(25)

    if (header[:8] != b'\211PNG\r\n\032\n' and (header[12:16] != b'IHDR')):
        raise Exception('%s not a png image' % filename)        
    w, h = struct.unpack('>LL', header[16:24])
    width = int(w)
    height = int(h)
    return width, height

    
class SchematicOpts(Opts):

    def __init__(self):

        super (SchematicOpts, self).__init__(
            {'draw_nodes': 'primary',
             'label_values': True,
             'label_ids': True,
             'label_nodes': 'primary',
             'scale' : 1.0,
             'dpi' : 150,             
             'cpt_size' : 1.5,
             'node_spacing' : 2.0,
             'help_lines' : 0.0,
             'style' : 'american'})


class EngFormat(object):

    def __init__(self, value, unit=''):

        self.value = value
        self.unit = unit

    def latex_math(self):
        """Make latex math-mode string."""

        return '$' + self.latex() + '$'
        
    def latex(self, trim=True, hundreds=False):
        """If hundreds True format like 100 pF rather than 0.1 nF"""

        prefixes = ('f', 'p', 'n', '$\mu$', 'm', '', 'k', 'M', 'G', 'T')

        sfmax = 3

        value = self.value
        if value == 0:
            return '0' + '\,' + r'\mbox{' + self.unit + r'}'

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
        # Primary 1, 2, a, a_3.  Not primary _1, _2, _a, _a_3, 0_3, a_b_c
        self.primary = (name[0] != '_' and len(parts) <= 2) and not (name[0].isdigit() and len(parts) != 1)

        # See MX1.sch for example.
        if '._' in name:
            self.primary = False
        
        self.elt_list = []
        self.pos = 'unknown'
        # Sanitised name
        self.s = name.replace('.', '@')
        self.label = latex_format_node_label(self.name)
        self.labelpos = None        
        self.pin = False
        self.pinlabel = ''
        self.namepos = None                
        self.pinname = ''
        self.pinpos = None
        self.clock = False
        self.auxiliary = None
        # Reference to node that defines the pin
        self.ref = None
        # Node is connected to a ground, etc. (an implicit power
        # supply connection)
        self.implicit = False

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
        """Number of elements (including wires but not open-circuits)
        connected to the node"""

        return self._count

    def belongs(self, cpt_name):
        return self.cptname == cpt_name
    
    def visible(self, draw_nodes):
        """Return true if node drawn.
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
        if self.count == 1:
            return True

        if draw_nodes in ('connections', 'connected'):
            return self.count > 2

        if draw_nodes == 'primary':        
            return self.primary
        
        raise ValueError('Unknown argument %s for draw_nodes' % draw_nodes)

    @property
    def port(self):
        """Return true if node is a port"""

        return self._port or self.count == 1

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
                                                self.count, self.pin))
    
    
class Schematic(NetfileMixin):

    def __init__(self, filename=None, allow_anon=False, **kwargs):

        self.elements = OrderedDict()
        self.nodes = {}
        self.hints = False
        self._init_parser(schemcpts, allow_anon)
        self.cpt_size = 1.2
        self.node_spacing = 2.0
        self.scale = 1.0
        self.dummy_node = 0
        self.context = None
        self.debug = False

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

    def _node_add(self, nodename, elt, auxiliary=False):

        if nodename not in self.nodes:
            self.nodes[nodename] = Node(nodename)
        node = self.nodes[nodename]
            
        node.append(elt)

        if node.auxiliary is None:
            node.auxiliary = auxiliary
        elif node.auxiliary and not auxiliary:
            # An auxiliary node can be made a proper node.
            node.auxiliary = False
        
        return node

    def _cpt_add(self, cpt):

        if cpt.offset != 0 and len(cpt.node_names) == 2:
            # Change component and add wires.

            n1 = cpt.node_names[0]
            n2 = cpt.node_names[1]
            # If have only two parallel components, can used fixed node names.
            self.dummy_node += 1                            
            on1 = n1 + 'off%d' % self.dummy_node
            self.dummy_node += 1                            
            on2 = n2 + 'off%d' % self.dummy_node

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

        if cpt.type in ('S', 'SW', 'U'):
            value_label = ''

        if cpt.args != ():

            # TODO, extend for mechanical and acoustical components.
            units_map = {'V': 'V', 'I': 'A', 'R': '$\Omega$',
                         'C': 'F', 'L': 'H'}

            expr = cpt.args[0]
            if cpt.classname in ('Vstep', 'Istep'):
                expr = '(%s) * Heaviside(t)' % expr
                value_label = Expr(expr, cache=False).latex_math()
            elif cpt.classname in ('Vs', 'Is'):
                value_label = Expr(expr, cache=False).latex_math()
            elif cpt.classname == 'TF':
                expr = sym.sympify(expr)
                if expr.is_Pow and expr.args[1] == -1:
                    value_label = '%s:1' % (1 / expr)
                else:
                    value_label = '1:%s' % expr
            elif cpt.type in ('F', 'H') and len(cpt.args) > 1:
                # This is hard to give a reasonable label since the
                # control current is specified by a voltage source.
                # The user will have to override manually.
                expr = cpt.args[1]
                value_label = Expr(expr, cache=False).latex_math()
            elif cpt.classname not in ('TP',):
                try:
                    # Handle things like 9/1000 that can occur
                    # when substituting cpt values.  
                    value = float(sym.Rational(expr))
                    if cpt.type in units_map:
                        value_label = EngFormat(
                            value, units_map[cpt.type]).latex_math()
                    else:
                        value_label = Expr(expr, cache=False).latex_math()

                except (ValueError, TypeError):
                    # This catches non numeric arg.
                    value_label = Expr(expr, cache=False).latex_math()

        # Currently, we only annnotated the component with the value,
        # expression, or symbol.  If this is not specified, it
        # defaults to the component identifier.  Note, some objects
        # we do not want to label, such as wires and ports.
        cpt.id_label = '' if id_label is None else latex_format_label(id_label)
        cpt.value_label = cpt.id_label if value_label is None else latex_format_label(value_label)

        if cpt.opts_string != '':
            self.hints = True

        self._invalidate()

        if cpt.name in self.elements:
            print('Overriding component %s' % cpt.name)
            # Need to search lists and update component.

        self.elements[cpt.name] = cpt

        for node in cpt.auxiliary_node_names:
            self._node_add(node, cpt, auxiliary=True)

        # Note, an auxiliary node can be trumped...
        for node in cpt.required_node_names:
            self._node_add(node, cpt, auxiliary=False)

    def make_graphs(self, debug=None):

        if debug is None:
            debug = self.debug
        
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

        self.xgraph = Graph('horizontal', self.nodes, debug)
        self.ygraph = Graph('vertical', self.nodes, debug)

        # Use components in orthogonal directions as constraints.  The
        # nodes of orthogonal components get combined into a
        # common node.
        for m, elt in enumerate(self.elements.values()):

            if elt.offset != 0:
                raise ValueError('offset field should be removed')
            if elt.directive or elt.ignore:
                continue
            
            elt.xlink(self.xgraph)
            elt.ylink(self.ygraph)

        # Now form forward and reverse directed graph using components
        # in the desired directions.
        # Note, this must be done after the linking step.
        for m, elt in enumerate(self.elements.values()):
            if elt.directive or elt.ignore:
                continue            
            elt.xplace(self.xgraph)
            elt.yplace(self.ygraph)
            
    def _positions_calculate(self):

        self.make_graphs()

        xpos, self.width = self.xgraph.analyse()
        ypos, self.height = self.ygraph.analyse()

        scale = self.node_spacing
        for n, node in self.nodes.items():
            node.pos = Pos(xpos[n] * scale, ypos[n] * scale)

    def _setup(self):
        # This is called before node positions are assigned.

        for elt in self.elements.values():
            elt.setup()

        for nodename, node in self.nodes.items():
            if node.ref is not None:
                continue
            if node.cptname is not None and not node.implicit:
                raise ValueError('Unreferenced pin connection %s for %s' % (node.name, node.elt_list))

    def _tikz_draw(self, style_args='', **kwargs):

        self.debug = kwargs.pop('debug', False)
        
        self._setup()
        
        self._positions_calculate()

        # Note, scale does not scale the font size.
        opts = ['scale=%.2f' % self.scale,
                'transform shape',
                '/tikz/circuitikz/bipoles/length=%.2fcm' % self.cpt_size]
        opts.append(style_args)
        if 'options' in kwargs:
            opts.append(kwargs.pop('options'))

        if 'font' in kwargs:
            font = kwargs.pop('font')            
            opts.append('font=' + font)
        
        s = r'\begin{tikzpicture}[%s]''\n' % ', '.join(opts)

        # Add preamble        
        if 'preamble' in kwargs:
            s += '  ' + kwargs.pop('preamble') + '\n'
        
        help = float(kwargs.pop('help_lines', 0))
        color = kwargs.pop('color', 'blue')
        if help != 0:
            start = Pos(-0.5, -0.5) * self.node_spacing
            stop = Pos(self.width + 0.5, self.height + 0.5) * self.node_spacing

            s += r'  \draw[help lines, %s] (%s) grid [xstep=%s, ystep=%s] (%s);''\n' % (
                color, start, help, help, stop)

        # Write coordinates.  TODO, not all coordinates are needed
        # so those can be weeded out to simplify the generated file.
        for n in self.nodes.values():
            s += r'  \coordinate (%s) at (%s);''\n' % (n.s, n.pos)

        # keyword args for second pass
        kwargs2 = kwargs.copy()
            
        # Draw components
        for elt in self.elements.values():
            if elt.ignore:
                continue
            if elt.directive:            
                for key, val in elt.opts.items():
                    kwargs[key] = val 

            s += elt.draw(**kwargs)
            s += elt.draw_nodes(**kwargs)
            s += elt.draw_pins()

        # Add the node labels
        for elt in self.elements.values():
            if elt.ignore:
                continue                        
            if elt.directive:
                for key, val in elt.opts.items():
                    kwargs2[key] = val
            s += elt.draw_node_labels(**kwargs2)            

        # Add postamble
        if 'postamble' in kwargs:
            s += '  ' + kwargs.pop('postamble') + '\n'        

        s += r'\end{tikzpicture}''\n'

        return s

    def tikz_draw(self, filename, **kwargs):

        root, ext = path.splitext(filename)

        style = kwargs.pop('style', 'american')
        self.dpi = float(kwargs.pop('dpi', 150))        
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

        self.circuitikz_date, self.circuitikz_version = circuitikz_version()
        if self.circuitikz_date is None:
            raise RuntimeError('circuitikz is not installed')

        content = self._tikz_draw(style_args=style_args, **kwargs)
        
        if nosave:
            return

        if self.debug:
            print('circuitikz version %s (%s)' % (self.circuitikz_version,
                                                  self.circuitikz_date))
            print('width=%d, height=%d, dpi=%d, cpt_size=%.2f, node_spacing=%.2f, scale=%.2f'
                  % (self.width, self.height, self.dpi, 
                     self.cpt_size, self.node_spacing, self.scale))
            print('Nodes:')
            for node in self.nodes.values():
                node.debug()
            # print(self.xgraph.cnodes)
            # print(self.ygraph.cnodes)

        if ext in ('.pytex', '.schtex', '.pgf'):
            open(filename, 'w').write(content)
            return

        # Need amsmath for operatorname
        template = ('\\documentclass[a4paper]{standalone}\n'
                    '\\usepackage{amsmath}\n'
                    '\\usepackage{circuitikz}\n'
                    '\\usetikzlibrary{fit, shapes, arrows}\n'
                    '\\begin{document}\n%s\\end{document}')
        content = template % content

        tex_filename = filename.replace(ext, '.tex')
        open(tex_filename, 'w').write(content)

        if ext == '.tex':
            return

        pdf_filename = tex_filename.replace('.tex', '.pdf')
        run_latex(tex_filename)
        if not self.debug:
            latex_cleanup(tex_filename, pdf_filename)

        if not path.exists(pdf_filename):
            raise RuntimeError('Could not generate %s with pdflatex' % 
                               pdf_filename)

        if ext == '.pdf':
            return

        if ext == '.svg':
            convert_pdf_svg(pdf_filename, root + '.svg')
            if not self.debug:
                remove(pdf_filename)
            return

        if ext == '.png':
            convert_pdf_png(pdf_filename, root + '.png', self.dpi)
            if not self.debug:
                remove(pdf_filename)
            return

        raise RuntimeError('Cannot create file of type %s' % ext)

    def draw(self, filename=None, **kwargs):
        """
        filename specifies the name of the file to produce.  If None,
        the schematic is displayed on the screen.

        Note, if using Jupyter, then need to first issue command %matplotlib inline

        kwargs include:
           label_ids: True to show component ids
           label_values: True to display component values
           draw_nodes: True to show all nodes,
             False or 'none' to show no nodes, 
             'primary' to show primary nodes,
             'connections' to show nodes that connect more than two components,
             'all' to show all nodes
           label_nodes: True to label all nodes,
             False or 'none' to label no nodes, 
             'primary' to label primary nodes,
             'alpha' to label nodes starting with a letter,
             'pins' to label nodes that are pins on a chip,
             'all' to label all nodes,
             'none' to label no nodes
           style: 'american', 'british', or 'european'
           scale: schematic scale factor, default 1.0
           node_spacing: spacing between component nodes, default 2.0
           cpt_size: size of a component, default 1.5
           dpi: dots per inch for png files
           help_lines: distance between lines in grid, default 0.0 (disabled)
           debug: True to display debug information
        """

        # None means don't care, so remove.
        for key in list(kwargs.keys()):
            if kwargs[key] is None:
                kwargs.pop(key)

        # Remove options that may be overridden
        for elt in self.elements.values():
            for key in list(elt.opts.keys()):
                if key in kwargs:
                    elt.opts.remove(key)

        # Default options
        opts = SchematicOpts()
        for key, val in opts.items():
            if key not in kwargs:            
                kwargs[key] = val            

        # Global options (at end of list for historical reasons)
        for eltname in reversed(self.elements):
            elt = self.elements[eltname]
            if not elt.directive:
                break
            for key, val in elt.opts.items():
                # val is a str
                kwargs[key] = val

        def in_ipynb():
            try:
                ip = get_ipython()
                cfg = ip.config

                kernapp = cfg['IPKernelApp']

                # Check if processing ipynb file for Jupyter notebook.
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

        if in_ipynb() and filename is None:
            png = 'png' in kwargs and kwargs.pop('png')
            svg = 'svg' in kwargs and kwargs.pop('svg')
            
            if not png and not svg:
                svg = False

            if svg:
                try:
                    from IPython.display import SVG, display_svg

                    svgfilename = tmpfilename('.svg')
                    self.tikz_draw(svgfilename, **kwargs)
                    
                    # Create and display SVG image object.
                    # Note, there is a problem displaying multiple SVG
                    # files since the later ones inherit the namespace of
                    # the first ones.
                    display_svg(SVG(filename=svgfilename)) 
                    return
                except:
                    pass

            from IPython.display import Image, display_png

            pngfilename = tmpfilename('.png')
            self.tikz_draw(pngfilename, **kwargs)
            
            # Create and display PNG image object.
            # There are two problems:
            # 1. The image metadata (width, height) is ignored
            #    when the ipynb file is loaded.
            # 2. The image metadata (width, height) is not stored
            #    when the ipynb file is written non-interactively.
            width, height = png_image_size(pngfilename)
            # width, height specify the image dimension in pixels
            display_png(Image(filename=pngfilename,
                              width=width, height=height))
            return

        if filename is None:
            filename = tmpfilename('.png')
            # Thicken up lines to reduce aliasing causing them to
            # disappear, especially when using pdftoppm.
            self.tikz_draw(filename=filename, 
                           options='bipoles/thickness=2',
                           **kwargs)
            display_matplotlib(filename, self.dpi)
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
