"""
This module performs schematic drawing using circuitikz from a netlist::

    >>> from lcapy import Schematic
    >>> sch = Schematic('''
    ... P1 1 0.1; down
    ... R1 3 1; right
    ... L1 2 3; right
    ... C1 3 0; down
    ... P2 2 0.2; down
    ... W 0 0.1; right
    ... W 0.2 0.2; right''')
    >>> sch.draw()

Copyright 2014--2022 Michael Hayes, UCECE
"""

# Strings starting with ;; are schematic options.  They are parsed in
# netfile.py and added to the opts attribute of the netlist.  They get
# passed to the draw method of the Schematic class.


from __future__ import print_function
from .config import autoground_default
from .latex import latex_format_label
from .expr import Expr
from . import schemcpts
import sympy as sym
from .engformatter import EngFormatter
from .schemmisc import Pos
from .schemnode import Node
from .schemplacer import schemplacer
from .opts import Opts
from .netfile import NetfileMixin
from .system import LatexRunner, PDFConverter, tmpfilename
from os import path, remove
from collections import OrderedDict
from warnings import warn


__all__ = ('Schematic', )

# Default dots per inch for png files
PNG_DPI = 300

CIRCUITIKZ_MIN_VERSION = '1.4.5'


def display_matplotlib(filename, dpi=PNG_DPI, scale=2):

    from matplotlib.pyplot import figure
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

        super(SchematicOpts, self).__init__(
            {'draw_nodes': 'primary',
             'label_values': True,
             'label_ids': True,
             'annotate_values': False,
             'label_nodes': 'primary',
             'node_label_anchor': 'south east',
             'autoground': 'none',
             'scale': 1.0,
             'dpi': PNG_DPI,
             'cpt_size': 1.5,
             'node_spacing': 2.0,
             'help_lines': 0.0,
             'style': 'american',
             'voltage_dir': 'RP'})


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
        self.debug = 0

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

    def _format_name(self, cpt_type, cpt_id):

        name = cpt_type
        subscript = cpt_id

        if cpt_type == 'REL':
            name = r'\mathcal{R}'
        elif len(name) > 1:
            name = r'\mathrm{%s}' % name

        if subscript != '':
            if len(subscript) > 1:
                subscript = r'\mathrm{%s}' % subscript
            name = name + '_{%s}' % subscript
        return latex_format_label('$' + name + '$')

    def _format_expr(self, expr):

        return Expr(expr, cache=False).latex_math()

    def _format_value_units(self, value, units):

        return EngFormatter().latex_math(value, units)

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

            w1 = 'W %s %s; rotate=%s, size=%s' % (
                n1, on1, cpt.angle + angle, size)
            w2 = 'W %s %s; rotate=%s, size=%s' % (
                n2, on2, cpt.angle + angle, size)

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

        # There are two possible labels for a component:
        # 1. Component name, e.g., R1
        # 2. Component value, expression, or symbol

        id_label = self._format_name(cpt.type, cpt.id)
        value_label = None

        if cpt.type in ('O', 'P', 'W') or id_label.find('#') != -1:
            id_label = None

        if cpt.type in ('S', 'SW', 'U'):
            value_label = ''

        unify = False
        if cpt.args != ():

            # TODO, extend for mechanical and acoustical components.
            units_map = {'V': 'V', 'I': 'A', 'R': '$\Omega$',
                         'C': 'F', 'L': 'H'}

            expr = cpt.args[0]
            if cpt.classname in ('Vstep', 'Istep'):
                expr = '(%s) * Heaviside(t)' % expr
                value_label = self._format_expr(expr)
            elif cpt.classname in ('Vs', 'Is'):
                value_label = self._format_expr(expr)
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
                value_label = self._format_expr(expr)
            elif cpt.classname not in ('TP',):
                try:
                    # Handle things like 9/1000 that can occur
                    # when substituting cpt values.
                    value = float(sym.Rational(expr))
                    if cpt.type in units_map:
                        value_label = self._format_value_units(
                            value, units_map[cpt.type])
                    else:
                        value_label = self._format_expr(expr)

                except (ValueError, TypeError):
                    # This catches non numeric arg.
                    value_label = self._format_expr(expr)

            # Ensure labels are the same when the value is not specified.
            # This will prevent printing the name and value.
            unify = expr == cpt.type + cpt.id

            # Currently, we only annnotate the component with the value,
            # expression, or symbol.  If this is not specified, it
            # defaults to the component identifier.  Note, some objects
            # we do not want to label, such as wires and ports.
        cpt.id_label = '' if id_label is None else latex_format_label(id_label)
        cpt.value_label = cpt.id_label if value_label is None else latex_format_label(
            value_label)

        if unify:
            cpt.value_label = cpt.id_label

        if cpt.opts_string != '':
            self.hints = True

        if cpt.name in self.elements:
            warn('Overriding component %s' % cpt.name)
            # Need to search lists and update component.

        self.elements[cpt.name] = cpt

        # Perhaps don't show nodes if cpt invisible?
        if not cpt.ignore:

            for node in cpt.auxiliary_node_names:
                self._node_add(node, cpt, auxiliary=True)

            # Note, an auxiliary node can be trumped...
            for node in cpt.required_node_names:
                self._node_add(node, cpt, auxiliary=False)

    def _setup(self, autoground):
        # This is called before node positions are assigned.

        if autoground in ('true', 'True', True):
            autoground = autoground_default
        elif autoground in (False, None, 'false', 'none'):
            autoground = False

        if autoground:
            for elt in self.elements.values():
                elt.autoground(autoground)

        for elt in self.elements.values():
            elt.split_nodes()

        for elt in self.elements.values():
            elt.setup()

        for nodename, node in self.nodes.items():
            if node.ref is not None:
                continue
            if node.cptname is not None and not node.implicit:
                raise ValueError('Unreferenced pin connection %s for %s' % (
                    node.name, node.elt_list))

    def make_graphs(self, debug=0):

        placer = schemplacer(self.elements, self.nodes, 'graph', debug)
        placer._make_graphs()
        return placer.xgraph, placer.ygraph

    def _positions_calculate(self, method='graph', debug=False):

        if self.debug & 4:
            print('Creating graphs')
        placer = schemplacer(self.elements, self.nodes, method, debug)
        if self.debug & 4:
            print('Solving graphs')
        self.width, self.height = placer.solve(self.node_spacing)

    def _tikz_draw(self, style_args='', **kwargs):

        method = kwargs.pop('method', 'graph')
        autoground = kwargs.get('autoground', False)

        self._setup(autoground)

        self._positions_calculate(method, self.debug)

        # Note, scale does not scale the font size.
        opts = ['scale=%.2f' % self.scale,
                'transform shape',
                '/tikz/circuitikz/bipoles/length=%.2fcm' % self.cpt_size]
        opts.append(style_args)
        if 'options' in kwargs:
            opts.append(kwargs.pop('options'))

        global_options = ['font', 'voltage_dir', 'color']
        for opt in global_options:
            if opt in kwargs:
                opts.append(opt.replace('_', ' ') + '=' + kwargs.pop(opt))

        s = r'\begin{tikzpicture}[%s]''\n' % ', '.join(opts)

        # Add preamble
        if 'preamble' in kwargs:
            s += '  ' + kwargs.pop('preamble') + '\n'

        help_lines = float(kwargs.pop('help_lines', 0))
        help_lines_color = kwargs.pop('help_lines_color', 'blue')
        if help_lines != 0:
            start = Pos(-0.5, -0.5) * self.node_spacing
            stop = Pos(self.width + 0.5, self.height + 0.5) * self.node_spacing

            s += r'  \draw[help lines, %s] (%s) grid [xstep=%s, ystep=%s] (%s);''\n' % (
                help_lines_color, start, help_lines, help_lines, stop)

        # Write coordinates.  TODO, not all coordinates are needed
        # so those can be weeded out to simplify the generated file.
        for n in self.nodes.values():
            s += r'  \coordinate (%s) at (%s);''\n' % (n.s, n.pos)

        # Keyword args for second pass
        kwargs2 = kwargs.copy()

        # Pass 1: Draw components
        for elt in self.elements.values():
            if elt.ignore:
                continue
            if elt.directive:
                for key, val in elt.opts.items():
                    # Directive overrides kwargs
                    kwargs[key] = val

            s += elt.draw(**kwargs)
            s += elt.draw_nodes(**kwargs)
            s += elt.draw_pins()

        # Pass 2: Add the node labels
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

        self.debug = kwargs.pop('debug', 0)
        style = kwargs.pop('style', 'american')
        self.dpi = float(kwargs.pop('dpi', PNG_DPI))
        self.cpt_size = float(kwargs.pop('cpt_size', 1.2))
        self.node_spacing = float(kwargs.pop('node_spacing', 2.0))
        self.scale = float(kwargs.pop('scale', 1.0))
        png_converter = kwargs.pop('png_converter', None)

        include = ''
        if 'include' in kwargs:
            include += kwargs.pop('include')
        if 'includefile' in kwargs:
            file = open(kwargs.pop('includefile'))
            include += file.read()
            file.close()

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

        latexrunner = LatexRunner((self.debug & 2) != 0)

        content = self._tikz_draw(style_args=style_args, **kwargs)

        if nosave:
            return

        if self.debug & 1:
            date, version = latexrunner.find_circuitikz_version()
            if date is None:
                raise RuntimeError('circuitikz is not installed')

            print('circuitikz version %s (%s)' % (version, date))
            print('width=%d, height=%d, dpi=%d, cpt_size=%.2f, node_spacing=%.2f, scale=%.2f'
                  % (self.width, self.height, self.dpi,
                     self.cpt_size, self.node_spacing, self.scale))
            print('Nodes:')
            for node in self.nodes.values():
                node.debug()

        if ext in ('.pytex', '.schtex', '.pgf'):
            if include != '':
                warn('Include option ignored for non-standalone file')
            open(filename, 'w').write(content)
            return

        # Need amsmath for operatorname
        template = ('\\documentclass[a4paper]{standalone}\n'
                    '\\usepackage{amsmath}\n'
                    '\\usepackage{circuitikz}\n'
                    '\\usetikzlibrary{fit, shapes, arrows, patterns, decorations.text, decorations.markings}\n'
                    '%s\n'
                    '\\begin{document}\n%s\\end{document}')
        content = template % (include, content)
        tex_filename = filename.replace(ext, '.tex')
        open(tex_filename, 'w').write(content)

        if ext == '.tex':
            return

        pdf_filename = tex_filename.replace('.tex', '.pdf')
        latexrunner.run(tex_filename)

        date, version = latexrunner.extract_circuitikz_version(tex_filename)
        if date is None:
            raise RuntimeError('circuitikz is not installed')

        if version < CIRCUITIKZ_MIN_VERSION:
            warn('Have circuitikz version %s; should upgrade to %s or later'
                 % (version, CIRCUITIKZ_MIN_VERSION))

        if not (self.debug & 1):
            latexrunner.cleanup(tex_filename, pdf_filename)

        if not path.exists(pdf_filename):
            raise RuntimeError('Could not generate %s with pdflatex' %
                               pdf_filename)

        if ext == '.pdf':
            return

        if ext == '.svg':
            pdfconverter = PDFConverter((self.debug & 2) != 0)
            pdfconverter.to_svg(pdf_filename, root + '.svg')
            if not (self.debug & 1):
                remove(pdf_filename)
            return

        if ext == '.png':
            pdfconverter = PDFConverter((self.debug & 2) != 0)
            pdfconverter.to_png(pdf_filename, root + '.png', dpi=self.dpi,
                                method=png_converter)
            if not (self.debug & 1):
                remove(pdf_filename)
            return

        raise RuntimeError('Cannot create file of type %s' % ext)

    def draw(self, filename=None, **kwargs):
        """
        filename specifies the name of the file to produce.  If None,
        the schematic is displayed on the screen.

        Note, if using Jupyter, then need to first issue command %matplotlib inline

        kwargs include:
           'label_ids': True to show component ids
           'label_values': True to display component values
           'annotate_values': True to display component values as separate label
           'draw_nodes': True to show all nodes,
             False or 'none' to show no nodes,
             'primary' to show primary nodes,
             'connections' to show nodes that connect more than two components,
             'all' to show all nodes
           'label_nodes': True to label all nodes,
             False or 'none' to label no nodes,
             'primary' to label primary nodes (nodes without an underscore),
             'alpha' to label nodes starting with a letter,
             'pins' to label nodes that are pins on a chip,
             'all' to label all nodes,
             'none' to label no nodes
           'node_label_anchor': where to position node label (default south east)
           'include': name of file to include before \\begin{document}
           'style': 'american', 'british', or 'european'
           'scale': schematic scale factor, default 1.0
           'node_spacing': spacing between component nodes, default 2.0
           'cpt_size': size of a component, default 1.5
           'dpi': dots per inch for png files, default 300
           'help_lines': distance between lines in grid, default 0 (disabled)
           'debug': non-zero to display debug information
        """

        # None means don't care, so remove
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

        # Global options (usually at end of netlist for historical reasons)
        # These can be overwritten
        for eltname in self.elements:
            elt = self.elements[eltname]
            if not elt.directive:
                continue
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
            warn('No schematic drawing hints provided!')

        png = kwargs.pop('png', False)
        svg = kwargs.pop('svg', False)

        if in_ipynb() and filename is None:

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

            png_filename = tmpfilename('.png')
            self.tikz_draw(png_filename, **kwargs)

            # Create and display PNG image object.
            # There are two problems:
            # 1. The image metadata (width, height) is ignored
            #    when the ipynb file is loaded.
            # 2. The image metadata (width, height) is not stored
            #    when the ipynb file is written non-interactively.
            width, height = png_image_size(png_filename)
            # width, height specify the image dimension in pixels
            display_png(Image(filename=png_filename,
                              width=width, height=height))
            if not (self.debug & 1):
                remove(png_filename)
            return

        if filename is None:
            filename = tmpfilename('.png')
            # Thicken up lines to reduce aliasing causing them to
            # disappear, especially when using pdftoppm.
            self.tikz_draw(filename=filename,
                           options='/tikz/circuitikz/bipoles/thickness=2',
                           **kwargs)
            display_matplotlib(filename, self.dpi)
            if not (self.debug & 1):
                remove(filename)
            return

        self.tikz_draw(filename=filename, **kwargs)

    def pdb(self):
        """Enter the python debugger."""

        import pdb
        pdb.set_trace()
        return self


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
