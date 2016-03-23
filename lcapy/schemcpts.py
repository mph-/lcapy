"""
This module defines and draws the schematic components using
circuitikz.   The components are defined at the bottom of this file.

Copyright 2015, 2016 Michael Hayes, UCECE
"""


from __future__ import print_function
from lcapy.latex import latex_str, format_label
from lcapy.schemmisc import Pos, Opts
import numpy as np
import sys

module = sys.modules[__name__]


class Cpt(object):

    voltage_keys = ('v', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<',
                    'v<', 'v>')
    current_keys = ('i', 'i_', 'i^', 'i_>',  'i_<', 'i^>', 'i^<',
                    'i>_', 'i<_', 'i>^', 'i<^', 'i>', 'i<')
    label_keys = ('l', 'l_', 'l^')
    implicit_keys =  ('implicit', 'ground', 'sground', 'rground')
    # The following keys do not get passed through to circuitikz.
    misc_keys = ('left', 'right', 'up', 'down', 'rotate', 'size',
                 'mirror', 'scale', 'invisible', 'variable', 'fixed')

    can_rotate = True
    can_scale = False
    can_mirror = False
    can_stretch = True

    def anon(self, cpt_type):

        sch = self.sch
        if cpt_type not in sch.anon:
            sch.anon[cpt_type] = 0
        sch.anon[cpt_type] += 1        
        return str(sch.anon[cpt_type])

    def __init__(self, sch, cpt_type, cpt_id, string, opts_string, nodes, *args):

        self.sch = sch
        self.type = cpt_type
        self.id = cpt_id

        if cpt_id == '' and sch is not None:
            cpt_id = 'anon' + self.anon(cpt_type)

        name = self.type + cpt_id

        self.net = string.split(';')[0]
        self.opts_string = opts_string
        # There are four sets of nodes:
        # 1. nodes are the names of the electrical nodes for a cpt.
        # 2. vnodes are the subset of the electrical nodes for a cpt that are drawn.
        # 3. dvnodes includes vnodes plus inherited nodes from other cpts (K).  The
        # lookup of this attribute is deferred until cpts are drawn.
        # 4. dnodes are the subset of dvnodes that are shown.
        self.nodes = nodes
        self.name = name
        self.args = args
        self.classname = self.__class__.__name__

        # Drawing hints
        self.opts = Opts(opts_string)

    def __repr__(self):
        return self.__str__()

    def __str__(self):

        if self.opts == {}:
            return self.net
        return self.net + '; ' + str(self.opts)

    @property
    def size(self):
        """Component size between its nodes"""
        if 'size' in self.opts:
            val = self.opts['size']
        elif self.right:
            val = self.opts['right']
        elif self.down:
            val = self.opts['down']
        elif self.left:
            val = self.opts['left']
        elif self.up:
            val = self.opts['up']
        else:
            val = 1
        if val == '':
            val = 1
        return float(val)

    @property
    def scale(self):
        """This is only used for scaling an opamp"""
        return float(self.opts.get('scale', 1.0))

    @property
    def down(self):
        return 'down' in self.opts

    @property
    def up(self):
        return 'up' in self.opts

    @property
    def left(self):
        return 'left' in self.opts

    @property
    def right(self):
        return 'right' in self.opts

    @property
    def horizontal(self):
        return self.left or self.right

    @property
    def vertical(self):
        return self.up or self.down

    def boolattr(self, opt):

        if opt not in self.opts:
            return False
        if self.opts[opt] == '':
            return True
        return self.opts[opt]

    @property
    def mirror(self):
        return self.boolattr('mirror')

    @property
    def invisible(self):
        return self.boolattr('invisible')

    @property
    def variable(self):
        return self.boolattr('variable')

    @property
    def fixed(self):
        return self.boolattr('fixed')

    @property
    def stretch(self):
        return self.can_stretch and not self.fixed

    @property
    def angle(self):
        """Return rotation angle"""
        if self.right:
            angle = 0   
        elif self.down:
            angle = -90
        elif self.left:
            angle = 180
        elif self.up:
            angle = 90
        else:
            angle = -90 if self.type in ('O', 'P') else 0
        
        if 'rotate' in self.opts:
            angle += float(self.opts['rotate'])
        return angle

    @property
    def R(self):
        """Return rotation matrix"""
        angle = self.angle
        
        Rdict = {0: ((1, 0), (0, 1)),
                 90: ((0, 1), (-1, 0)),
                 180: ((-1, 0), (0, -1)),
                 -180: ((-1, 0), (0, -1)),
                 -90: ((0, -1), (1, 0))}
        if angle in Rdict:
            return np.array(Rdict[angle])
        
        t = angle / 180.0 * np.pi
        return np.array(((np.cos(t), np.sin(t)),
                         (-np.sin(t), np.cos(t))))

    @property
    def vnodes(self):
        '''Visible node names'''
        return self.nodes

    @property
    def dvnodes(self):
        '''Nodes used to construct schematic (these are deferred vnodes)'''
        return [self.sch.nodes[n] for n in self.vnodes]

    @property
    def dnodes(self):
        '''Nodes that are drawn'''
        return self.dvnodes

    @property
    def coords(self):
        """Return coordinates of each of the nodes"""
        raise NotImplementedError('coords method not implemented for %s' % self)

    @property
    def tcoords(self):
        """Transformed coordinates for each of the nodes"""
        if hasattr(self, '_tcoords'):
            return self._tcoords

        self._tcoords = np.dot(np.array(self.coords), self.R)
        return self._tcoords

    @property
    def xvals(self):
        return self.tcoords[:, 0]

    @property
    def yvals(self):
        return self.tcoords[:, 1]

    def xlink(self, graphs):

        xvals = self.xvals
        for m1, n1 in enumerate(self.dvnodes):
            for m2, n2 in enumerate(self.dvnodes[m1 + 1:], m1 + 1):
                if xvals[m2] == xvals[m1]:
                    graphs.link(n1.name, n2.name)

    def ylink(self, graphs):

        yvals = self.yvals
        for m1, n1 in enumerate(self.dvnodes):
            for m2, n2 in enumerate(self.dvnodes[m1 + 1:], m1 + 1):
                if yvals[m2] == yvals[m1]:
                    graphs.link(n1.name, n2.name)

    def place(self, graphs, vals):
        
        size = self.size
        idx = np.argsort(vals)[::-1]
        for i in range(len(idx) - 1):
            m1 = idx[i]
            m2 = idx[i + 1]
            n1 = self.dvnodes[m1]
            n2 = self.dvnodes[m2]
            value = (vals[m2] - vals[m1]) * size
            graphs.add(self, n1.name, n2.name, value, self.stretch)

    def xplace(self, graphs):
        self.place(graphs, self.xvals)

    def yplace(self, graphs):
        self.place(graphs, self.yvals)

    def midpoint(self, node1, node2):
        return (node1.pos + node2.pos) * 0.5

    def _node_str(self, node1, **kwargs):

        draw_nodes = kwargs.get('draw_nodes', True)

        s = ''
        if node1.visible(draw_nodes) and not node1.pin:
            s = 'o' if node1.port else '*'
        return s

    def _node_pair_str(self, node1, node2, **kwargs):

        # Create o-o o-* *-* etc.
        s = self._node_str(node1, **kwargs)
        s += '-'
        s += self._node_str(node2, **kwargs)

        if s == '-':
            s = ''
        
        return s

    def _draw_node(self, n, **kwargs):

        draw_nodes = kwargs.get('draw_nodes', True)        

        s = ''
        if not draw_nodes:
            return s

        if not n.visible(draw_nodes) or n.pin:
            return s

        if n.port:
            s = r'  \draw (%s) node[ocirc] {};''\n' % n.name
        else:
            s = r'  \draw (%s) node[circ] {};''\n' % n.name

        return s

    def _draw_nodes(self, **kwargs):

        s = ''
        for n in self.dnodes:
            s += self._draw_node(n, **kwargs)
        return s

    def draw(self, **kwargs):
        raise NotImplementedError('draw method not implemented for %s' % self)

    def opts_str(self, choices):

        def fmt(key, val):
            return '%s=%s' % (key, format_label(val))

        return ','.join([fmt(key, val) 
                         for key, val in self.opts.items()
                         if key in choices])

    @property
    def voltage_str(self):

        return self.opts_str(self.voltage_keys)

    @property
    def current_str(self):

        return self.opts_str(self.current_keys)

    @property
    def label_str(self):

        return self.opts_str(self.label_keys)

    @property
    def args_str(self):

        def fmt(key, val):
            return '%s=%s' % (key, format_label(val))

        return ','.join([fmt(key, val) 
                         for key, val in self.opts.items()
                         if key not in self.voltage_keys + self.current_keys + self.label_keys + self.misc_keys + self.implicit_keys])

    def label(self, **kwargs):

        label_values = kwargs.get('label_values', True)
        label_str = self.default_label if label_values else ''

        # Override label if specified.  There are no placement options.
        str =  ','.join([format_label(val)
                         for key, val in self.opts.items()
                         if key in ('l', )])

        if str != '':
            label_str = str
        return label_str

    def check(self):
        """Check schematic options and return True if component is to be drawn"""

        if not self.can_rotate and self.angle != 0:
            raise ValueError('Cannot rotate component %s' % self.name)

        if not self.can_scale and self.scale != 1:
            raise ValueError('Cannot scale component %s' % self.name)

        if not self.can_mirror and self.mirror:
            raise ValueError('Cannot mirror component %s' % self.name)

        if self.left + self.right + self.up + self.down > 1:
            raise ValueError('Mutually exclusive drawing directions for %s' % self.name)
        
        return not self.invisible

    def tf(self, centre, offset):
        """Transform coordinate"""

        if isinstance(offset[0], tuple):
            return [self.tf(centre, offset1) for offset1 in offset]

        return centre + np.dot(offset, self.R) * self.scale


    def xtf(self, centre, offset):
        """Transform coordinate but with x-scaling only"""
        
        if isinstance(offset[0], tuple):
            return [self.xtf(centre, offset1) for offset1 in offset]

        offset = (offset[0] * self.scale, offset[1])
        return centre + np.dot(offset, self.R)


class Transistor(Cpt):
    """Transistor"""
    
    npos = ((1, 1.5), (0, 0.75), (1, 0))
    ppos = ((1, 0), (0, 0.75), (1, 1.5))

    can_mirror = True
    can_stretch = False
    can_scale = True

    @property
    def coords(self):
        if self.classname in ('Qpnp', 'Mpmos', 'Jpjf'):
            return self.npos if self.mirror else self.ppos
        else:
            return self.ppos if self.mirror else self.npos

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3 = self.dvnodes
        centre = (n1.pos + n3.pos) * 0.5

        s = r'  \draw (%s) node[%s, %s, scale=%s, rotate=%d] (%s) {};''\n' % (
            centre, self.tikz_cpt, self.args_str, 2 * self.scale,
            self.angle, self.name)
        s += r'  \draw (%s) node[] {%s};''\n'% (centre, self.label(**kwargs))

        # Add additional wires.  These help to compensate for the
        # slight differences in sizes of the different transistors.
        if self.tikz_cpt in ('pnp', 'npn'):
            s += r'  \draw (%s.C) -- (%s) (%s.B) -- (%s) (%s.E) -- (%s);''\n' % (
                self.name, n1.s, self.name, n2.s, self.name, n3.s)
        else:
            s += r'  \draw (%s.D) -- (%s) (%s.G) -- (%s) (%s.S) -- (%s);''\n' % (
                self.name, n1.s, self.name, n2.s, self.name, n3.s)

        s += self._draw_nodes(**kwargs)
        return s


class JFET(Transistor):
    """Transistor"""
    
    npos = ((1, 1.5), (0, 0.48), (1, 0))
    ppos = ((1, 0), (0, 1.02), (1, 1.5))


class MOSFET(Transistor):
    """Transistor"""
    
    npos = ((0.85, 1.52), (-0.25, 0.76), (0.85, 0))
    ppos = ((0.85, 0), (-0.25, 0.76), (0.85, 1.52))


class TwoPort(Cpt):
    """Two-port"""

    can_rotate = False

    @property
    def coords(self):
        return ((1, 1), (1, 0), (0, 1), (0, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        # TODO, fix positions if component rotated.

        n1, n2, n3, n4 = self.dvnodes
        width = n2.pos.x - n4.pos.x
        centre = (n1.pos + n2.pos + n3.pos + n4.pos) * 0.25
        top = Pos(centre.x, n1.pos.y + 0.15)

        titlestr = "%s-parameter two-port" % self.args[2]

        s = r'  \draw[thick] (%s) -- (%s) -- (%s) -- (%s) -- (%s);''\n' % (
            n4, n3, n1, n2, n4)
        s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
            centre, width, titlestr, self.name)
        s += r'  \draw (%s) node[minimum width=%.1f] {%s};''\n' % (
            top, width, self.label(**kwargs))

        s += self._draw_nodes(**kwargs)
        return s


class TL(Cpt):
    """Transmission line"""

    can_scale = True
    can_rotate = False

    @property
    def coords(self):
        return ((1.25, 0.5), (1.25, 0), (0, 0.5), (0, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.dvnodes

        centre = (n1.pos + n3.pos) * 0.5
        q = self.xtf(centre, ((-1.25, 0), (0.65, 0),
                              (0.525, -0.29), (-0.7, -0.29)))

        xs = self.scale
        # Rotation creates an ellipse!
        s = r'  \draw (%s) node[tlinestub,xscale=%s] {};''\n' % (
            centre + Pos(-1.3 * xs, 0), self.scale)
        s += r'  \draw (%s) node[] {%s};''\n'% (centre, self.label(**kwargs))
        s += r'  \draw (%s) -- (%s);''\n' % (q[0], n3.s)
        s += r'  \draw (%s) -- (%s);''\n' % (q[1], n1.s)
        s += r'  \draw (%s) |- (%s);''\n' % (q[2], n2.s)
        s += r'  \draw (%s) |- (%s);''\n' % (q[3], n4.s)
        s += self._draw_nodes(**kwargs)
        return s


class TF1(TwoPort):
    """Transformer"""

    can_rotate = True

    @property
    def coords(self):
        return ((0.5, 1), (0.5, 0), (0, 1), (0, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        link = kwargs.get('link', True)

        p = [node.pos for node in self.dvnodes]

        centre = (p[0] + p[1] + p[2] + p[3]) * 0.25
        q = self.tf(centre, ((-0.6, 0.6), (0.6, 0.6), (0, 0.65)))
        primary_dot = q[0]
        secondary_dot = q[1]
        labelpos = q[2]

        s = r'  \draw (%s) node[circ] {};''\n' % primary_dot
        s += r'  \draw (%s) node[circ] {};''\n' % secondary_dot
        s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
            labelpos, 0.5, self.name, self.label(**kwargs))

        if link:
            # TODO: allow for rotation
            width = p[0].x - p[2].x
            arcpos = Pos((p[0].x + p[2].x) / 2, secondary_dot.y - width / 2 + 0.2)

            s += r'  \draw [<->] ([shift=(45:%.2f)]%s) arc(45:135:%.2f);''\n' % (
                width / 2, arcpos, width / 2)

        if self.classname in ('TFcore', 'TFtapcore'):
            # Draw core
            q = self.tf(centre, ((-0.1, -0.4), (-0.1, 0.4),
                                 (0.1, -0.4), (0.1, 0.4)))
            s += r'  \draw[thick] (%s) -- (%s);''\n' % (q[0], q[1])
            s += r'  \draw[thick] (%s) -- (%s);''\n' % (q[2], q[3])

        return s


class Transformer(TF1):
    """Transformer"""

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.dvnodes

        s = r'  \draw (%s) to [inductor] (%s);''\n' % (n3.s, n4.s)
        s += r'  \draw (%s) to [inductor] (%s);''\n' % (n2.s, n1.s)

        s += super(Transformer, self).draw(link=False, **kwargs)
        s += self._draw_nodes(**kwargs)
        return s


class TFtap(TF1):
    """Transformer"""
    can_rotate = True


    @property
    def coords(self):
        return ((0.5, 1), (0.5, 0), (0, 1), (0, 0), (-0.125, 0.55), (0.625, 0.55))

    @property
    def dnodes(self):
        # Do not draw the taps.
        return self.dvnodes[0:4]

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4, n5, n6 = self.dvnodes

        s = r'  \draw (%s) to [inductor] (%s);''\n' % (n3.s, n4.s)
        s += r'  \draw (%s) to [inductor] (%s);''\n' % (n2.s, n1.s)

        s += super(TFtap, self).draw(link=False, **kwargs)
        s += self._draw_nodes(**kwargs)
        return s


class K(TF1):
    """Mutual coupling"""

    def __init__(self, sch, cpt_type, cpt_id, string, opts_string, nodes, *args):

        self.Lname1 = args[0]
        self.Lname2 = args[1]
        super (K, self).__init__(sch, cpt_type, cpt_id, string, opts_string, nodes, *args[2:])

    @property
    def dvnodes(self):

        # L1 and L2 need to be previously defined so we can find their nodes.
        L1 = self.sch.elements[self.Lname1]
        L2 = self.sch.elements[self.Lname2]
        return [self.sch.nodes[n] for n in L1.nodes + L2.nodes]


class OnePort(Cpt):
    """OnePort"""

    can_mirror = True
    can_scale = True

    @property
    def coords(self):
        return ((0, 0), (1, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        label_values = kwargs.get('label_values', True)
        label_ids = kwargs.get('label_ids', True)

        n1, n2 = self.dvnodes

        tikz_cpt = self.tikz_cpt
        if self.variable:
            if self.type in ('C', 'R', 'L'):
                tikz_cpt = 'v' + tikz_cpt
            else:
                raise Error('Component %s not variable' % self.name)

        label_pos = '_'
        voltage_pos = '^'
        if self.type in ('V', 'I', 'E', 'F', 'G', 'H', 'BAT'):

            # circuitikz expects the positive node first, except for
            # voltage and current sources!  So swap the nodes
            # otherwise they are drawn the wrong way around.
            n1, n2 = n2, n1

            if self.right or self.up:
                # Draw label on LHS for vertical cpt and below
                # for horizontal cpt.
                label_pos = '^'
                voltage_pos = '_'
        else:
            if self.left or self.down:
                # Draw label on LHS for vertical cpt and below
                # for horizontal cpt.
                label_pos = '^'
                voltage_pos = '_'

        # Add modifier to place voltage label on other side
        # from component identifier label.
        if 'v' in self.opts:
            self.opts['v' + voltage_pos] = self.opts.pop('v')

        # Reversed voltage.
        if 'vr' in self.opts:
            self.opts['v' + voltage_pos + '>'] = self.opts.pop('vr')

        current_pos = label_pos
        # Add modifier to place current label on other side
        # from voltage marks.
        if 'i' in self.opts:
            self.opts['i' + current_pos] = self.opts.pop('i')

        # Reversed current.
        if 'ir' in self.opts:
            self.opts['i' + current_pos + '<'] = self.opts.pop('ir')

        if 'l' in self.opts:
            self.opts['l' + label_pos] = self.opts.pop('l')

        node_pair_str = self._node_pair_str(n1, n2, **kwargs)

        args_str1 = ','.join([self.args_str])
        args_str2 = ','.join([self.voltage_str, self.current_str])

        if self.mirror:
            args_str += ',mirror'

        if self.scale != 1.0:
            args_str2 += ',bipoles/length=%scm' % (self.sch.cpt_size * self.scale)

        # Generate default label.
        if (label_ids and label_values and self.id_label != '' 
            and self.value_label):
            label_str = r'l%s={%s}{=%s}' % (label_pos, self.id_label,
                                            self.value_label)
        elif label_ids and self.id_label != '':
            label_str = r'l%s=%s' % (label_pos, self.id_label)
        elif label_values and self.value_label != '':
            label_str = r'l%s=%s' % (label_pos, self.value_label)
        else:
            label_str = ''

        # Override label if specified.
        if self.label_str != '':
            label_str = self.label_str
            
        s = r'  \draw[%s] (%s) to [%s,%s,%s,%s,n=%s] (%s);''\n' % (
            args_str1, n1.s, tikz_cpt, label_str, args_str2,
            node_pair_str, self.name, n2.s)
        return s


class VCS(OnePort):
    """Voltage controlled source"""

    can_rotate = True

    @property
    def vnodes(self):
        return self.nodes[0:2]


class CCS(OnePort):
    """Current controlled source"""

    can_rotate = True


class Opamp(Cpt):

    # The Nm node is not used (ground).
    ppos = ((2.4, 0.0), (0, 0.5), (0, -0.5))
    npos = ((2.4, 0.0), (0, -0.5), (0, 0.5))

    can_scale = True
    can_mirror = True
    can_stretch = False

    @property
    def vnodes(self):
        return (self.nodes[0], ) + self.nodes[2:]

    @property
    def coords(self):
        return self.npos if self.mirror else self.ppos

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n3, n4 = self.dvnodes

        centre = (n3.pos + n4.pos) * 0.25 + n1.pos * 0.5

        yscale = 2 * 1.019 * self.scale
        if not self.mirror:
            yscale = -yscale

        # Note, scale scales by area, xscale and yscale scale by length.
        s = r'  \draw (%s) node[op amp, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            centre, self.args_str, 2 * 1.01 * self.scale, yscale,
            -self.angle, self.name)
        s += r'  \draw (%s.out) |- (%s);''\n' % (self.name, n1.s)
        s += r'  \draw (%s.+) |- (%s);''\n' % (self.name, n3.s)
        s += r'  \draw (%s.-) |- (%s);''\n' % (self.name, n4.s)
        # Draw label separately to avoid being scaled by 2.
        s += r'  \draw (%s) node[] {%s};''\n' % (centre, self.label(**kwargs))
        
        s += self._draw_nodes(**kwargs)
        return s


class FDOpamp(Cpt):

    npos = ((2.05, 1), (2.05, 0), (0, 0), (0, 1))
    ppos = ((2.05, -1), (2.05, 0), (0, 0), (0, -1))

    can_scale = True
    can_mirror = True

    @property
    def coords(self):
        return self.npos if self.mirror else self.ppos

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.dvnodes

        centre = self.tf((n1.pos + n2.pos + n3.pos + n4.pos) * 0.25, (0.15, 0))

        yscale = 2 * 1.02 * self.scale
        if not self.mirror:
            yscale = -yscale

        s = r'  \draw (%s) node[fd op amp, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            centre, self.args_str, 2 * 1.01 * self.scale, yscale,
            -self.angle, self.name)
        s += r'  \draw (%s.out +) |- (%s);''\n' % (self.name, n1.s)
        s += r'  \draw (%s.out -) |- (%s);''\n' % (self.name, n2.s)
        s += r'  \draw (%s.+) |- (%s);''\n' % (self.name, n3.s)
        s += r'  \draw (%s.-) |- (%s);''\n' % (self.name, n4.s)
        # Draw label separately to avoid being scaled by 2.
        s += r'  \draw (%s) node[] {%s};''\n' % (centre, self.label(**kwargs))
        
        s += self._draw_nodes(**kwargs)
        return s


class SPDT(Cpt):
    """SPDT switch"""

    can_stretch = False

    @property
    def coords(self):
        return ((0, 0.1565), (0.59, 0.313), (0.59, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3 = self.dvnodes

        centre = n1.pos * 0.5 + (n2.pos + n3.pos) * 0.25
        s = r'  \draw (%s) node[spdt, %s, rotate=%d] (%s) {};''\n' % (
            centre, self.args_str, self.angle, self.name)
        
        # TODO, fix label position.
        centre = (n1.pos + n3.pos) * 0.5 + Pos(0, -0.5)
        s += r'  \draw (%s) node[] {%s};''\n' % (centre, self.label(**kwargs))
        s += self._draw_nodes(**kwargs)
        return s


class Chip(Cpt):
    """General purpose chip"""

    can_stretch = False

    def __init__(self, sch, cpt_type, cpt_id, string, opts_string, nodes, *args):

        super (Chip, self).__init__(sch, cpt_type, cpt_id, string,
                                    opts_string, nodes, *args[2:])

        pins = []
        for node in self.nodes:
            pins.append(self.name + '@' + node)
        self.nodes = pins

    # TODO, tweak coord if pin name ends in \ using pinpos to
    # accomodate inverting circle.  This will require stripping of the
    # \ from the label. Alternatively, do not use inverting circle and
    # add overline to symbol name when printing.

    @property
    def centre(self):
        N = len(self.dvnodes)
        return self.midpoint(self.dvnodes[0], self.dvnodes[N // 2])

    @property
    def width(self):
        return self.w * self.size * self.sch.node_spacing

    @property
    def height(self):
        return self.h * self.size * self.sch.node_spacing

    def draw(self, **kwargs):

        if not self.check():
            return ''

        for m, n in enumerate(self.dvnodes):
            n.pinpos = self.pinpos[m]

        centre = self.centre

        w, h = self.width, self.height
        c1 = centre + Pos(-0.5 * w, 0.5 * h)
        c2 = centre + Pos(0.5 * w, 0.5 * h)
        c3 = centre + Pos(0.5 * w, -0.5 * h)
        c4 = centre + Pos(-0.5 * w, -0.5 * h)

        s = r'  \draw[thick] (%s) -- (%s) -- (%s) -- (%s) -- (%s);''\n' % (
            c1, c2, c3, c4, c1)
        s += r'  \draw (%s) node[] {%s};''\n'% (centre, self.label(**kwargs))

        s += self._draw_nodes(**kwargs)
        return s


class Chip1310(Chip):
    """Chip of size 1 3 1 0"""

    w = 1.5
    h = 1
    pinpos = ('l', 'b', 'b', 'b', 'r')

    @property
    def centre(self):
        return self.midpoint(self.dvnodes[0], self.dvnodes[4])

    @property
    def coords(self):
        w, h = self.width, self.height

        return ((-0.75, 0), (-0.5, -0.5), (0, -0.5), (0.5, -0.5),
                (0.75, 0))


class Chip2121(Chip):
    """Chip of size 2 1 2 1"""

    w = 2
    h = 2
    pinpos = ('l', 'l', 'b', 'r', 'r', 't')

    @property
    def coords(self):
        w, h = self.width, self.height

        return ((0, 0.5), (0, -0.5), (1, -1), 
                (2, -0.5), (2, 0.5), (1, 1))


class Chip3131(Chip):
    """Chip of size 3 1 3 1"""

    w = 2
    h = 3
    pinpos = ('l', 'l', 'l', 'b', 'r', 'r', 'r', 't')

    @property
    def coords(self):
        w, h = self.width, self.height

        return ((0, 1), (0, 0), (0, -1), (1, -1.5), 
                (2, -1), (2, 0), (2, 1), (1, 1.5))


class Chip4141(Chip):
    """Chip of size 4 1 4 1"""

    w = 2
    h = 4
    pinpos = ('l', 'l', 'l', 'l', 'b', 'r', 'r', 'r', 'r', 't')

    @property
    def coords(self):
        w, h = self.width, self.height

        return ((0, 1.5), (0, 0.5), (0, -0.5), (0, -1.5),
                (1, -2), 
                (2, -1.5), (2, -0.5), (2, 0.5), (2, 1.5),
                (1, 2))

class Ubuffer(Chip):
    """Buffer with power supplies"""

    can_scale = True

    pinpos = ('l', 'b', 'r', 't')

    @property
    def coords(self):
        return ((0, 0), (0.5, -0.25), (1.0, 0), (0.5, 0.25))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        for m, n in enumerate(self.dvnodes):
            n.pinpos = self.pinpos[m]

        n1, n2, n3, n4 = self.dvnodes
        centre = (n1.pos + n3.pos) * 0.5

        # TODO, create pgf shape
        q = self.tf(centre, ((-1, 0), (1, 0), (0, 0.5), (0, -0.5),
                             (-1, 1), (-1, -1)))

        s = r'  \draw[thick] (%s) -- (%s) -- (%s) -- (%s);''\n' % (
            q[4], q[1], q[5], q[4])
        s += r'  \draw (%s) node[] (%s) {%s};''\n' % (
            centre, self.name, self.label(**kwargs))
        s += self._draw_nodes(**kwargs)
        return s

class Uinverter(Chip):
    """Inverter with power supplies"""

    can_scale = True

    pinpos = ('l', 'b', 'r', 't')

    @property
    def coords(self):
        return ((0, 0), (0.5, -0.22), (1.0, 0), (0.5, 0.22))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        for m, n in enumerate(self.dvnodes):
            n.pinpos = self.pinpos[m]

        n1, n2, n3, n4 = self.dvnodes
        centre = (n1.pos + n3.pos) * 0.5

        # TODO, create pgf shape
        w = 0.1

        q = self.tf(centre, ((-1, 1), (-1, -1), (1 - 2 * w, 0), (1 - w, 0)))

        s = r'  \draw[thick] (%s) -- (%s) -- (%s) -- (%s);''\n' % (
            q[0], q[2], q[1], q[0])
        s += r'  \draw[thick] (%s) node[ocirc, scale=%s] {};''\n' % (q[3], 1.8 * self.scale)
        s += r'  \draw (%s) node[] (%s) {%s};''\n' % (
            centre, self.name, self.label(**kwargs))
        s += self._draw_nodes(**kwargs)
        return s

class Wire(OnePort):

    def __init__(self, sch, cpt_type, cpt_id, string, opts_string, nodes, *args):

        super (Wire, self).__init__(sch, cpt_type, cpt_id, string,
                                    opts_string, nodes, *args)

        if self.implicit:
            # Rename second node since this is spatially different from
            # other nodes of the same name.  Add underscore so node
            # not drawn.
            self.nodes = (self.nodes[0], self.name + '@_' + self.nodes[1])

    @property
    def implicit(self):

        for key in self.implicit_keys:
            if key in self.opts:
                return True
        return False

    @property
    def coords(self):

        return ((0, 0), (1, 0))

    @property
    def vnodes(self):

        return self.nodes

    def draw_implicit(self, **kwargs):
        """Draw implict wires, i.e., connections to ground, etc."""

        kind = None
        for key in self.implicit_keys:
            if key in self.opts:
                kind = key
                break;

        # I like the sground symbol for power supplies but rground symbol
        # is also common.
        if kind == 'implicit':
            kind = 'sground'
        anchor = 'south west'
        if self.down:
            anchor = 'north west'

        n1, n2 = self.dvnodes
        s = r'  \draw (%s) -- (%s);''\n' % (n1.s, n2.s)
        s += r'  \draw (%s) node[%s,scale=0.5,rotate=%d] {};''\n' % (
            n2.s, kind, self.angle + 90)

        if 'l' in self.opts:
            lpos = self.tf(n2.pos, (0.25, 0))
            s += r'  \draw [anchor=%s] (%s) node {%s};''\n' % (
                anchor, lpos, self.label(**kwargs))
        return s

    def draw(self, **kwargs):

        if not self.check():
            return ''

        if self.implicit:
            return self.draw_implicit(**kwargs)
            
        def arrow_map(name):

            try:
                return {'tee': '|', 'otri' : 'open triangle 60',
                        'tri' : 'triangle 60'}[name]
            except:
                return name

        n1, n2 = self.dvnodes

        # W 1 2; up, arrow=tri, l=V_{dd}
        # W 1 3; right, arrow=otri
        # W 1 4; down, arrow=tee, l=0V
        # W 1 5; left, startarrow=tri, endarrow=open triangle 90, bus=8

        startarrow = self.opts.pop('startarrow', '')
        endarrow = self.opts.pop('arrow', '')
        endarrow = self.opts.pop('endarrow', endarrow)

        bus = self.opts.pop('bus', False)
        style = ''
        if bus:
            # TODO if bus has numeric arg, indicate number of lines with slash.
            style = 'ultra thick'

        # TODO, add arrow shapes for earth symbol.

        s = r'  \draw[%s-%s, %s, %s] (%s) to (%s);''\n' % (
            arrow_map(startarrow), arrow_map(endarrow), style,
            self.args_str, n1.s, n2.s)
        s += self._draw_nodes(**kwargs)

        if self.voltage_str != '':
            print('There is no voltage drop across an ideal wire!')

        if self.current_str != '' or 'l' in self.opts:
            # FIXME, we don't really want the wire drawn since this
            # can clobber the arrow.  We just want the current
            # annotation and/or the label.
            s += r'  \draw[%s] (%s) [short, %s, %s] to (%s);''\n' % (
                self.args_str, n1.s, self.current_str, self.label_str, n2.s)
        return s


class XT(Cpt):
    """Crystal"""

    can_scale = True

    @property
    def coords(self):
        return ((-0.5, 0), (0.5, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.dvnodes

        centre = (n1.pos + n2.pos) * 0.5
        q = self.xtf(centre, ((-0.3, 0), (-0.3, 0.3), (-0.3, -0.3),
                              (0.3, 0), (0.3, 0.3), (0.3, -0.3),
                              (-0.12, 0.3), (0.12, 0.3),
                              (0.12, -0.3), (-0.12, -0.3),
                              (0.0, -0.6)))

        s = r'  \draw[thick] (%s) -- (%s);''\n' % (q[1], q[2])
        s += r'  \draw[thick] (%s) -- (%s);''\n' % (q[4], q[5])
        s += r'  \draw[thick] (%s) -- (%s) -- (%s) -- (%s) -- (%s);''\n' % (
            q[6], q[7], q[8], q[9], q[6])
        s += r'  \draw (%s) -- (%s);''\n' % (q[0], n1.s)
        s += r'  \draw (%s) -- (%s);''\n' % (q[3], n2.s)
        s += r'  \draw (%s) node[] {%s};''\n'% (q[10], self.label(**kwargs))
        s += self._draw_nodes(**kwargs)
        return s


classes = {}

def defcpt(name, base, docstring, cpt=None):
    
    if isinstance(base, str):
        base = classes[base]

    newclass = type(name, (base, ), {'__doc__': docstring})

    if cpt is not None:
        newclass.tikz_cpt = cpt
    classes[name] = newclass


def make(classname, parent, cpt_type, cpt_id,
         string, opts_string, nodes, *args):

    # Create instance of component object
    try:
        newclass = getattr(module, classname)
    except:
        newclass = classes[classname]

    cpt = newclass(parent, cpt_type, cpt_id, string, opts_string, 
                   nodes, *args)
    # Add named attributes for the args?   Lname1, etc.
        
    return cpt

# Dynamically create classes.
defcpt('AM', OnePort, 'Ammeter', 'ammeter')

defcpt('BAT', OnePort, 'Battery', 'battery')

defcpt('C', OnePort, 'Capacitor', 'C')

defcpt('D', OnePort, 'Diode', 'D')
defcpt('Dled', 'D', 'LED', 'leD')
defcpt('Dphoto', 'D', 'Photo diode', 'pD')
defcpt('Dschottky', 'D', 'Schottky diode', 'zD')
defcpt('Dtunnel', 'D', 'Tunnel diode', 'tD')
defcpt('Dzener', 'D', 'Zener diode', 'zD')

defcpt('E', VCS, 'VCVS', 'american controlled voltage source')
defcpt('Eopamp', Opamp, 'Opamp')
defcpt('Efdopamp', FDOpamp, 'Fully differential opamp')
defcpt('F', VCS, 'CCCS', 'american controlled current source')
defcpt('G', CCS, 'VCCS', 'american controlled current source')
defcpt('H', CCS, 'CCVS', 'american controlled voltage source')

defcpt('I', OnePort, 'Current source', 'I')
defcpt('sI', OnePort, 'Current source', 'I')
defcpt('Isin', 'I', 'Sinusoidal current source', 'sI')
defcpt('Idc', 'I', 'DC current source', 'I')
defcpt('Iac', 'I', 'AC current source', 'sI')

defcpt('J', JFET, 'N JFET transistor', 'njfet')
defcpt('Jnjf', 'J', 'N JFET transistor', 'njfet')
defcpt('Jpjf', 'J', 'P JFET transistor', 'pjfet')

defcpt('L', OnePort, 'Inductor', 'L')

defcpt('M', MOSFET, 'N MOSJFET transistor', 'nmos')
defcpt('Mnmos', 'M', 'N channel MOSJFET transistor', 'nmos')
defcpt('Mpmos', 'M', 'P channel MOSJFET transistor', 'pmos')

defcpt('O', OnePort, 'Open circuit', 'open')
defcpt('P', OnePort, 'Port', 'open')

defcpt('Q', Transistor, 'NPN transistor', 'npn')
defcpt('Qpnp', 'Q', 'PNP transistor', 'pnp')
defcpt('Qnpn', 'Q', 'NPN transistor', 'npn')

defcpt('R', OnePort, 'Resistor', 'R')

defcpt('SW', OnePort, 'Switch', 'closing switch')
defcpt('SWno', 'SW', 'Normally open switch', 'closing switch')
defcpt('SWnc', 'SW', 'Normally closed switch', 'opening switch')
defcpt('SWpush', 'SW', 'Pushbutton switch', 'push button')
defcpt('SWspdt', SPDT, 'SPDT switch', 'spdt')

defcpt('TF', Transformer, 'Transformer', 'transformer')
defcpt('TFcore', Transformer, 'Transformer with core', 'transformer core')
defcpt('TFtapcore', TFtap, 'Tapped transformer with core', 'transformer core')
defcpt('TP', TwoPort, 'Two port', '')

defcpt('Uchip1310', Chip1310, 'General purpose chip')
defcpt('Uchip2121', Chip2121, 'General purpose chip')
defcpt('Uchip3131', Chip3131, 'General purpose chip')
defcpt('Uchip4141', Chip4141, 'General purpose chip')

defcpt('V', OnePort, 'Voltage source', 'V')
defcpt('sV', OnePort, 'Voltage source', 'V')
defcpt('Vsin', 'V', 'Sinusoidal voltage source', 'sV')
defcpt('Vdc', 'V', 'DC voltage source', 'V')
defcpt('Vstep', 'V', 'Step voltage source', 'V')
defcpt('Vac', 'V', 'AC voltage source', 'sV')

defcpt('VM', OnePort, 'Voltmeter', 'voltmeter')

defcpt('W', Wire, 'Wire', 'short')
defcpt('Y', OnePort, 'Admittance', 'european resistor')
defcpt('Z', OnePort, 'Impedance', 'european resistor')

# Perhaps AM for ammeter, VM for voltmeter, VR for variable resistor?
# Currently, a variable resistor is supported with the variable
# option.
