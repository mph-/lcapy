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

# There are two types of component (Cpt).
#
# 1. The stretchable components (such as resistors and capacitors) have
# wires that can be stretched.  The size attribute controls the
# spacing between the nodes but does not affect the component size.
# The component size can be changed with the scale attribute (this changes
# the dipole length).  The aspect ratio is not easy to change (need to use
# dipole/resistor/height).
#
# 2.  The fixed components (such as ICs) do not have wires and cannot
# be stretched.  The (unrotated) width is set by the size attribute and
# the (unrotated) height is set by the aspect attribute.  The scale attribute
# is not used.

# There are two paradigms used for specifying node coordinates:
#
# 1.  The old model.  required_node_names returns subset of
# explicit_node_names as a list.
#
# 2.  The new model.  node_map specifies the subset of required nodes.


class Cpt(object):

    voltage_keys = ('v', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<',
                    'v<', 'v>')
    current_keys = ('i', 'i_', 'i^', 'i_>',  'i_<', 'i^>', 'i^<',
                    'i>_', 'i<_', 'i>^', 'i<^', 'i>', 'i<')
    label_keys = ('l', 'l_', 'l^')
    implicit_keys =  ('implicit', 'ground', 'sground', 'rground')
    # The following keys do not get passed through to circuitikz.
    misc_keys = ('left', 'right', 'up', 'down', 'rotate', 'size',
                 'mirror', 'scale', 'invisible', 'variable', 'fixed',
                 'aspect', 'pins', 'image')

    can_rotate = True
    can_scale = False
    can_mirror = False
    can_stretch = True
    default_width = 1.0
    default_aspect = 1.0
    node_map = ()
    required_anchors = ()
    anchors = {}

    @property
    def s(self):
        """Sanitised name"""
        return self.name.replace('.', '@')

    def __init__(self, sch, name, cpt_type, cpt_id, string,
                 opts_string, node_names, keyword, *args):

        self.sch = sch
        self.type = cpt_type
        self.id = cpt_id
        self.name = name

        self.net = string.split(';')[0]
        self.opts_string = opts_string

        self.args = args
        self.classname = self.__class__.__name__

        # Drawing hints
        self.opts = Opts(opts_string)

        self.explicit_node_names = node_names

        extra_node_names = []
        for anchor in self.required_anchors:
            extra_node_names.append(name + '.' + anchor)
        self.extra_node_names = tuple(extra_node_names)
        
        anchor_node_names = []
        for anchor in self.anchors.keys():
            anchor_node_names.append(name + '.' + anchor)
        self.anchor_node_names = tuple(anchor_node_names)

        self.node_names = self.required_node_names + self.anchor_node_names

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
            val = self.default_width
        if val == '':
            val = self.default_width
        return float(val)

    @property
    def scale(self):
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
            angle = -90 if self.type in ('P', ) else 0
        
        if 'rotate' in self.opts:
            angle += float(self.opts['rotate'])
        return angle

    @property
    def w(self):
        """Normalised width"""
        return 1.0

    @property
    def h(self):
        """Normalised height"""
        return self.w / self.aspect

    @property
    def aspect(self):
        return float(self.opts.get('aspect', self.default_aspect))

    def R(self, angle_offset=0):
        """Return rotation matrix"""
        angle = self.angle + angle_offset
        
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
    def required_node_names(self):
        """Subset of explicit_node_names"""

        # Old model.   The number of node names in the list
        # must match the number of entries in coords.
        if self.node_map == ():
            return self.explicit_node_names

        # New model.
        node_names = []
        for anchor, node_name in zip(self.node_map, self.explicit_node_names):
            if anchor != '':
                node_names.append(node_name)
        return tuple(node_names)

    def node(self, anchor):
        """Return node by anchor"""

        if anchor in self.node_map:
            index = self.node_map.index(anchor)
            node_name = self.explicit_node_names[index]
        else:
            node_name = self.name + '.' + anchor
        for node in self.nodes:
            if node.name == node_name:
                return node
        raise ValueError('Unknown anchor %s for %s' % (anchor, self))

    @property
    def nodes(self):
        """Nodes used to draw the current element."""

        if hasattr(self, '_nodes'):
            return self._nodes

        # Perhaps determine coords here as well and cache them?
        
        node_names = self.required_node_names + self.extra_node_names + self.anchor_node_names
        
        rnodes = []
        for n in node_names:
            if n in self.sch.nodes:
                rnodes.append(self.sch.nodes[n])
        self._nodes = rnodes
        return rnodes

    @property
    def drawn_nodes(self):
        """Nodes that are drawn"""
        return self.nodes

    @property
    def coords(self):

        # Determine the required coords.
        rcoords = []
        for node in self.nodes:
            node_name = node.name
            if node_name in self.explicit_node_names:
                index = self.explicit_node_names.index(node_name)
                anchor = self.node_map[index]
            elif node_name in self.sch.nodes:
                anchor = node_name.split('.')[-1]
            else:
                raise ValueError('Unknown node %s' % node_name)

            if anchor != '':
                rcoords.append(self.anchors[anchor])
        return rcoords

    @property
    def scoords(self):
        """Scaled coordinates for each of the nodes"""

        c = np.array(self.coords)
        S = np.array(((self.w, self.h),) * c.shape[0])
        return c * S

    @property
    def tcoords(self):
        """Transformed coordinates for each of the nodes"""
        if hasattr(self, '_tcoords'):
            return self._tcoords

        self._tcoords = np.dot(self.scoords, self.R())
        return self._tcoords

    @property
    def xvals(self):
        return self.tcoords[:, 0]

    @property
    def yvals(self):
        return self.tcoords[:, 1]

    def xlink(self, graphs):

        xvals = self.xvals
        for m1, n1 in enumerate(self.nodes):
            for m2, n2 in enumerate(self.nodes[m1 + 1:], m1 + 1):
                if xvals[m2] == xvals[m1]:
                    graphs.link(n1.name, n2.name)

    def ylink(self, graphs):

        yvals = self.yvals
        for m1, n1 in enumerate(self.nodes):
            for m2, n2 in enumerate(self.nodes[m1 + 1:], m1 + 1):
                if yvals[m2] == yvals[m1]:
                    graphs.link(n1.name, n2.name)

    def place(self, graphs, vals):
        
        size = self.size
        idx = np.argsort(vals)[::-1]
        for i in range(len(idx) - 1):
            m1 = idx[i]
            m2 = idx[i + 1]
            n1 = self.nodes[m1]
            n2 = self.nodes[m2]
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

    def draw_node(self, n, **kwargs):

        draw_nodes = kwargs.get('draw_nodes', True)        

        s = ''
        if not draw_nodes:
            return s

        if not n.visible(draw_nodes) or n.pin:
            return s

        if n.port:
            s = r'  \draw (%s) node[ocirc] {};''\n' % n.s
        else:
            s = r'  \draw (%s) node[circ] {};''\n' % n.s

        return s

    def draw_nodes(self, **kwargs):

        s = ''
        for n in self.drawn_nodes:
            s += self.draw_node(n, **kwargs)
        return s

    def draw(self, **kwargs):
        raise NotImplementedError('draw method not implemented for %s' % self)


    def opts_str_list(self, choices):
        """Format voltage, current, or label string as a key-value pair
        and return list of strings"""

        def fmt(key, val):
            label = format_label(val)
            if label == '':
                label = '{}'
            if not (label[0] == '{' and label[-1] == '}'):
                label = '{' + label + '}'

            return '%s=%s' % (key, label)

        return [fmt(key, val) for key, val in self.opts.items()
                if key in choices]

    def opts_str(self, choices):
        """Format voltage, current, or label string as a key-value pair"""

        return ','.join(self.opts_str_list(choices))

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
    def label_str_list(self):

        return self.opts_str_list(self.label_keys)

    @property
    def args_str(self):

        def fmt(key, val):
            return '%s=%s' % (key, format_label(val))

        return ','.join([fmt(key, val) 
                         for key, val in self.opts.items()
                         if key not in self.voltage_keys + self.current_keys + self.label_keys + self.misc_keys + self.implicit_keys])

    def label(self, **kwargs):

        label_values = kwargs.get('label_values', True)
        label_ids = kwargs.get('label_ids', True)        

        label_str = ''
        if label_ids:
            label_str = self.id_label
        if label_values and self.value_label != '':
            label_str = self.value_label        

        # Override label if specified.  There are no placement options.
        string =  ','.join([format_label(val)
                            for key, val in self.opts.items()
                            if key in ('l', )])

        if string != '':
            label_str = string
        # Remove curly braces.
        if len(label_str) > 1 and label_str[0] == '{' and label_str[-1] == '}':
            label_str = label_str[1:-1]
        return label_str

    def label_make(self, label_pos='', **kwargs):

        # TODO merge with label

        label_values = kwargs.get('label_values', True)
        label_ids = kwargs.get('label_ids', True)

        # Generate default label.
        if (label_ids and label_values and self.id_label != '' 
            and self.value_label and self.id_label != self.value_label):
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

    def tf(self, centre, offset, angle_offset=0.0):
        """Transform coordinate"""

        raise NotImplementedError('tf method not implemented for %s' % self)

    def draw_path(self, points, style='', join='--', closed=False):

        path = (' %s ' % join).join(['(%s)' % point for point in points])
        if closed:
            path += ' %s cycle' % join

        args_str = self.args_str
        if style == '':
            s = args_str
        elif args_str == '':
            s = style
        else:
            s = style + ', ' + args_str
        if s != '':
            s = '[%s]' % s

        return r'  \draw%s %s;''\n' % (s, path)

    def draw_label(self, pos, **kwargs):

        return r'  \draw[%s] (%s) node[] {%s};''\n'% (
            self.args_str, pos, self.label(**kwargs))


class StretchyCpt(Cpt):

    can_stretch = True

    def xtf(self, centre, offset, angle_offset=0.0):
        """Transform coordinate."""

        # Note the size attribute is not used.   This only scales the x coords.
        if isinstance(offset[0], tuple):
            return [self.xtf(centre, offset1, angle_offset) for offset1 in offset]

        return centre + np.dot((offset[0] * self.w * self.scale, offset[1] * self.h), self.R(angle_offset)) * self.sch.node_spacing


    def tf(self, centre, offset, angle_offset=0.0):
        """Transform coordinate."""

        # Note the size attribute is not used.
        if isinstance(offset[0], tuple):
            return [self.tf(centre, offset1, angle_offset) for offset1 in offset]

        return centre + np.dot((offset[0] * self.w, offset[1] * self.h), self.R(angle_offset)) * self.scale * self.sch.node_spacing


class FixedCpt(Cpt):

    can_stretch = False

    def tf(self, centre, offset, angle_offset=0.0):
        """Transform coordinate."""

        if isinstance(offset[0], tuple):
            return [self.tf(centre, offset1, angle_offset) for offset1 in offset]

        return centre + np.dot((offset[0] * self.w, offset[1] * self.h), self.R(angle_offset)) * self.size * self.scale * self.sch.node_spacing


class Transistor(FixedCpt):
    """Transistor"""
    
    npos = ((1, 1.5), (0, 0.75), (1, 0))
    ppos = ((1, 0), (0, 0.75), (1, 1.5))

    can_mirror = True
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

        n1, n2, n3 = self.nodes
        centre = (n1.pos + n3.pos) * 0.5

        s = r'  \draw (%s) node[%s, %s, scale=%s, rotate=%d] (%s) {};''\n' % (
            centre, self.tikz_cpt, self.args_str, 2 * self.scale,
            self.angle, self.s)
        s += self.draw_label(centre, **kwargs)

        # Add additional wires.  These help to compensate for the
        # slight differences in sizes of the different transistors.
        if self.tikz_cpt in ('pnp', 'npn'):
            s += r'  \draw (%s.C) -- (%s) (%s.B) -- (%s) (%s.E) -- (%s);''\n' % (
                self.s, n1.s, self.s, n2.s, self.s, n3.s)
        else:
            s += r'  \draw (%s.D) -- (%s) (%s.G) -- (%s) (%s.S) -- (%s);''\n' % (
                self.s, n1.s, self.s, n2.s, self.s, n3.s)

        s += self.draw_nodes(**kwargs)
        return s


class JFET(Transistor):
    """Transistor"""

    #node_map = ('c', 'b', 'e')
    npos = ((1, 1.5), (0, 0.48), (1, 0))
    ppos = ((1, 0), (0, 1.02), (1, 1.5))


class MOSFET(Transistor):
    """Transistor"""

    #node_map = ('d', 'g', 's')    
    npos = ((0.85, 1.52), (-0.25, 0.76), (0.85, 0))
    ppos = ((0.85, 0), (-0.25, 0.76), (0.85, 1.52))


class TwoPort(FixedCpt):
    """Two-port"""

    # TODO
    can_rotate = False

    @property
    def coords(self):
        return ((1.5, 1), (1.5, 0), (0, 1), (0, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.nodes
        width = n2.pos.x - n4.pos.x
        centre = (n1.pos + n2.pos + n3.pos + n4.pos) * 0.25

        q = self.tf(centre, ((-1.5, -1.5), (-1.5, 1.5), (1.5, 1.5), (1.5, -1.5),
                             (0, 1.15)))

        top = q[4]

        titlestr = ''
        if len(self.args) > 0:
            titlestr = "%s-parameter two-port" % self.args[0]

        s = self.draw_path(q[0:4], closed=True)
        s += r'  \draw (%s) node[text width=%.1fcm, align=center] (%s) {%s};''\n' % (
            centre, width, titlestr, self.s)
        s += r'  \draw (%s) node[text width=%.1fcm, align=center, %s] {%s};''\n' % (
            top, width, self.args_str, self.label(**kwargs))

        s += self.draw_nodes(**kwargs)
        return s


class MX(FixedCpt):
    """Mixer"""

    # Dubious
    can_scale = True

    @property
    def coords(self):
        return ((0.25, 0.25), (-0.25, 0.25), (0, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3 = self.nodes

        centre = (n1.pos + n2.pos) * 0.5
        q = self.tf(centre, ((0, 0.35)))

        s = r'  \draw (%s) node[mixer,xscale=%s] {};''\n' % (
            centre, self.scale * self.size)
        s += self.draw_label(q, **kwargs)
        return s


class SP(FixedCpt):
    """Summing point"""

    # Dubious
    can_scale = True
    can_mirror = True

    @property
    def coords(self):
        if self.mirror:
            return ((-0.25, 0), (0, 0.25), (0.25, 0), (0, -0.25))
        return ((-0.25, 0), (0, -0.25), (0.25, 0), (0, 0.25))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n = self.nodes

        centre = (n[0].pos + n[2].pos) * 0.5
        q = self.tf(centre, ((0.3, 0.3), (-0.125, 0), (0, -0.125),
                             (0, 0.125), (0, 0.125)))
        xscale = self.scale * self.size
        yscale = self.scale * self.size       
        if self.mirror:
            yscale = -yscale

        s = r'  \draw (%s) node[mixer, xscale=%s, yscale=%s, rotate=%s] {};''\n' % (
            centre, xscale, yscale, self.angle)
        s += self.draw_label(q[0], **kwargs)
        s += r'  \draw (%s) node[] {$%s$};''\n'% (q[1], self.labels[0])

        if self.mirror:
            s += r'  \draw (%s) node[] {$%s$};''\n'% (q[4], self.labels[0])
        else:
            s += r'  \draw (%s) node[] {$%s$};''\n'% (q[2], self.labels[1])
        if len(self.labels) > 2:
            s += r'  \draw (%s) node[] {$%s$};''\n'% (q[3], self.labels[2])
        return s


class SP3(SP):
    """Summing point"""

    @property
    def coords(self):
        if self.mirror:
            return ((-0.25, 0), (0, 0.25), (0.25, 0))
        return ((-0.25, 0), (0, -0.25), (0.25, 0))


class SPpp(SP3):
    """Summing point"""

    labels = '++'


class SPpm(SP3):
    """Summing point"""

    labels = '+-'


class SPppp(SP):
    """Summing point"""

    labels = '+++'


class SPpmm(SP):
    """Summing point"""

    labels = '+--'


class SPppm(SP):
    """Summing point"""

    labels = '++-'


class TL(StretchyCpt):
    """Transmission line"""

    # Dubious.  Perhaps should stretch this component in proportion to size?
    can_scale = True

    @property
    def coords(self):
        return ((1.25, 0.5), (1.25, 0), (0, 0.5), (0, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.nodes

        centre = (n1.pos + n3.pos) * 0.5
        q = self.xtf(centre, ((0.32, 0),
                              (0.27, -0.145),
                              (-0.35, 0),
                              (-0.35, -0.145),
                              (-0.65, 0)))

        # Rotation creates an ellipse!
        s = r'  \draw (%s) node[tlinestub, xscale=%s, rotate=%s] {};''\n' % (
            q[4], self.scale, self.angle)
        s += self.draw_label(centre, **kwargs)

        s += self.draw_path((q[0], n1.s))
        s += self.draw_path((q[1], n2.s), join='|-')
        s += self.draw_path((q[2], n3.s))
        s += self.draw_path((q[3], n4.s), join='|-')
        s += self.draw_nodes(**kwargs)
        return s


class TF1(FixedCpt):
    """Transformer"""

    @property
    def coords(self):
        return ((0.5, 1), (0.5, 0), (0, 1), (0, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        link = kwargs.get('link', True)

        p = [node.pos for node in self.nodes]

        centre = (p[0] + p[1] + p[2] + p[3]) * 0.25
        q = self.tf(centre, ((-0.35, 0.3), (0.35, 0.3), (0, 0.375)))
        primary_dot = q[0]
        secondary_dot = q[1]
        labelpos = q[2]

        s = r'  \draw (%s) node[circ] {};''\n' % primary_dot
        s += r'  \draw (%s) node[circ] {};''\n' % secondary_dot
        s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
            labelpos, 0.5, self.s, self.label(**kwargs))

        if link:
            # TODO: allow for rotation
            width = p[0].x - p[2].x
            arcpos = Pos((p[0].x + p[2].x) / 2, secondary_dot.y - width / 2 + 0.3)

            s += r'  \draw [<->] ([shift=(45:%.2f)]%s) arc(45:135:%.2f);''\n' % (
                width / 2, arcpos, width / 2)

        if self.classname in ('TFcore', 'TFtapcore'):
            # Draw core
            q = self.tf(centre, ((-0.05, -0.2), (-0.05, 0.2),
                                 (0.05, -0.2), (0.05, 0.2)))
            s += self.draw_path(q[0:2], style='thick')
            s += self.draw_path(q[2:4], style='thick')

        return s


class Transformer(TF1):
    """Transformer"""

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.nodes

        s = r'  \draw (%s) to [inductor] (%s);''\n' % (n3.s, n4.s)
        s += r'  \draw (%s) to [inductor] (%s);''\n' % (n2.s, n1.s)

        s += super(Transformer, self).draw(link=False, **kwargs)
        s += self.draw_nodes(**kwargs)
        return s


class TFtap(TF1):
    """Transformer"""

    @property
    def coords(self):
        return ((0.5, 1), (0.5, 0), (0, 1), (0, 0), (-0.125, 0.55), (0.625, 0.55))

    @property
    def drawn_nodes(self):
        # Do not draw the taps.
        return self.nodes[0:4]

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4, n5, n6 = self.nodes

        s = r'  \draw (%s) to [inductor] (%s);''\n' % (n3.s, n4.s)
        s += r'  \draw (%s) to [inductor] (%s);''\n' % (n2.s, n1.s)

        s += super(TFtap, self).draw(link=False, **kwargs)
        s += self.draw_nodes(**kwargs)
        return s


class K(TF1):
    """Mutual coupling"""

    def __init__(self, sch, name, cpt_type, cpt_id, string,
                 opts_string, node_names, keyword, *args):

        self.Lname1 = args[0]
        self.Lname2 = args[1]
        super (K, self).__init__(sch, name, cpt_type, cpt_id, string,
                                 opts_string, node_names, keyword, *args[2:])

    @property
    def nodes(self):

        # CHECKME
        
        # L1 and L2 need to be previously defined so we can find their nodes.
        L1 = self.sch.elements[self.Lname1]
        L2 = self.sch.elements[self.Lname2]
        return [self.sch.nodes[n] for n in L1.node_names + L2.node_names]


class OnePort(StretchyCpt):
    """OnePort"""

    can_mirror = True
    can_scale = True

    @property
    def coords(self):
        return ((0, 0), (1, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.nodes

        tikz_cpt = self.tikz_cpt
        if self.variable:
            if self.type in ('C', 'R', 'L'):
                tikz_cpt = 'v' + tikz_cpt
            else:
                raise Error('Component %s not variable' % self.name)

        label_pos = '_'
        voltage_pos = '^'
        if (self.type in ('V', 'I', 'E', 'F', 'G', 'H', 'BAT') and
            self.sch.circuitikz_version < '2016/01/01'):

            # Old versions of circuitikz expects the positive node
            # first, except for voltage and current sources!  So
            # swap the nodes otherwise they are drawn the wrong
            # way around.
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

        args_str = self.args_str
        args_str2 = ','.join([self.voltage_str, self.current_str])

        if self.mirror:
            args_str += ',mirror'

        if self.scale != 1.0:
            args_str2 += ',bipoles/length=%scm' % (self.sch.cpt_size * self.scale)

        label_str = self.label_make(label_pos, **kwargs)
            
        s = r'  \draw[%s] (%s) to [%s,%s,%s,%s,n=%s] (%s);''\n' % (
            args_str, n1.s, tikz_cpt, label_str, args_str2,
            node_pair_str, self.s, n2.s)
        return s


class VCS(OnePort):
    """Voltage controlled source"""

    @property
    def required_node_names(self):
        return self.explicit_node_names[0:2]


class CCS(OnePort):
    """Current controlled source"""
    pass


class Opamp(FixedCpt):

    can_scale = True
    can_mirror = True

    required_anchors = ('mid', )

    # The Nm node is not used (ground).
    node_map = ('out', '', '+', '-')
    
    panchors = {'out' : (2.4, 0.0),
                '+' : (0.0, 0.5),
                '-' : (0.0, -0.5),
                'mid' : (1.2, 0.0),
                'vdd' : (1.2, 0.5),
                'vdd2' : (0.8, 0.745),
                'vss2' : (0.8, -0.745),
                'vss' : (1.2, -0.5),
                'ref' : (1.6, -0.255),
                'r+' : (0.35, 0.25),
                'r-' : (0.35, -0.25)}
    
    nanchors = {'out' : (2.4, 0.0),
                '+' : (0.0, -0.5),
                '-' : (0.0, 0.5),
                'mid' : (1.2, 0.0),
                'vdd' : (1.2, 0.5),
                'vdd2' : (0.8, 0.745),
                'vss2' : (0.8, -0.745),
                'vss' : (1.2, -0.5),
                'ref' : (1.6, -0.255),
                'r+' : (0.35, 0.25),
                'r-' : (0.35, -0.25)}

    @property
    def anchors(self):
        return self.nanchors if self.mirror else self.panchors
    
    def draw(self, **kwargs):

        if not self.check():
            return ''

        yscale = 2 * 1.019 * self.scale
        if not self.mirror:
            yscale = -yscale

        centre = self.node('mid')

        # Note, scale scales by area, xscale and yscale scale by length.
        s = r'  \draw (%s) node[op amp, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            centre.s,
            self.args_str, 2 * 1.01 * self.scale, yscale,
            -self.angle, self.s)
        s += r'  \draw (%s.out) |- (%s);''\n' % (self.s, self.node('out').s)
        s += r'  \draw (%s.+) |- (%s);''\n' % (self.s, self.node('+').s)
        s += r'  \draw (%s.-) |- (%s);''\n' % (self.s, self.node('-').s)
        s += self.draw_label(centre.s, **kwargs)
        s += self.draw_nodes(**kwargs)
        return s


class FDOpamp(FixedCpt):

    can_scale = True
    can_mirror = True

    required_anchors = ('mid', )

    node_map = ('out+', 'out-', '+', '-')
    
    panchors = {'out+' : (2.05, -0.5),
                'out-' : (2.05, 0.5),                
                '+' : (0.0, 0.5),
                '-' : (0.0, -0.5),
                'mid' : (1.2, 0.0),
                'vdd' : (1.0, 0.62),
                'vss' : (1.0, -0.62),
                'r+' : (0.35, 0.25),
                'r-' : (0.35, -0.25)}
    
    nanchors = {'out+' : (2.05, 0.5),
                'out-' : (2.05, -0.5),
                '+' : (0.0, -0.5),
                '-' : (0.0, 0.5),
                'mid' : (1.2, 0.0),
                'vdd' : (1.0, 0.62),
                'vss' : (1.0, -0.62),
                'r+' : (0.35, 0.25),
                'r-' : (0.35, -0.25)}

    @property
    def anchors(self):
        return self.nanchors if self.mirror else self.panchors
    
    def draw(self, **kwargs):

        if not self.check():
            return ''

        centre = self.node('mid')

        yscale = 2 * 1.02 * self.scale
        if not self.mirror:
            yscale = -yscale

        s = r'  \draw (%s) node[fd op amp, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            centre.s, self.args_str, 2 * 1.01 * self.scale, yscale,
            -self.angle, self.s)
        s += r'  \draw (%s.out +) |- (%s);''\n' % (self.s, self.node('out+').s)
        s += r'  \draw (%s.out -) |- (%s);''\n' % (self.s, self.node('out-').s)
        s += r'  \draw (%s.+) |- (%s);''\n' % (self.s, self.node('+').s)
        s += r'  \draw (%s.-) |- (%s);''\n' % (self.s, self.node('-').s)
        s += self.draw_label(centre.s, **kwargs)
        s += self.draw_nodes(**kwargs)
        return s


class SPDT(StretchyCpt):
    """SPDT switch"""

    @property
    def coords(self):
        return ((0, 0.1565), (0.59, 0.313), (0.59, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3 = self.nodes

        centre = n1.pos * 0.5 + (n2.pos + n3.pos) * 0.25
        s = r'  \draw (%s) node[spdt, %s, rotate=%d] (%s) {};''\n' % (
            centre, self.args_str, self.angle, self.s)
        
        # TODO, fix label position.
        centre = (n1.pos + n3.pos) * 0.5 + Pos(0, -0.5)
        s += self.draw_label(centre, **kwargs)
        s += self.draw_nodes(**kwargs)
        return s


class Shape(FixedCpt):
    """General purpose shape"""

    default_aspect = 1.0
    can_mirror = True


    @property
    def centre(self):
        if hasattr(self, 'anchors'):
            # Look for centre anchor.
            for node in self.nodes:
                if node.name.split('.')[-1] == 'mid':
                    return node.pos
        
        N = len(self.nodes)
        return self.midpoint(self.nodes[0], self.nodes[N // 2])

    @property
    def width(self):
        return self.w * self.size * self.sch.node_spacing

    @property
    def height(self):
        return self.h * self.size * self.sch.node_spacing

    def draw(self, **kwargs):

        if not self.check():
            return ''

        label = self.label(**kwargs)
        if 'image' in self.opts:
            # Override label with image
            label = r'\includegraphics[width=%scm]{%s}' % (self.width - 0.5,
                                                           self.opts['image'])

        text_width = self.width * 0.8

        # shape border rotate rotates the box but not the text
        s = r'  \draw (%s) node[%s, thick, inner sep=0pt, minimum width=%scm, minimum height=%scm, text width=%scm, align=center, shape border rotate=%s, draw, %s] (%s) {%s};''\n'% (
            self.centre, self.shape, self.width, self.height, 
            text_width, self.angle, self.args_str, self.s, label)

        return s


class Box2(Shape):
    """Square box,  A rectangle is created by defining aspect."""

    shape = 'rectangle'

    @property
    def coords(self):
        return ((-0.5, 0), (0.5, 0))


class Box4(Shape):
    """Box4"""

    shape = 'rectangle'

    @property
    def coords(self):
        return ((-0.5, 0), (0, -0.5), (0.5, 0), (0, 0.5))


class Box12(Shape):
    """Box12"""

    shape = 'rectangle'

    @property
    def coords(self):
        return ((-0.5, 0.25), (-0.5, 0), (-0.5, -0.25),
                (-0.25, -0.5), (0, -0.5), (0.25, -0.5),
                (0.5, -0.25), (0.5, 0), (0.5, 0.25),
                (0.25, 0.5), (0, 0.5), (-0.25, 0.5))


class Box(Shape):
    """Box"""

    shape = 'rectangle'    
    anchors = {'nw' : (-0.5, 0.5), 'wnw' : (-0.5, 0.25),
               'w' : (-0.5, 0), 'wsw' : (-0.5, -0.25), 
               'sw' : (-0.5, -0.5), 'ssw' : (-0.25, -0.5),
               's' : (0, -0.5), 'sse' : (0.25, -0.5),
               'se' : (0.5, -0.5), 'ese' : (0.5, -0.25),
               'e' : (0.5, 0), 'ene' : (0.5, 0.25),
               'ne' : (0.5, 0.5), 'nne' : (0.25, 0.5),
               'n' : (0, 0.5), 'nnw' : (-0.25, 0.5),
               'mid' : (0.0, 0.0)}
    required_anchors = ('mid', )    


class Ellipse(Shape):
    """Ellipse"""

    # Ellipse needs the tikz shapes library.
    shape = 'ellipse'
    anchors = {'nw' : (-0.3536, 0.3536), 'wnw' : (-0.4619, 0.1913),
               'w' : (-0.5, 0), 'wsw' : (-0.4619, -0.1913), 
               'sw' : (-0.3536, -0.3536), 'ssw' : (-0.1913, -0.4619),
               's' : (0, -0.5), 'sse' : (0.1913, -0.4619),
               'se' : (0.3536, -0.3536), 'ese' : (0.4619, -0.1913),
               'e' : (0.5, 0), 'ene' : (0.4619, 0.1913),
               'ne' : (0.3536, 0.35365), 'nne' : (0.1913, 0.4619),
               'n' : (0, 0.5), 'nnw' : (-0.1913, 0.4619),
               'mid' : (0.0, 0.0)}
    required_anchors = ('mid', )


class Circle(Ellipse):
    """Circle"""

    shape = 'circle'


class Circle2(Shape):
    """Circle"""

    shape = 'circle'

    @property
    def coords(self):
        return ((-0.5, 0), (0.5, 0))


class Circle4(Shape):
    """Circle4"""

    shape = 'circle'

    @property
    def coords(self):
        return ((-0.5, 0), (0, -0.5), (0.5, 0), (0, 0.5))


class Triangle(Shape):
    """Equilateral triangle, The triangle shape can be altered by defining
    defining aspect."""    

    shape = 'triangle'
    # 1 / sqrt(3) approx 0.5774, 1 / (2 * sqrt(3)) approx 0.2887
    anchors = {'c1' : (0.0, 0.5774),
               'c2' : (-0.5, -0.2887),
               'c3' : (0.5, -0.2887),
               'mid' : (0.0, 0.0)}
    required_anchors = ('mid', 'c1', 'c2', 'c3')

    def draw(self, **kwargs):

        if not self.check():
            return ''

        s = self.draw_path([self.node('c1').pos, self.node('c2').pos,
                            self.node('c3').pos], closed=True)
        s += self.draw_label(self.node('mid').pos)

        return s

class TR(Box2):
    """Transfer function"""

    default_width = 1.5
    default_aspect = 1.5


class Chip(Shape):
    """General purpose chip"""

    default_width = 2.0

    # Could allow can_scale but not a lot of point since nodes
    # will not be on the boundary of the chip.

    # TODO, tweak coord if pin name ends in \ using pinpos to
    # accomodate inverting circle.  This will require stripping of the
    # \ from the label. Alternatively, do not use inverting circle and
    # add overline to symbol name when printing.

    @property
    def path(self):
        return ((-0.5, 0.5), (0.5, 0.5), (0.5, -0.5), (-0.5, -0.5))

    def pinmap(self, pos):

        pinmap = ['l', 't', 'r', 'b']
        if pos not in pinmap:
            return pos
            
        index = pinmap.index(pos)
        angle = int(self.angle)
        if angle < 0:
            angle += 360

        angles = (0, 90, 180, 270)
        if angle not in angles:
            raise ValueError('Cannot rotate pinpos %s by %s' % (pos, self.angle))

        index += angles.index(angle)
        pos = pinmap[index % len(pinmap)]
        return pos

    def name_pins(self):

        pins = self.opts.get('pins', '')
        if pins != '' and pins != 'auto':
            if pins[0] != '{':
                raise ValueError('Expecting { for pins in %s' % self)
            if pins[-1] != '}':
                raise ValueError('Expecting } for pins in %s' % self)
            pins = pins[1:-1]
                
            pins = pins.split(',')
            if len(pins) != len(self.nodes):
                raise ValueError('Expecting %d pin names, got %s in %s' % (
                    len(self.nodes), len(pins), self))

        for m, n in enumerate(self.nodes):
            n.pinpos = self.pinmap(self.pinpos[m])
            label = ''
            if pins == 'auto':
                label = n.name.split('.')[-1]
                if label[0] == '_':
                    label = ''
            elif pins != '':
                label = pins[m].strip()

            n.clock = label != '' and label[0] == '>'
            if n.clock:
                # Remove clock designator
                label = label[1:]

            n.label = label

    def draw(self, **kwargs):

        if not self.check():
            return ''

        self.name_pins()
            
        centre = self.centre
        q = self.tf(centre, self.path)

        s = self.draw_path(q, closed=True, style='thick')
        s += r'  \draw (%s) node[text width=%scm, align=center, %s] {%s};''\n'% (
            centre, self.width - 0.5, self.args_str, self.label(**kwargs))

        for m, n in enumerate(self.nodes):
            if n.clock:
                # TODO, tweak for pinpos
                q = self.tf(n.pos, ((0, 0.125 * 0.707), (0.125, 0), 
                                    (0, -0.125 * 0.707)))
                s += self.draw_path(q[0:3], style='thick')
                
        return s


class Uchip1310(Chip):
    """Chip of size 1 3 1 0"""

    default_aspect = 4.0 / 3.0
    pinpos = ('l', 'b', 'b', 'b', 'r')

    @property
    def centre(self):
        return self.midpoint(self.nodes[0], self.nodes[4])

    @property
    def coords(self):
        return ((0, 0), (0.25, -0.5), (0.5, -0.5), (0.75, -0.5), (1, 0))


class Uchip2121(Chip):
    """Chip of size 2 1 2 1"""

    pinpos = ('l', 'l', 'b', 'r', 'r', 't')

    @property
    def coords(self):
        return ((0, 0.25), (0, -0.25),
                (0.5, -0.5), 
                (1.0, -0.25), (1.0, 0.25),
                (0.5, 0.5))


class Uchip3131(Chip):
    """Chip of size 3 1 3 1"""

    pinpos = ('l', 'l', 'l', 'b', 'r', 'r', 'r', 't')

    @property
    def coords(self):
        return ((-0.5, 0.25), (-0.5, 0), (-0.5, -0.25), (0.0, -0.5), 
                (0.5, -0.25), (0.5, 0), (0.5, 0.25), (0.0, 0.5))

    @property
    def path(self):
        return ((-0.5, 0.375), (0.5, 0.375), (0.5, -0.375), (-0.5, -0.375))


class Uchip4141(Chip):
    """Chip of size 4 1 4 1"""

    default_width = 2
    default_aspect = 0.5
    pinpos = ('l', 'l', 'l', 'l', 'b', 'r', 'r', 'r', 'r', 't')

    @property
    def coords(self):
        return ((-0.5, 0.375), (-0.5, 0.125), (-0.5, -0.125), (-0.5, -0.375),
                (0.0, -0.5), 
                (0.5, -0.375), (0.5, -0.125), (0.5, 0.125), (0.5, 0.375),
                (0, 0.5))


class Uadc(Chip):
    """ADC"""

    # in, vref, vss, clk, data, fs, vdd, vref
    pinpos = ('l', 'b', 'b', 'r', 'r', 'r', 't', 't')

    @property
    def coords(self):
        return ((0, 0.0), (0.5, -0.5), (0.75, -0.5), 
                (1.0, -0.25), (1.0, 0), (1.0, 0.25), (0.75, 0.5), (0.5, 0.5))

    @property
    def path(self):
        return ((-0.5, 0.0), (-0.25, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.25, 0.5))

    @property
    def centre(self):
        return self.midpoint(self.nodes[0], self.nodes[4])


class Udac(Chip):
    """DAC"""

    # fs, data, clk, vss, out, vdd, vref
    pinpos = ('l', 'l', 'l', 'b', 'b', 'r', 't', 't')

    @property
    def coords(self):
        return ((0, 0.25), (0, 0), (0, -0.25), (0.25, -0.5), (0.5, -0.5),
                (1.0, 0), (0.5, 0.5), (0.25, 0.5))

    @property
    def path(self):
        return ((-0.5, -0.5), (0.25, -0.5), (0.5, 0), (0.25, 0.5), (-0.5, 0.5))

    @property
    def centre(self):
        return self.midpoint(self.nodes[1], self.nodes[5])


class Udiffamp(Chip):
    """Amplifier"""

    default_width = 1.0
    pinpos = ('l', 'l', 'b', 'r', 't')

    @property
    def coords(self):
        return ((-0.5, 0.25), (-0.5, -0.25), (0.0, -0.25), (0.5, 0), (0.0, 0.25))

    @property
    def path(self):
        return ((-0.5, 0.5), (-0.5, -0.5), (0.5, 0))

    @property
    def centre(self):
        n1, n2, n3, n4, n5 = self.nodes
        return (n1.pos + n2.pos) * 0.25 + n4.pos * 0.5


class Ubuffer(Chip):
    """Buffer with power supplies"""

    default_width = 1.0
    pinpos = ('l', 'b', 'r', 't')

    @property
    def coords(self):
        return ((0, 0), (0.5, -0.25), (1.0, 0), (0.5, 0.25))

    @property
    def path(self):
        return ((-0.5, 0.5), (0.5, 0), (-0.5, -0.5))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        self.name_pins()

        n1, n2, n3, n4 = self.nodes
        centre = (n1.pos + n3.pos) * 0.5

        q = self.tf(centre, self.path)
        s = self.draw_path(q[0:3], closed=True, style='thick')
        s += self.draw_label(centre, **kwargs)
        s += self.draw_nodes(**kwargs)
        return s


class Uinverter(Chip):
    """Inverter with power supplies"""

    default_width = 1.0
    pinpos = ('l', 'b', 'r', 't')

    @property
    def coords(self):
        return ((0, 0), (0.5, -0.22), (1.0, 0), (0.5, 0.22))

    @property
    def path(self):
        w = 0.05
        return ((-0.5, 0.5), (0.5 -2 * w, 0), (-0.5, -0.5), (0.5 - w, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        self.name_pins()

        n1, n2, n3, n4 = self.nodes
        centre = (n1.pos + n3.pos) * 0.5

        q = self.tf(centre, self.path)
        s = self.draw_path(q[0:3], closed=True, style='thick')
        s += r'  \draw[thick] (%s) node[ocirc, scale=%s] {};''\n' % (
            q[3], 1.8 * self.size * self.scale)
        s += self.draw_label(centre, **kwargs)
        s += self.draw_nodes(**kwargs)
        return s


class Wire(OnePort):

    def __init__(self, sch, name, cpt_type, cpt_id, string,
                 opts_string, node_names, keyword, *args):

        implicit = False
        for key in self.implicit_keys:
            if key in opts_string:
                implicit = True
                break

        if implicit:
            # Rename second node since this is spatially different from
            # other nodes of the same name.  Add underscore so node
            # not drawn.
            node_names = (node_names[0], name + '@_' + node_names[1])
        
        super (Wire, self).__init__(sch, name, cpt_type, cpt_id, string,
                                    opts_string, node_names, keyword, *args)
        self.implicit = implicit

    @property
    def coords(self):
        return ((0, 0), (1, 0))

    def draw_implicit(self, **kwargs):
        """Draw implicit wires, i.e., connections to ground, etc."""

        kind = None
        for key in self.implicit_keys:
            if key in self.opts:
                kind = key
                break;

        # I like the sground symbol for power supplies but rground symbol
        # is also common.
        if (kind is None) or (kind == 'implicit'):
            kind = 'sground'
        anchor = 'south west'
        if self.down:
            anchor = 'north west'

        n1, n2 = self.nodes
        s = self.draw_path((n1.s, n2.s))
        s += r'  \draw (%s) node[%s, scale=0.5, rotate=%d] {};''\n' % (
            n2.s, kind, self.angle + 90)

        if 'l' in self.opts:
            lpos = self.tf(n2.pos, (0.125, 0))
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

        n1, n2 = self.nodes

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
        s += self.draw_nodes(**kwargs)

        if self.voltage_str != '':
            print('There is no voltage drop across an ideal wire!')

        if self.current_str != '' or self.label_str != '':
            # FIXME, we don't really want the wire drawn since this
            # can clobber the arrow.  We just want the current
            # annotation and/or the label.

            # To handle multiple labels, we need to draw separate wires.
            for label_str in self.label_str_list:
                s += r'  \draw[%s] (%s) [short, %s, %s] to (%s);''\n' % (
                    self.args_str, n1.s, self.current_str, label_str, n2.s)
            if self.label_str_list == []:
                s += r'  \draw[%s] (%s) [short, %s] to (%s);''\n' % (
                    self.args_str, n1.s, self.current_str, n2.s)
        return s


class FB(StretchyCpt):
    """Ferrite bead"""

    can_scale = True

    @property
    def coords(self):
        return ((-0.5, 0), (0.5, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.nodes

        centre = (n1.pos + n2.pos) * 0.5
        w = 0.125
        h = 0.4
        
        q = self.tf(centre, ((-0.5 * w, -0.5 * h), (-0.5 * w, 0.5 * h),
                             (0.5 * w, 0.5 * h), (0.5 * w, -0.5 * h),
                             (0, h)), -30)
        q2 = self.tf(centre, ((-0.53 * w, 0), (0.53 * w, 0), (0, h)))

        s = self.draw_path(q[0:4], closed=True, style='thick')
        s += self.draw_path((n1.s, q2[0]))
        s += self.draw_path((q2[1], n2.s))
        s += self.draw_label(q[4], **kwargs)
        s += self.draw_nodes(**kwargs)
        return s


class XT(StretchyCpt):
    """Crystal"""

    can_scale = True

    @property
    def coords(self):
        return ((-0.5, 0), (0.5, 0))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.nodes

        centre = (n1.pos + n2.pos) * 0.5
        q = self.tf(centre, ((-0.15, 0), (-0.15, 0.15), (-0.15, -0.15),
                             (0.15, 0), (0.15, 0.15), (0.15, -0.15),
                             (-0.06, 0.15), (0.06, 0.15),
                             (0.06, -0.15), (-0.06, -0.15),
                             (0.0, -0.3)))

        s = self.draw_path((q[1], q[2]), style='thick')
        s += self.draw_path((q[4], q[5]), style='thick')
        s += self.draw_path(q[6:10], closed=True, style='thick')
        s += self.draw_path((q[0], n1.s), style='thick')
        s += self.draw_path((q[3], n2.s), style='thick')
        s += self.draw_label(q[10], **kwargs)
        s += self.draw_nodes(**kwargs)
        return s


classes = {}

def defcpt(name, base, docstring, cpt=None):
    
    if isinstance(base, str):
        base = classes[base]

    newclass = type(name, (base, ), {'__doc__': docstring})

    if cpt is not None:
        newclass.tikz_cpt = cpt
    classes[name] = newclass


def make(classname, parent, name, cpt_type, cpt_id,
         string, opts_string, node_names, *args):

    # Create instance of component object
    try:
        newclass = getattr(module, classname)
    except:
        newclass = classes[classname]

    cpt = newclass(parent, name, cpt_type, cpt_id, string, opts_string, 
                   node_names, *args)
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
defcpt('Istep', 'I', 'Step current source', 'I')
defcpt('Iac', 'I', 'AC current source', 'sI')
defcpt('Inoise', 'I', 'Noise current source', 'sI')

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

defcpt('Sbox', Box, 'Box shape')
defcpt('Scircle', Circle, 'Circle shape')
defcpt('Sellipse', Ellipse, 'Ellipse shape')
defcpt('Striangle', Triangle, 'Triangle shape')

defcpt('SW', OnePort, 'Switch', 'closing switch')
defcpt('SWno', 'SW', 'Normally open switch', 'closing switch')
defcpt('SWnc', 'SW', 'Normally closed switch', 'opening switch')
defcpt('SWpush', 'SW', 'Pushbutton switch', 'push button')
defcpt('SWspdt', SPDT, 'SPDT switch', 'spdt')

defcpt('TF', Transformer, 'Transformer', 'transformer')
defcpt('TFcore', Transformer, 'Transformer with core', 'transformer core')
defcpt('TFtapcore', TFtap, 'Tapped transformer with core', 'transformer core')
defcpt('TP', TwoPort, 'Two port', '')

defcpt('Ubox', Box2, 'Box')
defcpt('Ucircle', Circle2, 'Circle')
defcpt('Ubox4', Box4, 'Box')
defcpt('Ubox12', Box12, 'Box')
defcpt('Ucircle4', Circle4, 'Circle')


defcpt('V', OnePort, 'Voltage source', 'V')
defcpt('sV', OnePort, 'Voltage source', 'V')
defcpt('Vsin', 'V', 'Sinusoidal voltage source', 'sV')
defcpt('Vdc', 'V', 'DC voltage source', 'V')
defcpt('Vstep', 'V', 'Step voltage source', 'V')
defcpt('Vac', 'V', 'AC voltage source', 'sV')
defcpt('Vnoise', 'V', 'Noise voltage source', 'sV')

defcpt('VM', OnePort, 'Voltmeter', 'voltmeter')

defcpt('W', Wire, 'Wire', 'short')
defcpt('Y', OnePort, 'Admittance', 'european resistor')
defcpt('Z', OnePort, 'Impedance', 'european resistor')

# Perhaps AM for ammeter, VM for voltmeter, VR for variable resistor?
# Currently, a variable resistor is supported with the variable
# option.
