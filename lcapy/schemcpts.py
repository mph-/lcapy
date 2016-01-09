"""
This module defines and draws the schematic components using
circuitikz.   The components are defined at the bottom of this file.

Copyright 2015, 2016 Michael Hayes, UCECE
"""


from __future__ import print_function
from lcapy.latex import latex_str
from lcapy.schemmisc import Pos
import numpy as np
import sys

module = sys.modules[__name__]


class Cpt(object):

    def __init__(self, sch, cpt_type, cpt_id, string, opts_string, nodes, *args):

        self.sch = sch
        self.type = cpt_type
        self.id = cpt_id

        if cpt_id == '' and sch is not None:
            if cpt_type not in sch.anon:
                sch.anon[cpt_type] = 0
            sch.anon[cpt_type] += 1
            cpt_id = '#%d' % sch.anon[cpt_type]

        name = self.type + cpt_id

        self.string = string
        self.opts_string = opts_string
        # There are three sets of nodes:
        # 1. nodes are the names of the electrical nodes for a cpt.
        # 2. vnodes are the subset of the electrical nodes for a cpt that are drawn.
        # 3. dnodes includes vnodes plus inherited nodes from other cpts (K).  The
        # lookup of this attribute is deferred until cpts are drawn.
        self.nodes = nodes
        self.name = name
        self.args = args
        self.classname = self.__class__.__name__

    def __repr__(self):

        if hasattr(self, 'string'):
            return self.string
        
        return type(self)

    @property
    def size(self):
        return float(self.opts['size'])

    @property
    def down(self):
        return self.opts['dir']  == 'down'

    @property
    def up(self):
        return self.opts['dir']  == 'up'

    @property
    def left(self):
        return self.opts['dir']  == 'left'

    @property
    def right(self):
        return self.opts['dir']  == 'right'

    @property
    def horizontal(self):
        return self.opts['dir'] in ('left', 'right')

    @property
    def vertical(self):
        return self.opts['dir'] in ('left', 'down')

    @property
    def mirror(self):
        return 'mirror' in self.opts

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
            raise ValueError('Unknown orientation: %s' % self.opts['dir'])
        
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
    def variable(self):
        return 'variable' in self.opts

    @property
    def vnodes(self):
        '''Visible nodes'''
        return self.nodes

    @property
    def dnodes(self):
        '''Nodes used to construct schematic'''
        return self.vnodes

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
        for m1, n1 in enumerate(self.dnodes):
            for m2, n2 in enumerate(self.dnodes[m1 + 1:], m1 + 1):
                if xvals[m2] == xvals[m1]:
                    graphs.link(n1, n2)

    def ylink(self, graphs):

        yvals = self.yvals
        for m1, n1 in enumerate(self.dnodes):
            for m2, n2 in enumerate(self.dnodes[m1 + 1:], m1 + 1):
                if yvals[m2] == yvals[m1]:
                    graphs.link(n1, n2)

    def xplace(self, graphs):

        size = self.size
        xvals = self.xvals
        for m1, n1 in enumerate(self.dnodes):
            for m2, n2 in enumerate(self.dnodes[m1 + 1:], m1 + 1):
                value = (xvals[m2] - xvals[m1]) * size
                graphs.add(n1, n2, value)

    def yplace(self, graphs):

        size = self.size
        yvals = self.yvals
        for m1, n1 in enumerate(self.dnodes):
            for m2, n2 in enumerate(self.dnodes[m1 + 1:], m1 + 1):
                value = (yvals[m2] - yvals[m1]) * size
                graphs.add(n1, n2, value)

    def _node_str(self, node1, node2, draw_nodes=True):

        node_str = ''
        if node1.visible(draw_nodes):
            node_str = 'o' if node1.port else '*'

        node_str += '-'

        if node2.visible(draw_nodes):
            node_str += 'o' if node2.port else '*'

        if node_str == '-':
            node_str = ''
        
        return node_str

    def _draw_node(self, n, draw_nodes=True):
        
        s = ''
        if not draw_nodes:
            return s

        node = self.sch.nodes[n]
        if not node.visible(draw_nodes):
            return s

        if node.port:
            s = r'  \draw (%s) node[ocirc] {};''\n' % node.pos
        else:
            s = r'  \draw (%s) node[circ] {};''\n' % node.pos

        return s

    def _draw_nodes(self, draw_nodes=True):

        s = ''
        for n in self.dnodes:
            s += self._draw_node(n, draw_nodes)
        return s

    def draw(self, **kwargs):
        raise NotImplementedError('draw method not implemented for %s' % self)


class Transistor(Cpt):
    """Transistor"""
    
    npos = ((1, 1.5), (0, 0.75), (1, 0))
    ppos = ((1, 0), (0, 0.75), (1, 1.5))

    @property
    def coords(self):
        if self.classname in ('Qpnp', 'Mpmos', 'Jpjf'):
            return self.npos if self.mirror else self.ppos
        else:
            return self.ppos if self.mirror else self.npos

    def draw(self, **kwargs):

        label_values = kwargs.get('label_values', True)
        draw_nodes = kwargs.get('draw_nodes', True)

        p1, p2, p3 = [self.sch.nodes[n].pos for n in self.dnodes]
        centre = (p1 + p3) * 0.5

        label_str = '$%s$' % self.default_label if label_values else ''
        args_str = ''
        for key, val in self.opts.iteritems():
            if key in ('color', ):
                args_str += '%s=%s, ' % (key, val)                

        s = r'  \draw (%s) node[%s, %s, scale=%.1f, rotate=%d] (T) {};''\n' % (
            centre, self.tikz_cpt, args_str, self.sch.scale * 2, self.angle)
        s += r'  \draw (%s) node [] {%s};''\n'% (centre, label_str)

        # Add additional wires.
        if self.tikz_cpt in ('pnp', 'npn'):
            s += r'  \draw (T.C) -- (%s) (T.B) -- (%s) (T.E) -- (%s);''\n' % self.dnodes
        else:
            s += r'  \draw (T.D) -- (%s) (T.G) -- (%s) (T.S) -- (%s);''\n' % self.dnodes

        s += self._draw_nodes(draw_nodes)
        return s


class JFET(Transistor):
    """Transistor"""
    
    npos = ((1, 1.5), (0, 0.48), (1, 0))
    ppos = ((1, 0), (0, 1.02), (1, 1.5))


class MOSFET(Transistor):
    """Transistor"""
    
    npos = ((0.85, 1.5), (-0.25, 0.75), (0.85, 0))
    ppos = ((0.85, 0), (-0.25, 0.75), (0.85, 1.5))


class TwoPort(Cpt):
    """Two-port"""

    @property
    def coords(self):
        return ((1, 1), (1, 0), (0, 1), (0, 0))

    def draw(self, **kwargs):

        label_values = kwargs.get('label_values', True)
        draw_nodes = kwargs.get('draw_nodes', True)


        # TODO, fix positions if component rotated.

        p1, p2, p3, p4 = [self.sch.nodes[n].pos for n in self.dnodes]
        width = p2.x - p4.x
        extra = 0.25
        p1.y += extra
        p2.y -= extra
        p3.y += extra
        p4.y -= extra
        centre = (p1 + p2 + p3 + p4) * 0.25
        top = Pos(centre.x, p1.y + 0.15)

        label_str = '$%s$' % self.default_label if label_values else ''
        titlestr = "%s-parameter two-port" % self.args[2]

        s = r'  \draw (%s) -- (%s) -- (%s) -- (%s) -- (%s);''\n' % (
            p4, p3, p1, p2, p4)
        s += r'  \draw (%s) node[minimum width=%.1f] {%s};''\n' % (
            centre, width, titlestr)
        s += r'  \draw (%s) node[minimum width=%.1f] {%s};''\n' % (
            top, width, label_str)

        s += self._draw_nodes(draw_nodes)
        return s


class TF1(TwoPort):
    """Transformer"""

    def draw(self, **kwargs):

        label_values = kwargs.get('label_values', True)
        link = kwargs.get('link', True)

        p1, p2, p3, p4 = [self.sch.nodes[n].pos for n in self.dnodes]

        xoffset = 0.06
        yoffset = 0.40

        # TODO, fix positions if component rotated.
        primary_dot = Pos(p3.x - xoffset, 0.5 * (p3.y + p4.y) + yoffset)
        secondary_dot = Pos(p1.x + xoffset, 0.5 * (p1.y + p2.y) + yoffset)

        centre = (p1 + p2 + p3 + p4) * 0.25
        labelpos = Pos(centre.x, primary_dot.y)

        label_str = '$%s$' % self.default_label if label_values else ''

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


class TF(TF1):
    """Transformer"""

    def draw(self, **kwargs):

        n1, n2, n3, n4 = self.dnodes

        s = r'  \draw (%s) to [inductor] (%s);''\n' % (n3, n4)
        s += r'  \draw (%s) to [inductor] (%s);''\n' % (n1, n2)

        s += super(TF, self).draw(link=True, **kwargs)

        draw_nodes = kwargs.get('draw_nodes', True)
        s += self._draw_nodes(draw_nodes)
        return s


class K(TF1):
    """Mutual coupling"""

    def __init__(self, sch, cpt_type, cpt_id, string, opts_string, nodes, *args):

        self.Lname1 = args[0]
        self.Lname2 = args[1]
        super (K, self).__init__(sch, cpt_type, cpt_id, string, opts_string, nodes, *args)

    @property
    def dnodes(self):

        # L1 and L2 need to be previously defined so we can find their nodes.
        L1 = self.sch.elements[self.Lname1]
        L2 = self.sch.elements[self.Lname2]
        return L1.nodes + L2.nodes


class OnePort(Cpt):
    """OnePort"""

    @property
    def coords(self):
        return ((0, 0), (1, 0))

    def draw(self, **kwargs):

        label_values = kwargs.get('label_values', True)
        draw_nodes = kwargs.get('draw_nodes', True)
        label_ids = kwargs.get('label_ids', True)

        n1, n2 = self.dnodes

        tikz_cpt = self.tikz_cpt
        if self.type == 'R' and self.variable:
            tikz_cpt = 'vR'

        id_pos = '_'
        voltage_pos = '^'
        if self.type in ('V', 'I', 'E', 'F', 'G', 'H'):

            # circuitikz expects the positive node first, except for
            # voltage and current sources!  So swap the nodes
            # otherwise they are drawn the wrong way around.
            n1, n2 = n2, n1

            if self.horizontal:
                # Draw label on LHS for vertical cpt and below
                # for horizontal cpt.
                id_pos = '^'
                voltage_pos = '_'
        else:
            if self.vertical:
                # Draw label on LHS for vertical cpt and below
                # for horizontal cpt.
                id_pos = '^'
                voltage_pos = '_'

        # Add modifier to place voltage label on other side
        # from component identifier label.
        if 'v' in self.opts:
            self.opts['v' + voltage_pos] = self.opts.pop('v')

        # Reversed voltage.
        if 'vr' in self.opts:
            self.opts['v' + voltage_pos + '>'] = self.opts.pop('vr')

        current_pos = id_pos
        # Add modifier to place current label on other side
        # from voltage marks.
        if 'i' in self.opts:
            self.opts['i' + current_pos] = self.opts.pop('i')

        # Reversed current.
        if 'ir' in self.opts:
            self.opts['i' + current_pos + '<'] = self.opts.pop('ir')

        # Current, voltage, label options.
        # It might be better to allow any options and prune out
        # dir and size.
        voltage_str = ''
        current_str = ''
        label_str = ''
        args_str = ''
        for key, val in self.opts.iteritems():
            if key in ('i', 'i_', 'i^', 'i_>', 'i_<', 'i^>', 'i^<',
                       'i>_', 'i<_', 'i>^', 'i<^'):
                current_str += '%s=${%s}$, ' % (key, latex_str(val))
            elif key in ('v', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<'):
                voltage_str += '%s=${%s}$, ' % (key, latex_str(val))
            elif key in ('l', 'l^', 'l_'):
                label_str += '%s=${%s}$, ' % (key, latex_str(val))
            elif key in ('color', ):
                args_str += '%s=%s, ' % (key, val)                

        node_str = self._node_str(self.sch.nodes[n1], self.sch.nodes[n2],
                                  draw_nodes)

        args_str += voltage_str + current_str

        # Generate default label unless specified.
        if label_str == '':
            if self.type not in ('O', 'W'):
                
                label_str = ', l%s=${%s}$' % (id_pos, self.default_label)
                
                if label_ids and self.value_label != '':
                    label_str = r', l%s=${%s=%s}$' % (
                        id_pos, self.id_label, self.value_label)
        else:
            label_str = ', ' + label_str

        if not label_values:
            label_str = ''

        if not label_values and label_ids:
            label_str = ', l%s=${%s}$' % (id_pos, self.id_label)

        s = r'  \draw (%s) to [align=right, %s%s, %s%s] (%s);''\n' % (
            n1, tikz_cpt, label_str, args_str, node_str, n2)
        return s


class VCS(OnePort):
    """Voltage controlled source"""

    @property
    def vnodes(self):
        return self.nodes[0:2]


class CCS(OnePort):
    """Current controlled source"""


class Opamp(Cpt):

    # The Nm node is not used (ground).
    ppos = ((2.4, 0.5), (0, 1), (0, 0))
    npos = ((2.4, 0.5), (0, 0), (0, 1))

    @property
    def vnodes(self):
        return (self.nodes[0], ) + self.nodes[2:]

    @property
    def coords(self):
        return self.npos if self.mirror else self.ppos

    def draw(self, **kwargs):

        draw_nodes = kwargs.get('draw_nodes', True)
        label_values = kwargs.get('label_values', True)

        p1, p3, p4 = [self.sch.nodes[n].pos for n in self.dnodes]

        centre = (p3 + p4) * 0.25 + p1 * 0.5

        label_str = '$%s$' % self.default_label if label_values else ''
        args_str = ''
        for key, val in self.opts.iteritems():
            if key in ('color', ):
                args_str += '%s=%s, ' % (key, val)                

        s = r'  \draw (%s) node[op amp, %s, scale=%.3f, rotate=%d] (opamp) {};' % (
            centre, args_str, self.sch.scale * 2 * 1.01, self.angle)
        # Draw label separately to avoid being scaled by 2.
        s += r'  \draw (%s) node [] {%s};' % (centre, label_str)
        
        s += self._draw_nodes(draw_nodes)
        return s


class FDOpamp(Cpt):

    @property
    def coords(self):
        return ((2, 1), (2, 0), (0, 1), (0, 0))

    def draw(self, **kwargs):

        draw_nodes = kwargs.get('draw_nodes', True)
        label_values = kwargs.get('label_values', True)

        p1, p2, p3, p4 = [self.sch.nodes[n].pos for n in self.dnodes]

        centre = (p1 + p2 + p3 + p4) * 0.25 + np.dot((0.18, 0), self.R)

        label_str = '$%s$' % self.default_label if label_values else ''
        args_str = ''
        for key, val in self.opts.iteritems():
            if key in ('color', ):
                args_str += '%s=%s, ' % (key, val)                

        s = r'  \draw (%s) node[fd op amp, %s, scale=%.3f, rotate=%d] (opamp) {};' % (
            centre, args_str, self.sch.scale * 2 * 1.015, self.angle)
        # Draw label separately to avoid being scaled by 2.
        s += r'  \draw (%s) node [] {%s};' % (centre, label_str)
        
        s += self._draw_nodes(draw_nodes)
        return s


class Wire(OnePort):

    @property
    def coords(self):

        if 'implicit' in self.opts or 'ground' in self.opts or 'sground' in self.opts:
            return ((0, 0), )
        return ((0, 0), (1, 0))

    @property
    def vnodes(self):

        if 'implicit' in self.opts or 'ground' in self.opts or 'sground' in self.opts:
            return (self.nodes[0], )
        return self.nodes

    def draw_implicit(self, **kwargs):
        """Draw implict wires, i.e., connections to ground, etc."""

        kind = ''
        if 'implicit' in self.opts:
            kind = 'sground'
        elif 'ground' in self.opts:
            kind = 'ground'
        elif 'sground' in self.opts:
            kind = 'sground'

        anchor = 'south west'
        if self.down:
            anchor = 'north west'

        n1 = self.dnodes[0]
        s = r'  \draw (%s) node[%s, rotate=%d] {};''\n' % (
            n1, kind, self.angle + 90)

        lpos = self.sch.nodes[n1].pos + np.dot((0.25, 0), self.R)

        if 'l' in self.opts:
            label_str = '${%s}$' % latex_str(self.opts['l'])
            s += r'  \draw {[anchor=%s] (%s) node {%s}};''\n' % (
                anchor, lpos, label_str)
        return s

    def draw(self, **kwargs):

        if 'implicit' in self.opts or 'ground' in self.opts or 'sground' in self.opts:
            return self.draw_implicit(**kwargs)
                                    
        return super(Wire, self).draw(**kwargs)


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

defcpt('SW', OnePort, 'SWitch', 'closing switch')
defcpt('SWno', 'SW', 'Normally open switch', 'closing switch')
defcpt('SWnc', 'SW', 'Normally closed switch', 'opening switch')
defcpt('SWpush', 'SW', 'Pushbutton switch', 'push button')

defcpt('TF', TwoPort, 'Transformer', 'transformer')
defcpt('TP', TwoPort, 'Two port', '')

defcpt('V', OnePort, 'Voltage source', 'V')
defcpt('sV', OnePort, 'Voltage source', 'V')
defcpt('Vsin', 'V', 'Sinusoidal voltage source', 'sV')
defcpt('Vdc', 'V', 'DC voltage source', 'V')
defcpt('Vstep', 'V', 'Step voltage source', 'V')
defcpt('Vac', 'V', 'AC voltage source', 'sV')

defcpt('W', Wire, 'Wire', 'short')
defcpt('Y', OnePort, 'Admittance', 'european resistor')
defcpt('Z', OnePort, 'Impedance', 'european resistor')

# Perhaps AM for ammeter, VM for voltmeter, VR for variable resistor?
# Currently, a variable resistor is supported with the variable
# option.
