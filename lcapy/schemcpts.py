from __future__ import print_function
from lcapy.latex import latex_str
from lcapy.schemmisc import Pos
import numpy as np


class Cpt(object):

    yscale = 1.0
    xscale = 1.0
    cpt_type_counter = 1

    pos = ((0, 0), (0, 1))


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
    def variable(self):
        return 'variable' in self.opts

    @property
    def vnodes(self):
        '''Drawn nodes'''
        return self.nodes

    @property
    def tpos(self):
        """Transformed node positions"""

        if hasattr(self, '_tpos'):
            return self._tpos

        # right + -
        #
        # down  +
        #       -
        #
        # left  - +
        #
        # up    -
        #       +

        tpos = np.array(self.pos)
        if self.left:            
            # Negate x.
            tpos = np.dot(tpos, np.array(((-1, 0), (0, 1))))
        elif self.up: 
            # Swap x/y. 
            tpos = np.dot(tpos, np.array(((0, 1), (1, 0))))
        elif self.down: 
            # Swap x/y and negate y. 
            tpos = np.dot(tpos, np.array(((0, -1), (1, 0))))
        elif not self.right: 
            raise ValueError('Unknown orientation: %s' % self.opts['dir'])
        self._tpos = tpos
        return tpos

    @property
    def xvals(self):
        return self.tpos[:, 0]

    @property
    def yvals(self):
        return self.tpos[:, 1]

    def __init__(self, name, *args, **kwargs):

        self.name = name
        self.args = args

    def __repr__(self):

        if hasattr(self, 'string'):
            return self.string
        
        return type(self)

    def fixup(self, sch):
        pass

    def xlink(self, graphs):

        xvals = self.xvals
        for m1, n1 in enumerate(self.vnodes):
            for m2, n2 in enumerate(self.vnodes[m1 + 1:], m1 + 1):
                if xvals[m2] == xvals[m1]:
                    graphs.link(n1, n2)

    def ylink(self, graphs):

        yvals = self.yvals
        for m1, n1 in enumerate(self.vnodes):
            for m2, n2 in enumerate(self.vnodes[m1 + 1:], m1 + 1):
                if yvals[m2] == yvals[m1]:
                    graphs.link(n1, n2)

    def xplace(self, graphs):

        size = self.size
        xvals = self.xvals
        for m1, n1 in enumerate(self.vnodes):
            for m2, n2 in enumerate(self.vnodes[m1 + 1:], m1 + 1):
                value = (xvals[m2] - xvals[m1]) * self.xscale * size
                graphs.add(n1, n2, value)

    def yplace(self, graphs):

        size = self.size
        yvals = self.yvals
        for m1, n1 in enumerate(self.vnodes):
            for m2, n2 in enumerate(self.vnodes[m1 + 1:], m1 + 1):
                value = (yvals[m2] - yvals[m1]) * self.yscale * size
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

    def _draw_node(self, sch, n, draw_nodes=True):
        
        s = ''
        if not draw_nodes:
            return s

        node = sch.nodes[n]
        if not node.visible(draw_nodes):
            return s

        if node.port:
            s = r'  \draw (%s) node[ocirc] {};''\n' % node.pos
        else:
            s = r'  \draw (%s) node[circ] {};''\n' % node.pos

        return s

    def _draw_nodes(self, sch, draw_nodes=True):

        s = ''
        for n in self.vnodes:
            s += self._draw_node(sch, n, draw_nodes)
        return s

    def draw(self, sch, **kwargs):
        raise NotImplementedError('draw method not implemented')


class Transistor(Cpt):
    """Transistor"""
    
    yscale = 1.5
    xscale = 0.85
    npos = ((1, 1), (0, 0.5), (1, 0))
    ppos = ((1, 0), (0, 0.5), (1, 1))

    def fixup(self, sch):

        if self.classname in ('Qpnp', 'Mpmos', 'Jpjf'):
            self.pos = self.npos if self.mirror else self.ppos
        else:
            self.pos = self.ppos if self.mirror else self.npos

    def draw(self, sch, **kwargs):

        label_values = kwargs.get('label_values', True)
        draw_nodes = kwargs.get('draw_nodes', True)

        n1, n2, n3 = self.vnodes
        p1, p2, p3 = sch.nodes[n1].pos, sch.nodes[n2].pos, sch.nodes[n3].pos

        centre = (p1 + p3) * 0.5

        label_str = '$%s$' % self.default_label if label_values else ''
        args_str = ''
        if self.mirror:
            if self.vertical:
                args_str += ', xscale=-1'
            else:
                args_str += ', yscale=-1'
        for key, val in self.opts.iteritems():
            if key in ('color', ):
                args_str += '%s=%s, ' % (key, val)                

        angle = 0
        if self.down:
            angle = -90
        if self.up:
            angle = 90

        s = r'  \draw (%s) node[%s, %s, scale=%.1f, rotate=%d] (T) {};''\n' % (
            centre, self.tikz_cpt, args_str, sch.scale * 2, angle)
        s += r'  \draw (%s) node [] {%s};''\n'% (centre, label_str)

        # Add additional wires.
        if self.tikz_cpt in ('pnp', 'npn'):
            s += r'  \draw (T.C) -- (%s) (T.B) -- (%s) (T.E) -- (%s);''\n' % (n1, n2, n3)
        else:
            s += r'  \draw (T.D) -- (%s) (T.G) -- (%s) (T.S) -- (%s);''\n' % (n1, n2, n3)

        s += self._draw_nodes(sch, draw_nodes)
        return s


class TwoPort(Cpt):
    """Two-port"""

    pos = ((1, 1), (1, 0), (0, 1), (0, 0))


    def draw(self, sch, **kwargs):

        label_values = kwargs.get('label_values', True)
        draw_nodes = kwargs.get('draw_nodes', True)

        if not self.right:
            raise ValueError('Cannot draw twoport network %s in direction %s'
                             % (self.name, self.opts['dir']))

        p1, p2, p3, p4 = [sch.nodes[n].pos for n in self.nodes]
        width = p2.x - p4.x
        # height = p1.y - p2.y
        extra = 0.25
        p1.y += extra
        p2.y -= extra
        p3.y += extra
        p4.y -= extra
        centre = Pos(0.5 * (p3.x + p1.x), 0.5 * (p2.y + p1.y))
        top = Pos(centre.x, p1.y + 0.15)

        label_str = '$%s$' % self.default_label if label_values else ''
        titlestr = "%s-parameter two-port" % self.args[2]

        s = r'  \draw (%s) -- (%s) -- (%s) -- (%s) -- (%s);''\n' % (
            p4, p3, p1, p2, p4)
        s += r'  \draw (%s) node[minimum width=%.1f] {%s};''\n' % (
            centre, width, titlestr)
        s += r'  \draw (%s) node[minimum width=%.1f] {%s};''\n' % (
            top, width, label_str)

        s += self._draw_nodes(sch, draw_nodes)
        return s


class JFET(Transistor):
    """Transistor"""
    
    npos = ((1, 1), (0, 0.32), (1, 0))
    ppos = ((1, 0), (0, 0.68), (1, 1))


class TF1(TwoPort):
    """Transformer"""

    def draw(self, sch, **kwargs):

        label_values = kwargs.get('label_values', True)
        link = kwargs.get('link', True)

        p1, p2, p3, p4 = [sch.nodes[n].pos for n in self.nodes]

        xoffset = 0.06
        yoffset = 0.40

        primary_dot = Pos(p3.x - xoffset, 0.5 * (p3.y + p4.y) + yoffset)
        secondary_dot = Pos(p1.x + xoffset, 0.5 * (p1.y + p2.y) + yoffset)

        centre = Pos(0.5 * (p3.x + p1.x), 0.5 * (p2.y + p1.y))
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


    def draw(self, sch, **kwargs):

        if not self.right:
            raise ValueError('Cannot draw transformer %s in direction %s'
                             % (self.name, self.opts['dir']))

        n1, n2, n3, n4 = self.nodes

        s = r'  \draw (%s) to [inductor] (%s);''\n' % (n3, n4)
        s += r'  \draw (%s) to [inductor] (%s);''\n' % (n1, n2)

        s += super(TF, self).draw(sch, link=True, **kwargs)

        draw_nodes = kwargs.get('draw_nodes', True)
        s += self._draw_nodes(sch, draw_nodes)
        return s


class K(TF1):
    """Mutual coupling"""

    def fixup(self, sch):
        
        L1 = sch.elements[self.Lname1]
        L2 = sch.elements[self.Lname2]

        self.nodes = L2.nodes + L1.nodes

    def draw(self, sch, **kwargs):

        if not self.right:
            raise ValueError('Cannot draw mutual coupling %s in direction %s'
                             % (self.name, self.opts['dir']))

        return super(K, self).draw(sch, link=True, **kwargs)


class OnePort(Cpt):
    """OnePort"""

    pos = ((0, 0), (1, 0))

    def draw(self, sch, **kwargs):

        label_values = kwargs.get('label_values', True)
        draw_nodes = kwargs.get('draw_nodes', True)
        label_ids = kwargs.get('label_ids', True)

        n1, n2 = self.nodes[0:2]

        tikz_cpt = self.tikz_cpt
        if self.type == 'R' and self.variable:
            tikz_cpt = 'vR'

        id_pos = '_'
        voltage_pos = '^'
        if self.type in ('V', 'Vdc', 'Vstep', 'Vac', 'Vacstep', 'Vimpulse', 'v',
                            'I', 'Idc', 'Istep', 'Iac', 'Iacstep', 'Iimpulse', 'i',
                            'E', 'F', 'G', 'H', 'Vs', 'Is'):

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

        node_str = self._node_str(sch.nodes[n1], sch.nodes[n2],
                                  draw_nodes)

        args_str += voltage_str + current_str

        # Generate default label unless specified.
        if label_str == '':
            if self.type not in ('O', 'W'):
                
                label_str = ', l%s=${%s}$' % (id_pos, self.default_label)
                
                if label_ids and self.value_label != '':
                    label_str = r', l%s=${%s=%s}$' % (id_pos, self.id_label, self.value_label)
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

    def fixup(self, sch):
        
        self.pos = self.npos if self.mirror else self.ppos

    def draw(self, sch, **kwargs):

        draw_nodes = kwargs.get('draw_nodes', True)
        label_values = kwargs.get('label_values', True)

        if not self.right:
            raise ValueError('Cannot draw opamp %s in direction %s'
                             % (self.name, self.opts['dir']))

        p1, p2, p3, p4 = [sch.nodes[n].pos for n in self.nodes]

        centre = Pos(0.5 * (p3.x + p1.x), p1.y)

        label_str = '$%s$' % self.default_label if label_values else ''
        args_str = '' if self.mirror else 'yscale=-1'
        for key, val in self.opts.iteritems():
            if key in ('color', ):
                args_str += '%s=%s, ' % (key, val)                

        s = r'  \draw (%s) node[op amp, %s, scale=%.1f] (opamp) {};' % (
            centre, args_str, sch.scale * 2)
        # Draw label separately to avoid being scaled by 2.
        s += r'  \draw (%s) node [] {%s};' % (centre, label_str)
        
        s += self._draw_nodes(sch, draw_nodes)
        return s


class FDOpamp(Cpt):

    pos = ((2, 1), (2, 0), (0, 1), (0, 0))

    def draw(self, sch, **kwargs):

        draw_nodes = kwargs.get('draw_nodes', True)
        label_values = kwargs.get('label_values', True)

        if not self.right:
            raise ValueError('Cannot draw fdopamp %s in direction %s'
                             % (self.name, self.opts['dir']))

        p1, p2, p3, p4 = [sch.nodes[n].pos for n in self.nodes]

        centre = (p1 + p4) * 0.5 + 0.2

        label_str = '$%s$' % self.default_label if label_values else ''
        args_str = '' if self.mirror else 'yscale=-1'
        for key, val in self.opts.iteritems():
            if key in ('color', ):
                args_str += '%s=%s, ' % (key, val)                

        s = r'  \draw (%s) node[fd op amp, %s, scale=%.1f] (opamp) {};' % (
            centre, args_str, sch.scale * 2)
        # Draw label separately to avoid being scaled by 2.
        s += r'  \draw (%s) node [] {%s};' % (centre, label_str)
        
        s += self._draw_nodes(sch, draw_nodes)
        return s


class Wire(OnePort):


    def draw_implicit(self, sch, **kwargs):

        # Draw implict wires, i.e., connections to ground, etc.

        kind = ''
        if 'implicit' in self.opts:
            kind = 'sground'
        elif 'ground' in self.opts:
            kind = 'ground'
        elif 'sground' in self.opts:
            kind = 'sground'

        args = [kind]
        if self.up:
            args.append('yscale=-1')
        if self.left:
            args.append('xscale=-1')
        args_str = ','.join(args)

        offset = (0.0, 0.0)
        anchor = 'south west'
        p = sch.nodes[n2].pos
        if self.down:
            offset = (0.0, 0.25)
            anchor = 'north west'
        elif self.up:
            offset = (0.0, -0.25)
        pos = p + offset

        s = r'  \draw (%s) to [short] (%s);''\n' % (n1, pos)
        s += r'  \draw (%s) node[%s] {};''\n' % (pos, args_str)

        if 'l' in self.opts:
            label_str = '${%s}$' % latex_str(self.opts['l'])
            s += r'  \draw {[anchor=%s] (%s) node {%s}};''\n' % (anchor, n2, label_str)
        return s

    def draw(self, sch, **kwargs):

        if ('implicit' in self.opts or 'ground' in self.opts
            or 'sground' in self.opts):
            return self.draw_implicit(nodes)
                                    
        return (super, Wire).draw(**kwargs)


classes = {}

def defcpt(name, base, docstring, cpt=None):
    
    if isinstance(base, str):
        base = classes[base]

    newclass = type(name, (base, ), {'__doc__': docstring})

    if cpt is not None:
        newclass.tikz_cpt = cpt
    classes[name] = newclass

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
defcpt('F', VCS, 'VCCS', 'american controlled current source')
defcpt('G', CCS, 'CCVS', 'american controlled voltage source')
defcpt('H', CCS, 'CCCS', 'american controlled current source')

defcpt('I', OnePort, 'Current source', 'I')
defcpt('Isin', 'I', 'Sinusoidal current source', 'sI')
defcpt('Idc', 'I', 'DC current source', 'I')
defcpt('Iac', 'I', 'AC current source', 'sI')

defcpt('J', JFET, 'N JFET transistor', 'njfet')
defcpt('Jnjf', 'J', 'N JFET transistor', 'njfet')
defcpt('Jpjf', 'J', 'P JFET transistor', 'pjfet')

defcpt('L', OnePort, 'Inductor', 'L')

defcpt('M', 'J', 'N MOSJFET transistor', 'nmos')
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
defcpt('Vsin', 'V', 'Sinusoidal voltage source', 'sV')
defcpt('Vdc', 'V', 'DC voltage source', 'V')
defcpt('Vac', 'V', 'AC voltage source', 'sV')

defcpt('W', OnePort, 'Wire', 'short')
defcpt('Y', OnePort, 'Admittance', 'european resistor')
defcpt('Z', OnePort, 'Impedance', 'european resistor')

# Add TP and TF.
# Perhaps AM for ammeter, VM for voltmeter, VR for variable resistor?
# Currently, a variable resistor is supported with the variable
# option.

