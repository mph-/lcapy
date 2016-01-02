from __future__ import print_function
from schemmisc import *
import numpy as np


class Cpt(object):

    pos = {}
    yscale = 1.0
    xscale = 1.0
    cpt_type_counter = 1

    pos = {1: (0, 0), 2: (0, 1)}

    @property
    def num_nodes(self):
        return len(self.pos)

    @property
    def xvals_flip(self):
        return self.xmin + self.xmax - self.xvals

    @property
    def yvals_flip(self):
        return self.ymin + self.ymax - self.yvals

    @property
    def xmin(self):
        return self.xvals.min()

    @property
    def xmax(self):
        return self.xvals.max()

    @property
    def ymin(self):
        return self.yvals.min()

    @property
    def ymax(self):
        return self.yvals.max()

    @property
    def xextent(self):
        return self.xmax - self.xmin

    @property
    def yextent(self):
        return self.ymax - self.ymin

    def __init__(self, name, *args, **kwargs):

        self.name = name
        self.args = args
        self.xvals = np.array(self.pos.values())[:, 0]
        self.yvals = np.array(self.pos.values())[:, 1]
            
        if True:
            return

        nodes = kwargs['nodes']
        if len(nodes) != len(self.xvals):
            raise ValueError('Expecting %d nodes, got %d' % 
                             (len(self.xvals), len(nodes)))
        self.nodes = nodes

        if hasattr(kwargs, 'opts'):
            opts = kwargs['opts']
        else:
            opts = {}
        if 'dir' not in opts:
            opts['dir'] = None
        if 'size' not in opts:
            opts['size'] = 1

        if opts['dir'] is None:
            opts['dir'] = 'down' if isinstance(self, (Open, Port)) else 'right'
        self.opts = opts

        if hasattr(self, 'check'):
            self.check()

    def __repr__(self):

        if hasattr(self, 'string'):
            return self.string
        
        return type(self)

    def draw(self, label_values=True, draw_nodes=True):
        pass


class Transistor(Cpt):
    """Transistor"""
    
    yscale = 1.5
    xscale = 0.85

    pos = {1: (1, 1), 2: (0, 0.5), 3: (1, 0)}


    def check(self):

        # For common base, will need to support up and down.
        if self.opts['dir'] not in ('left', 'right'):
            raise ValueError('Cannot draw transistor %s in direction %s'
                             '; try left or right'
                             % (self.name, self.opts['dir']))

    def draw(self, label_values=True, draw_nodes=True):

        n1, n2, n3 = self.nodes
        p1, p2, p3 = self.coords

        centre = Pos(p3.x, 0.5 * (p1.y + p3.y))

        label_str = '$%s$' % self.default_label if label_values else ''
        args_str = '' if self.opts['dir'] == 'right' else 'xscale=-1'
        if 'mirror' in self.opts:
            args_str += ', yscale=-1'
        for key, val in self.opts.iteritems():
            if key in ('color', ):
                args_str += '%s=%s, ' % (key, val)                

        s = r'  \draw (%s) node[%s, %s, scale=%.1f] (T) {};''\n' % (
            centre, tikz_cpt, args_str, self.scale * 2)
        s += r'  \draw (%s) node [] {%s};''\n'% (centre, label_str)

        if tikz_cpt in ('pnp', 'pmos', 'pjfet'):
            n1, n3 = n3, n1

        # Add additional wires.
        if tikz_cpt in ('pnp', 'npn'):
            s += r'  \draw (T.C) -- (%s) (T.B) -- (%s) (T.E) -- (%s);''\n' % (n1, n2, n3)
        else:
            s += r'  \draw (T.D) -- (%s) (T.G) -- (%s) (T.S) -- (%s);''\n' % (n1, n2, n3)


        # s += self._tikz_draw_nodes(self, draw_nodes)
        return s


class TwoPort(Cpt):
    """Two-port"""

    pos = {1: (0, 0), 2: (0, 1), 3: (1, 0), 4: (1, 1)}


class TF(TwoPort):
    """Transformer"""


class OnePort(Cpt):
    """OnePort"""

    # horiz, need to rotate for up/down
    pos = {1: (0, 0), 2: (1, 0)}


class VCS(OnePort):
    """Voltage controlled source"""


class CCS(OnePort):
    """Current controlled source"""


class K(Cpt):
    """K"""


classes = {}

def defcpt(name, base, docstring, cpt=None, pos=None):
    
    if isinstance(base, str):
        base = classes[base]

    newclass = type(name, (base, ), {'__doc__': docstring})

    if cpt is not None:
        newclass.tikz_cpt = cpt
    if pos is not None:
        newclass.pos= pos
    classes[name] = newclass

# Dynamically create classes.

defcpt('C', OnePort, 'Capacitor', 'C')

defcpt('D', OnePort, 'Diode', 'I')
defcpt('Dled', 'D', 'LED', 'leD')
defcpt('Dphoto', 'D', 'Photo diode', 'pD')
defcpt('Dshottky', 'D', 'Shottky diode', 'zD')
defcpt('Dtunnel', 'D', 'Tunnel diode', 'tD')
defcpt('Dzener', 'D', 'Zener diode', 'zD')

defcpt('E', VCS, 'VCVS', 'V')
defcpt('Eopamp', 'E', 'VCVS', 'V')
defcpt('F', VCS, 'VCCS', 'I')
defcpt('G', CCS, 'CCVS', 'V')
defcpt('H', CCS, 'CCCS', 'I')

defcpt('I', OnePort, 'Current source', 'I')
defcpt('Isin', 'I', 'Sinusoidal current source', 'sI')
defcpt('Idc', 'I', 'DC current source', 'I')
defcpt('Iac', 'I', 'AC current source', 'I')

defcpt('J', Transistor, 'N JFET transistor', 'njft', {1: (1, 1), 2: (0, 0.5), 3: (1, 0)})
defcpt('Jnjf', 'J', 'N JFET transistor', 'njf')
defcpt('Jpjf', 'J', 'P JFET transistor', 'pjf')

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

defcpt('SW', OnePort, 'SWitch', 'SW')
defcpt('SWno', 'SW', 'Normally open switch', 'SWno')
defcpt('SWnc', 'SW', 'Normally closed switch', 'SWnc')
defcpt('SWpush', 'SW', 'Pushbutton switch', 'SWpush')

defcpt('V', OnePort, 'Voltage source', 'V')
defcpt('Vsin', 'V', 'Sinusoidal voltage source', 'sV')
defcpt('Vdc', 'V', 'DC voltage source', 'V')
defcpt('Vac', 'V', 'AC voltage source', 'V')

defcpt('W', OnePort, 'Wire', 'short')

# Add TP and TF.
# Perhaps AM for ammeter, VM for voltmeter, VR for variable resistor?
# Currently, a variable resistor is supported with the variable
# option.

