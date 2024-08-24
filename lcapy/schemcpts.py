"""This module defines and draws the schematic components using
circuitikz.   The components are defined at the bottom of this file.

Copyright 2015--2024 Michael Hayes, UCECE
"""

from __future__ import print_function
from warnings import warn
import sys
import numpy as np
from .schemmisc import Pos, Steps
from .label import Label
from .labels import Labels
from .latex import latex_format_label
from .config import implicit_default
from .schematics.components.cpt import Cpt
from .schematics.utils import check_boolean
from .schematics.components.bipole import Bipole
from .schematics.components.cable import Cable
from .schematics.components.fixedcpt import FixedCpt
from .schematics.components.shape import Shape
from .schematics.components.stretchycpt import StretchyCpt
from .schematics.components.transistor import Transistor
from .schematics.components.twoport import TwoPort
from .schematics.components.unipole import Unipole
from .schematics.components.wire import Wire



module = sys.modules[__name__]

# There are two types of component (Cpt):
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

# The node_pinnames attribute specifies the subset of required nodes.
# It is a list of pinnames that are used to find the pin coordinates
# and thus the node coordinates.  Note, some components do not have
# any explicit nodes (shapes, chips, etc).

# The direction commands rotate the component:
# right (0 degrees)
# up    (90 degrees)
# left  (180 degrees)
# down  (-90 degrees)
#
# Lcapy uses this information to position the second node with respect
# to the first.   Circuitikz then infers the rotation.

# `node_names` is a list of the node names specified for each component.
# For example, given `R1 1 2` the node_names for R1 are ['1', '2]
# '
# `required_node_names` is a list comprising a subset of node_names,
# ignoring the nodes that are not drawn, say the ground node for an opamp.
#
# `auxiliary_node_names` is a list of node names for each component used
# to draw the component.
#
# `pin_node_names` is a list of the pin nodes for each component.
# For example, ['R1.p', 'R1.n'] or ['U1.t1', 'U1.l1', 'U1.b1', 'U1.r1'].
# This does not include the aliases such as `U1.vdd`.
#
# `attribute_node_names` is a list of node names specified as
# schematic attributes.  For example, given `R1 1 2; .p.l=foo` the
# list is ['.p'].  CHECKME
#
# `all_node_names` is the union of required_node_names,
# auxiliary_node_names, pin_node_names, and attribute_node_names
#
# `relative_node_names` is a list of node names defined in a component
# net with a dot prefix, for example, `R1 .a .b`.  This is shorthand
# for `R1 R1.a R1.b`.  The resulting list is ['R1.a', 'R1.b'].
#
# `U1 chip2121; W 1 U1.vdd` Here U1.vdd is an alias for U1.t1.
#
# The node linking compares xvals and yvals.  These are derived from
# transformed coords (tcoords).  coords uses required_pins which uses nodes.
# Finally, nodes is a subset of all_node_names selected if the node
# name is explicitly registered.

# TODO, this has evolved into a can of words and could do with a
# complete rewrite.

# Note, the tikz origin is in the lower left corner


class A(Unipole):
    """Annotation."""

    tikz_cpt = ''


class ANT(Unipole):
    """Antenna"""

    can_mirror = True
    can_invert = True
    tikz_cpt = 'antenna'
    kinds = {'rx': 'rxantenna', 'tx': 'txantenna'}


class BAT(Bipole):
    """Battery"""

    tikz_cpt = 'battery'
    kinds = {'cell1': 'battery1'}


class BL(Bipole):
    """Block"""

    tikz_cpt = 'twoport'

    kinds = {'vco': 'vco', 'bandpass': 'bandpass', 'bandstop':
             'bandstop', 'highpass': 'highpass', 'lowpass': 'lowpass',
             'allpass': 'allpass', 'highpass2': 'highpass2',
             'lowpass2': 'lowpass2', 'adc': 'adc', 'dac': 'dac',
             'dsp': 'dsp', 'fft': 'fft', 'amp': 'amp', 'vamp': 'vamp',
             'phaseshifter': 'phaseshifter', 'vphaseshifter':
             'vphaseshifter', 'piattenuator': 'piattenuator',
             'vpiattenuator': 'vpiattenuator', 'tattenuator':
             'tattenuator', 'vtattenuator': 'vtattenuator', 'dcdc':
             'sdcdc', 'dcac': 'sdcac', 'acdc': 'sacdc', 'detector':
             'detector', 'twoport': 'twoport', 'twoportsplit': 'twoportsplit'}


class C(Bipole):
    """Capacitor"""

    tikz_cpt = 'C'
    kinds = {'electrolytic': 'eC', 'polar': 'cC',
             'variable': 'vC', 'sensor': 'sC', 'ferroelectric': 'ferrocap',
             'curved': 'cC', 'tunable': 'vC, tunable end arrow={Bar}'}


class D(Bipole):
    """Diode"""

    tikz_cpt = 'D'
    kinds = {'led': 'leD', 'photo': 'pD', 'schottky': 'sD',
             'zener': 'zD', 'zzener': 'zzD', 'tunnel': 'tD', 'varcap': 'VC',
             'bidirectional': 'biD', 'tvs': 'tvsDo', 'laser': 'lasD'}
    styles = {'empty': '', 'full': '*', 'stroke': '-'}


class I(Bipole):
    """Current source"""

    tikz_cpt = 'I'
    kinds = {'dc': 'I', 'step': 'I', 'ac': 'sI',
             'square': 'sqI', 'triangle': 'tI', 'noise': 'nI'}


class L(Bipole):
    """Inductor"""

    tikz_cpt = 'L'
    kinds = {'variable': 'vL', 'choke': 'cute choke',
             'twolineschoke': 'cute choke, twolineschoke', 'sensor': 'sL',
             'tunable': 'vL, tunable end arrow={Bar}'}


class R(Bipole):
    """Resistor"""

    tikz_cpt = 'R'
    kinds = {'variable': 'vR', 'tunable': 'vR, tunable end arrow = {Bar}'}


class REL(Bipole):
    """Reluctance"""

    tikz_cpt = 'R'
    kinds = {'variable': 'vR', 'tunable': 'vR, tunable end arrow = {Bar}'}


class V(Bipole):
    """Voltage source"""

    tikz_cpt = 'V'
    kinds = {'dc': 'V', 'step': 'V', 'ac': 'sV',
             'square': 'sqV', 'triangle': 'tV', 'noise': 'nV'}


class Y(Bipole):
    """Admittance"""

    tikz_cpt = 'generic'
    kinds = {'variable': 'variable european resistor',
             'tunable': 'variable european resistor, tunable end arrow = {Bar}',
             'sensor': 'european resistive sensor'}


class Z(Bipole):
    """Impedance"""

    tikz_cpt = 'generic'
    kinds = {'variable': 'variable european resistor',
             'tunable': 'variable european resistor, tunable end arrow = {Bar}',
             'sensor': 'european resistive sensor'}


class BJT(Transistor):
    """BJT"""

    node_pinnames = ('c', 'b', 'e')
    aliases = {'emitter': 'e', 'base': 'b', 'collector': 'c'}
    normal_pins = {'e': ('lx', 0.55, 0),
                   'b': ('lx', 0, 0.5),
                   'c': ('lx', 0.55, 1)}
    mirror_pins = {'e': ('lx', 0.55, 1),
                   'b': ('lx', 0, 0.5),
                   'c': ('lx', 0.55, 0)}
    invert_pins = {'e': ('lx', 0, 0),
                   'b': ('lx', 0.55, 0.5),
                   'c': ('lx', 0, 1)}
    mirror_invert_pins = {'e': ('lx', 0, 1),
                          'b': ('lx', 0.55, 0.5),
                          'c': ('lx', 0, 0)}
    kinds = {'npn': 'npn',
             'pnp': 'pnp',
             'bodydiode': '-bodydiode',
             'nigbt': 'nigbt',
             'pigbt': 'pigbt',
             'nigbt-bodydiode': 'nigbt-bodydiode',
             'pigbt-bodydiode': 'pigbt-bodydiode',
             'Lnigbt': 'Lnigbt',
             'Lpigbt': 'Lpigbt',
             'Lnigbt-bodydiode': 'Lnigbt-bodydiode',
             'Lpigbt-bodydiode': 'Lpigbt-bodydiode'}


class JFET(Transistor):
    """JFET"""

    node_pinnames = ('d', 'g', 's')
    aliases = {'drain': 'd', 'gate': 'g', 'source': 's'}
    normal_pins = {'d': ('lx', 0.55, 1),
                   'g': ('lx', 0, 0.335),
                   's': ('lx', 0.55, 0)}
    mirror_pins = {'d': ('lx', 0.55, 0),
                   'g': ('lx', 0, 0.645),
                   's': ('lx', 0.55, 1)}
    invert_pins = {'d': ('lx', 0, 1),
                   'g': ('lx', 0.55, 0.335),
                   's': ('lx', 0, 0)}
    mirror_invert_pins = {'d': ('lx', 0, 0),
                          'g': ('lx', 0.55, 0.645),
                          's': ('lx', 0, 1)}


class MOSFET(Transistor):
    """MOSFET"""

    node_pinnames = ('d', 'g', 's')
    aliases = {'drain': 'd', 'gate': 'g', 'source': 's'}
    normal_pins = {'d': ('lx', 0.55, 1),
                   'g': ('lx', 0, 0.5),
                   's': ('lx', 0.55, 0)}
    mirror_pins = {'d': ('lx', 0.55, 0),
                   'g': ('lx', 0, 0.5),
                   's': ('lx', 0.55, 1)}
    invert_pins = {'d': ('lx', 0, 1),
                   'g': ('lx', 0.55, 0.5),
                   's': ('lx', 0, 0)}
    mirror_invert_pins = {'d': ('lx', 0, 0),
                          'g': ('lx', 0.55, 0.5),
                          's': ('lx', 0, 1)}
    normal_pins2 = {'d': ('lx', 0.55, 1),
                    'g': ('lx', 0, 0.355),
                    's': ('lx', 0.55, 0)}
    mirror_pins2 = {'d': ('lx', 0.55, 0),
                    'g': ('lx', 0, 0.645),
                    's': ('lx', 0.55, 1)}
    invert_pins2 = {'d': ('lx', 0, 1),
                    'g': ('lx', 0.55, 0.335),
                    's': ('lx', 0, 0)}
    mirror_invert_pins2 = {'d': ('lx', 0, 0),
                           'g': ('lx', 0.55, 0.645),
                           's': ('lx', 0, 1)}
    kinds = {'nmos': 'nmos', 'pmos': 'pmos',
             'nmosd': 'nmosd', 'pmosd': 'pmosd',
             'nfet': 'nfet', 'pfet': 'pfet',
             'nfetd': 'nfetd', 'pfetd': 'pfetd',
             'nfet-bodydiode': 'nfet-bodydiode',
             'pfet-bodydiode': 'pfet-bodydiode',
             'nfetd-bodydiode': 'nfetd-bodydiode',
             'pfetd-bodydiode': 'pfetd-bodydiode',
             'nigfetd': 'nigfetd',
             'pigfetd': 'pigfetd',
             'nigfetd-bodydiode': 'nigfetd-bodydiode',
             'pigfetd-bodydiode': 'pigfetd-bodydiode',
             'nigfete': 'nfigete',
             'pigfete': 'pigfete',
             'nigfete-bodydiode': 'nfigete-bodydiode',
             'pigfete-bodydiode': 'pigfete-bodydiode',
             'nigfetebulk': 'nigfetebulk', 'pigfetebulk': 'pigfetebulk',
             'hemt': 'hemt'}


class MT(Bipole):
    """Motor"""

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.nodes
        centre = (n1.pos + n2.pos) * 0.5

        s = r'  \draw (%s) node[elmech, %s, rotate=%d] (%s) {};''\n' % (
            centre, self.draw_args_str(**kwargs), self.angle + 90, self.s)
        # Draw label separately, shape border rotate does not seem to work
        s += self.draw_label(centre, **kwargs)
        s += r'  \draw (%s) |- (%s.north);''\n' % (n1.s, self.s)
        s += r'  \draw (%s.south) |- (%s);''\n' % (self.s, n2.s)
        return s


class MX(FixedCpt):
    """Mixer"""

    # Dubious
    can_scale = True

    node_pinnames = ('in1', 'in2', 'out')
    pins = {'in1': ('lx', 0.25, 0.25),
            'in2': ('lx', -0.25, 0.25),
            'out': ('rx', 0, 0)}

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

    node_pinnames = ('in1', 'in2', 'out', 'in3')
    normal_pins = {'in1': ('lx', -0.25, 0),
                   'in2': ('bx', 0, -0.25),
                   'out': ('rx', 0.25, 0),
                   'in3': ('tx', 0, 0.25)}
    mirror_pins = {'in1': ('lx', -0.25, 0),
                   'in2': ('tx', 0, 0.25),
                   'out': ('rx', 0.25, 0),
                   'in3': ('bx', 0, -0.25)}

    @property
    def pins(self):
        return self.mirror_pins if self.mirror else self.normal_pins

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
        s += r'  \draw (%s) node[] {$%s$};''\n' % (q[1], self.splabels[0])

        if self.mirror:
            s += r'  \draw (%s) node[] {$%s$};''\n' % (q[4], self.splabels[0])
        else:
            s += r'  \draw (%s) node[] {$%s$};''\n' % (q[2], self.splabels[1])
        if len(self.splabels) > 2:
            s += r'  \draw (%s) node[] {$%s$};''\n' % (q[3], self.splabels[2])
        return s


class SP3(SP):
    """Summing point"""

    node_pinnames = ('in1', 'in2', 'out')
    normal_pins = {'in1': ('lx', -0.25, 0),
                   'in2': ('bx', 0, -0.25),
                   'out': ('rx', 0.25, 0)}
    mirror_pins = {'in1': ('lx', -0.25, 0),
                   'in2': ('tx', 0, 0.25),
                   'out': ('rx', 0.25, 0)}

    @property
    def pins(self):
        return self.mirror_pins if self.mirror else self.normal_pins


class SPpp(SP3):
    """Summing point"""

    splabels = '++'


class SPpm(SP3):
    """Summing point"""

    splabels = '+-'


class SPppp(SP):
    """Summing point"""

    splabels = '+++'


class SPpmm(SP):
    """Summing point"""

    splabels = '+--'


class SPppm(SP):
    """Summing point"""

    splabels = '++-'


class TL(StretchyCpt):
    """Transmission line"""

    # Scaling is dubious.  Perhaps should stretch this
    # component in proportion to size?  Applying an xscale without a
    # corresponding yscale changes the ellipse.  This should be fixed
    # in circuitikz.
    can_scale = True
    w = 1

    node_pinnames = ('out1', 'out2', 'in1', 'in2')
    aliases = {'out+': 'out1', 'out-': 'out2', 'in+': 'in1', 'in-': 'in2'}
    pins = {'in1': ('lx', 0, 0.5),
            'in2': ('lx', 0, 0),
            'out1': ('rx', w, 0.5),
            'out2': ('rx', w, 0)}

    @property
    def drawn_nodes(self):

        if self.nowires:
            return self.nodes[0:1]
        return self.nodes

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.nodes
        centre = (n1.pos + n3.pos) * 0.5 + Pos(1 - self.w, 0)

        q = self.xtf(centre, ((0.23, 0),
                              (0.20, -0.1),
                              (-0.35, 0),
                              (-0.27, -0.11)))
        pos = n3.pos

        # Newer versions of circuitikz have a tline component with
        # wires on each end.  Rotation creates an ellipse!
        s = r'  \draw (%s) node[tlinestub, xscale=%s, rotate=%s] {};''\n' % (
            pos, self.scale, self.angle)

        s += self.draw_label(centre, **kwargs)

        # Output wire
        s += self.draw_path((q[0], n1.s))
        # Extra input wire
        s += self.draw_path((q[2], n3.s))
        if not self.nowires:
            # Output ground wire
            s += self.draw_path((q[1], n2.s), join='|-')
            # Input ground wire
            s += self.draw_path((q[3], n4.s), join='|-')
        return s


class TF1(FixedCpt):
    """Transformer"""

    can_scale = True
    w = 0.8
    default_aspect = w
    node_pinnames = ('s+', 's-', 'p+', 'p-')
    pins = {'s+': ('rx', w, 1),
            's-': ('rx', w, 0),
            'p+': ('lx', 0, 1),
            'p-': ('lx', 0, 0)}
    misc = {'pdot': (0.1 - 0.5 * w, 0.34),
            'sdot': (0.5 * w - 0.1, 0.34),
            'link': (0, 0.15),
            'label': (0, 0.48)}

    def draw(self, link=True, **kwargs):

        if not self.check():
            return ''

        centre = self.midpoint(self.nodes[0], self.nodes[3])
        pdot_pos = self.tf(centre, self.misc['pdot'])
        sdot_pos = self.tf(centre, self.misc['sdot'])
        label_pos = self.tf(centre, self.misc['label'])

        s = ''
        if not self.nodots:
            s += r'  \draw (%s) node[circ] {};''\n' % pdot_pos
            s += r'  \draw (%s) node[circ] {};''\n' % sdot_pos

        s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
            label_pos, 0.5, self.s, self.label(**kwargs))

        if link:
            # TODO: allow for rotation
            width = (sdot_pos - pdot_pos).x
            link_pos = self.tf(centre, self.misc['link'])

            s += r'  \draw[<->] ([shift=(45:%.2f)]%s) arc(45:135:%.2f);''\n' % (
                width / 2, link_pos, width / 2)

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

        n1, n2, n3, n4 = self.nodes[0:4]

        draw_args_str = self.draw_args_str(**kwargs)

        s = self.draw_cpt(n3.s, n4.s, cpt='inductor', args='scale=%s' %
                          self.scale, dargs=draw_args_str)
        s += self.draw_cpt(n2.s, n1.s, cpt='inductor', args='scale=%s' %
                          self.scale, dargs=draw_args_str)

        s += super(Transformer, self).draw(link=False, **kwargs)
        return s


class TFtap(TF1):
    """Transformer"""

    node_pinnames = ('s+', 's-', 'p+', 'p-', 'ptap', 'stap')
    w = 0.8
    pins = {'s+': ('rx', w, 1),
            's-': ('rx', w, 0),
            'p+': ('lx', 0, 1),
            'p-': ('lx', 0, 0),
            'ptap': ('lx', 0, 0.5),
            'stap': ('rx', w, 0.5)}

    @property
    def drawn_nodes(self):
        # Do not draw the taps.
        return self.nodes[0:4]

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.nodes[0:4]

        draw_args_str = self.draw_args_str(**kwargs)

        s = self.draw_cpt(n3.s, n4.s, cpt='inductor', args='scale=%s' %
                          self.scale, dargs=draw_args_str)
        s += self.draw_cpt(n2.s, n1.s, cpt='inductor', args='scale=%s' %
                          self.scale, dargs=draw_args_str)

        s += super(TFtap, self).draw(link=False, **kwargs)
        return s


class K(TF1):
    """Mutual coupling"""

    def __init__(self, sch, namespace, name, cpt_type, cpt_id, string,
                 opts_string, node_names, keyword, *args):

        self.Lname1 = args[0]
        self.Lname2 = args[1]
        super(K, self).__init__(sch, namespace, name,
                                cpt_type, cpt_id, string, opts_string,
                                node_names, keyword, *args[2:])

    @property
    def coords(self):
        return ((0.5, 1), (0.5, 0), (0, 1), (0, 0))

    @property
    def scales(self):
        return [self.scale] * len(self.coords)

    @property
    def nodes(self):

        # This needs to be more robust; currently it depends on the order
        # that the inductors are defined.

        # L1 and L2 need to be previously defined so we can find their nodes.
        L1 = self.sch.elements[self.Lname1]
        L2 = self.sch.elements[self.Lname2]
        # L1 is on the left; L2 is on the right
        nodes = [self.sch.nodes[n] for n in L2.node_names + L1.node_names]
        return nodes


class Gyrator(FixedCpt):
    """Gyrator"""

    node_pinnames = ('out+', 'out-', 'in+', 'in-')
    pins = {'out+': ('rx', 1, 1),
            'out-': ('rx', 1, 0),
            'in+': ('lx', 0, 1),
            'in-': ('lx', 0, 0)}

    def draw(self, **kwargs):

        if not self.check():
            return ''

        yscale = self.scale
        if not self.mirror:
            yscale = -yscale

        s = r'  \draw (%s) node[gyrator, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            self.midpoint(self.nodes[0], self.nodes[3]),
            self.draw_args_str(**kwargs), 0.95 * self.scale, 0.89 * yscale,
            -self.angle, self.s)

        s += self.draw_label(self.centre, **kwargs)
        return s


class Triode(FixedCpt):
    """Triode"""

    node_pinnames = ('a', 'g', 'k')
    aliases = {'anode': 'a', 'grid': 'g', 'cathode': 'k'}
    pins = {'a': ('l', 0.75, 0),
            'g': ('l', 0.25, 0.5),
            'k': ('l', -0.25, 0)}

    def draw(self, **kwargs):

        if not self.check():
            return ''

        yscale = self.scale
        if self.mirror:
            yscale = -yscale

        # stupid thing above sets distance between nodes.
        # get correct distance, then slide into place.
        mid = self.centre
        mid = mid + Pos(0, -0.5)

        s = r'  \draw (%s) node[triode, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            str(mid),
            self.draw_args_str(**kwargs), 1 * self.scale, 1 * yscale,
            0, self.s)

        s += self.draw_label(self.centre, **kwargs)
        return s


class Potentiometer(Bipole):
    """Potentiometer  Np, Nm, No"""

    # This is not really a bipole but circuitikz treats it as such

    can_stretch = False

    node_pinnames = ('p', 'n', 'wiper')
    aliases = {'+': 'p', '-': 'n'}
    pins = {'p': ('rx', 0, 0),
            'n': ('rx', 1, 0),
            'wiper': ('lx', 0.5, 0.3)}


class VCS(Bipole):
    """Voltage controlled source"""

    node_pinnames = ('+', '-', '', '')


class CCS(Bipole):
    """Current controlled source"""

    node_pinnames = ('+', '-', '', '')


class SPDT(FixedCpt):
    """SPDT switch"""

    can_mirror = True
    can_invert = True

    node_pinnames = ('p', 'n', 'common')
    aliases = {'+': 'p', '-': 'n'}
    normal_pins = {'p': ('lx', 0, 0.169),
                   'n': ('rx', 0.632, 0.338),
                   'common': ('lx', 0.632, 0)}
    invert_pins = {'p': ('lx', 0.632, 0.169),
                   'n': ('rx', 0, 0.338),
                   'common': ('lx', 0, 0)}

    @property
    def pins(self):
        # The mirror pins are at the same locations as the normal pins
        return self.invert_pins if self.invert else self.normal_pins

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3 = self.nodes
        centre = n1.pos * 0.5 + (n2.pos + n3.pos) * 0.25

        args = self.draw_args(self.opts, **kwargs)
        if self.angle != 0:
            args.append('rotate=%d' % self.angle)

        xscale = self.scale
        yscale = self.scale
        if self.mirror:
            yscale = -yscale
        if self.invert:
            xscale = -xscale
        if xscale != 1:
            args.append('xscale=%s' % xscale)
        if yscale != 1:
            args.append('yscale=%s' % yscale)
        args.insert(0, 'spdt')

        s = r'  \draw (%s) node[%s] (%s) {};''\n' % (
            centre, ', '.join(args), self.s)

        # TODO, fix label position.
        centre = (n1.pos + n3.pos) * 0.5 + Pos(0.5, -0.5)
        s += self.draw_label(centre, **kwargs)
        return s


class Box2(Shape):
    """Square box,  A rectangle is created by defining aspect."""

    shape = 'rectangle'
    pins = {'w': ('l', -0.5, 0),
            'e': ('r', 0.5, 0)}


class Box4(Shape):
    """Box4"""

    shape = 'rectangle'
    pins = {'w': ('l', -0.5, 0),
            's': ('b', 0, -0.5),
            'e': ('r', 0.5, 0),
            'n': ('t', 0, 0.5)}


class Box12(Shape):
    """Box12"""

    shape = 'rectangle'
    pins = {'wnw': ('l', -0.5, 0.25),
            'w': ('l', -0.5, 0),
            'wsw': ('l', -0.5, -0.25),
            'ssw': ('b', -0.25, -0.5),
            's': ('b', 0, -0.5),
            'sse': ('b', 0.25, -0.5),
            'ese': ('r', 0.5, -0.25),
            'e': ('r', 0.5, 0),
            'ene': ('r', 0.5, 0.25),
            'nne': ('t', 0.25, 0.5),
            'n': ('t', 0, 0.5),
            'nnw': ('t', -0.25, 0.5)}


class Box(Shape):
    """Box"""

    shape = 'rectangle'
    pins = {'nw': ('t', -0.5, 0.5), 'wnw': ('l', -0.5, 0.25),
            'w': ('l', -0.5, 0), 'wsw': ('l', -0.5, -0.25),
            'sw': ('b', -0.5, -0.5), 'ssw': ('b', -0.25, -0.5),
            's': ('b', 0, -0.5), 'sse': ('b', 0.25, -0.5),
            'se': ('b', 0.5, -0.5), 'ese': ('r', 0.5, -0.25),
            'e': ('r', 0.5, 0), 'ene': ('r', 0.5, 0.25),
            'ne': ('t', 0.5, 0.5), 'nne': ('t', 0.25, 0.5),
            'n': ('t', 0, 0.5), 'nnw': ('t', -0.25, 0.5)}


class Ellipse(Shape):
    """Ellipse"""

    # Ellipse needs the tikz shapes library.
    shape = 'ellipse'
    pins = {'nw': ('t', -0.3536, 0.3536), 'wnw': ('l', -0.4619, 0.1913),
            'w': ('l', -0.5, 0), 'wsw': ('l', -0.4619, -0.1913),
            'sw': ('b', -0.3536, -0.3536), 'ssw': ('b', -0.1913, -0.4619),
            's': ('b', 0, -0.5), 'sse': ('b', 0.1913, -0.4619),
            'se': ('r', 0.3536, -0.3536), 'ese': ('r', 0.4619, -0.1913),
            'e': ('r', 0.5, 0), 'ene': ('r', 0.4619, 0.1913),
            'ne': ('r', 0.3536, 0.35365), 'nne': ('t', 0.1913, 0.4619),
            'n': ('t', 0, 0.5), 'nnw': ('t', -0.1913, 0.4619)}


class Circle(Ellipse):
    """Circle"""

    shape = 'circle'


class Circle2(Shape):
    """Circle"""

    shape = 'circle'
    pins = {'w': ('l', -0.5, 0),
            'e': ('r', 0.5, 0)}


class Circle4(Shape):
    """Circle4"""

    shape = 'circle'
    pins = {'w': ('l', -0.5, 0),
            's': ('b', 0, -0.5),
            'e': ('r', 0.5, 0),
            'n': ('t', 0, 0.5)}


class Triangle(Shape):
    """Equilateral triangle. The triangle shape can be altered by defining
    aspect."""

    shape = 'triangle'

    # 1 / sqrt(3) approx 0.5774, 1 / (2 * sqrt(3)) approx 0.2887
    pins = {'n': ('t', 0.0, 0.5774),
            'ne': ('r', 0.25, 0.14435),
            'nw': ('l', -0.25, 0.14435),
            'w': ('l', -0.5, -0.2887),
            'e': ('r', 0.5, -0.2887),
            's': ('b', 0.0, -0.2887),
            'se': ('b', 0.25, -0.2887),
            'sw': ('b', -0.25, -0.2887),
            'ssw': ('b', -0.125, -0.2887),
            'sse': ('b', 0.125, -0.2887),
            'nne': ('r', 0.125, 0.355),
            'nnw': ('l', -0.125, 0.355),
            'wsw': ('b', -0.375, -0.2887),
            'ese': ('b', 0.375, -0.2887),
            'ene': ('r', 0.375, -0.075),
            'wnw': ('l', -0.375, -0.075)}

    auxiliary = {'mid': ('c', 0.0, 0.0),
                 'bl': ('l', -0.5, -0.2887),
                 'br': ('r', 0.5, -0.2887),
                 'top': ('t', 0, 0.5774),
                 'tl': ('l', -0.5, 0.5774),
                 'tr': ('r', 0.5, 0.5774)}
    required_auxiliary = ('top', 'bl', 'br', 'mid')

    def draw(self, **kwargs):

        if not self.check():
            return ''

        s = self.draw_path([self.node('top').pos, self.node('bl').pos,
                            self.node('br').pos], closed=True, style='thick')
        s += self.draw_label(self.node('mid').pos, **kwargs)

        return s


class TR(Box2):
    """Transfer function"""

    default_width = 1.5
    default_aspect = 1.5
    node_pinnames = ('w', 'e')


class Chip(Shape):
    """General purpose chip"""

    do_transpose = True
    default_width = 2.0

    # Could allow can_scale but not a lot of point since nodes
    # will not be on the boundary of the chip.

    # TODO, tweak coord if pin name starts with \ using pinpos to
    # accomodate inverting circle.  This will require stripping of the
    # \ from the label. Alternatively, do not use inverting circle and
    # add overline to symbol name when printing.

    @property
    def path(self):
        return ((-0.5, 0.5), (0.5, 0.5), (0.5, -0.5), (-0.5, -0.5))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        centre = self.node('mid')

        s = ''
        if isinstance(self.path, list):
            # For multiple paths, e.g., isoamp
            for path in self.path:
                q = self.tf(centre.pos, path)
                s += self.draw_path(q, closed=True, style='thick')
        else:
            q = self.tf(centre.pos, self.path)
            s = self.draw_path(q, closed=True, style='thick')

        label = self.label(**kwargs)
        if label != '':
            s += r'  \draw (%s) node[text width=%.2fcm, align=center, %s] {%s};''\n' % (
                centre.s, self.width - 0.5, self.draw_args_str(**kwargs), label)

        # Draw clock symbols
        for m, n in enumerate(self.nodes):
            if n.clock:
                q = self.tf(n.pos, ((0, 0.125 * 0.707), (0.125, 0),
                                    (0, -0.125 * 0.707)))
                s += self.draw_path(q[0:3], style='thick')
        return s


class Uchip1313(Chip):
    """Chip of size 1 3 1 3"""

    pins = {'l1': ('l', -0.5, 0),
            'b1': ('b', -0.25, -0.5),
            'b2': ('b', 0, -0.5),
            'b3': ('b', 0.25, -0.5),
            'r1': ('r', 0.5, 0),
            't1': ('t', -0.25, 0.5),
            't2': ('t', 0, 0.5),
            't3': ('t', 0.25, 0.5)}
    aliases = {'vss': 'b2', 'vdd': 't2', 'in': 'l1', 'out': 'r1'}


class Uchip2121(Chip):
    """Chip of size 2 1 2 1"""

    pins = {'l1': ('l', -0.5, 0.25),
            'l2': ('l', -0.5, -0.25),
            'b1': ('b', 0, -0.5),
            'r2': ('r', 0.5, -0.25),
            'r1': ('r', 0.5, 0.25),
            't1': ('t', 0, 0.5)}
    aliases = {'vss': 'b1', 'vdd': 't1'}

    pinlabels = {'vss': 'VSS', 'vdd': 'VDD'}


class Uchip3131(Chip):
    """Chip of size 3 1 3 1"""

    pins = {'l1': ('l', -0.5, 0.25),
            'l2': ('l', -0.5, 0),
            'l3': ('l', -0.5, -0.25),
            'b1': ('b', 0.0, -0.5),
            'r3': ('r', 0.5, -0.25),
            'r2': ('r', 0.5, 0),
            'r1': ('r', 0.5, 0.25),
            't1': ('t', 0.0, 0.5)}
    aliases = {'vss': 'b1', 'vdd': 't1'}

    pinlabels = {'vss': 'VSS', 'vdd': 'VDD'}

    @property
    def path(self):
        return ((-0.5, 0.5), (0.5, 0.5), (0.5, -0.5), (-0.5, -0.5))


class Uchip3333(Chip):
    """Chip of size 3 3 3 3"""

    pins = {'l1': ('l', -0.5, 0.25),
            'l2': ('l', -0.5, 0),
            'l3': ('l', -0.5, -0.25),
            'b1': ('b', -0.25, -0.5),
            'b2': ('b', 0.0, -0.5),
            'b3': ('b', 0.25, -0.5),
            'r3': ('r', 0.5, -0.25),
            'r2': ('r', 0.5, 0),
            'r1': ('r', 0.5, 0.25),
            't1': ('t', -0.25, 0.5),
            't2': ('t', 0.0, 0.5),
            't3': ('t', 0.25, 0.5)}
    aliases = {'vss': 'b2', 'vdd': 't2'}

    pinlabels = {'vss': 'VSS', 'vdd': 'VDD'}

    @property
    def path(self):
        return ((-0.5, 0.5), (0.5, 0.5), (0.5, -0.5), (-0.5, -0.5))


class Uchip2222(Chip):
    """Chip of size 2 2 2 2"""

    pins = {'l1': ('l', -0.5, 0.25),
            'l2': ('l', -0.5, -0.25),
            'b1': ('b', -0.25, -0.5),
            'b2': ('b', 0.25, -0.5),
            'r1': ('r', 0.5, 0.25),
            'r2': ('r', 0.5, -0.25),
            't1': ('t', -0.25, 0.5),
            't2': ('t', 0.25, 0.5)}


class Uchip4141(Chip):
    """Chip of size 4 1 4 1"""

    pins = {'l1': ('l', -0.5, 0.375),
            'l2': ('l', -0.5, 0.125),
            'l3': ('l', -0.5, -0.125),
            'l4': ('l', -0.5, -0.375),
            'b1': ('b', 0.0, -0.5),
            'r4': ('r', 0.5, -0.375),
            'r3': ('r', 0.5, -0.125),
            'r2': ('r', 0.5, 0.125),
            'r1': ('r', 0.5, 0.375),
            't1': ('t', 0.0, 0.5)}
    aliases = {'vss': 'b1', 'vdd': 't1'}


class Uchip4444(Chip):
    """Chip of size 4 4 4 4"""

    pins = {'l1': ('l', -0.5, 0.375),
            'l2': ('l', -0.5, 0.125),
            'l3': ('l', -0.5, -0.125),
            'l4': ('l', -0.5, -0.375),
            'b1': ('b', -0.375, -0.5),
            'b2': ('b', -0.125, -0.5),
            'b3': ('b', .125, -0.5),
            'b4': ('b', 0.375, -0.5),
            'r4': ('r', 0.5, -0.375),
            'r3': ('r', 0.5, -0.125),
            'r2': ('r', 0.5, 0.125),
            'r1': ('r', 0.5, 0.375),
            't1': ('t', -0.375, 0.5),
            't2': ('t', -0.125, 0.5),
            't3': ('t', .125, 0.5),
            't4': ('t', 0.375, 0.5)}


class Uchip5555(Chip):
    """Chip of size 5 5 5 5"""

    pins = {'l1': ('l', -0.5, 0.4),
            'l2': ('l', -0.5, 0.2),
            'l3': ('l', -0.5, 0.0),
            'l4': ('l', -0.5, -0.2),
            'l5': ('l', -0.5, -0.4),
            'b1': ('b', -0.4, -0.5),
            'b2': ('b', -0.2, -0.5),
            'b3': ('b', 0.0, -0.5),
            'b4': ('b', 0.2, -0.5),
            'b5': ('b', 0.4, -0.5),
            'r5': ('r', 0.5, -0.4),
            'r4': ('r', 0.5, -0.2),
            'r3': ('r', 0.5, 0.0),
            'r2': ('r', 0.5, 0.2),
            'r1': ('r', 0.5, 0.4),
            't1': ('t', -0.4, 0.5),
            't2': ('t', -0.2, 0.5),
            't3': ('t', 0, 0.5),
            't4': ('t', 0.2, 0.5),
            't5': ('t', 0.4, 0.5)}


class Uchip6666(Chip):
    """Chip of size 6 6 6 6"""

    pins = {'l1': ('l', -0.5, 0.375),
            'l2': ('l', -0.5, 0.225),
            'l3': ('l', -0.5, 0.075),
            'l4': ('l', -0.5, -0.075),
            'l5': ('l', -0.5, -0.225),
            'l6': ('l', -0.5, -0.375),
            'b1': ('b', -0.375, -0.5),
            'b2': ('b', -0.225, -0.5),
            'b3': ('b', -0.075, -0.5),
            'b4': ('b', 0.075, -0.5),
            'b5': ('b', 0.225, -0.5),
            'b6': ('b', 0.375, -0.5),
            'r6': ('r', 0.5, -0.375),
            'r5': ('r', 0.5, -0.225),
            'r4': ('r', 0.5, -0.075),
            'r3': ('r', 0.5, 0.075),
            'r2': ('r', 0.5, 0.225),
            'r1': ('r', 0.5, 0.375),
            't1': ('t', -0.375, 0.5),
            't2': ('t', -0.225, 0.5),
            't3': ('t', -0.075, 0.5),
            't4': ('t', 0.075, 0.5),
            't5': ('t', 0.225, 0.5),
            't6': ('t', 0.357, 0.5)}


class Uchip7777(Chip):
    """Chip of size 7 7 7 7"""

    pins = {'l1': ('l', -0.5, 0.375),
            'l2': ('l', -0.5, 0.25),
            'l3': ('l', -0.5, 0.125),
            'l4': ('l', -0.5, 0.0),
            'l5': ('l', -0.5, -0.125),
            'l6': ('l', -0.5, -0.25),
            'l7': ('l', -0.5, -0.375),
            'b1': ('b', -0.375, -0.5),
            'b2': ('b', -0.25, -0.5),
            'b3': ('b', -0.125, -0.5),
            'b4': ('b', 0.0, -0.5),
            'b5': ('b', 0.125, -0.5),
            'b6': ('b', 0.25, -0.5),
            'b7': ('b', 0.375, -0.5),
            'r7': ('r', 0.5, -0.375),
            'r6': ('r', 0.5, -0.25),
            'r5': ('r', 0.5, -0.125),
            'r4': ('r', 0.5, 0.0),
            'r3': ('r', 0.5, 0.125),
            'r2': ('r', 0.5, 0.25),
            'r1': ('r', 0.5, 0.375),
            't1': ('t', -0.375, 0.5),
            't2': ('t', -0.25, 0.5),
            't3': ('t', -0.125, 0.5),
            't4': ('t', 0.0, 0.5),
            't5': ('t', 0.125, 0.5),
            't6': ('t', 0.25, 0.5),
            't7': ('t', 0.375, 0.5)}


class Uchip8181(Chip):
    """Chip of size 8 1 8 1"""

    pins = {'l1': ('l', -0.5, 0.4375),
            'l2': ('l', -0.5, 0.3125),
            'l3': ('l', -0.5, 0.1875),
            'l4': ('l', -0.5, 0.0625),
            'l5': ('l', -0.5, -0.0625),
            'l6': ('l', -0.5, -0.1875),
            'l7': ('l', -0.5, -0.3125),
            'l8': ('l', -0.5, -0.4375),
            'b1': ('b', 0.0, -0.5),
            'r8': ('r', 0.5, -0.4375),
            'r7': ('r', 0.5, -0.3125),
            'r6': ('r', 0.5, -0.1875),
            'r5': ('r', 0.5, -0.0625),
            'r4': ('r', 0.5, 0.0625),
            'r3': ('r', 0.5, 0.1875),
            'r2': ('r', 0.5, 0.3125),
            'r1': ('r', 0.5, 0.4375),
            't1': ('t', 0.0, 0.5)}
    aliases = {'vss': 'b1', 'vdd': 't1'}


class Uchip8888(Chip):
    """Chip of size 8 8 8 8"""

    pins = {'l1': ('l', -0.5, 0.4375),
            'l2': ('l', -0.5, 0.3125),
            'l3': ('l', -0.5, 0.1875),
            'l4': ('l', -0.5, 0.0625),
            'l5': ('l', -0.5, -0.0625),
            'l6': ('l', -0.5, -0.1875),
            'l7': ('l', -0.5, -0.3125),
            'l8': ('l', -0.5, -0.4375),
            'b1': ('b', -0.4375, -0.5),
            'b2': ('b', -0.3125, -0.5),
            'b3': ('b', -0.1875, -0.5),
            'b4': ('b', -0.0625, -0.5),
            'b5': ('b', 0.0625, -0.5),
            'b6': ('b', 0.1875, -0.5),
            'b7': ('b', 0.3125, -0.5),
            'b8': ('b', 0.4375, -0.5),
            'r8': ('r', 0.5, -0.4375),
            'r7': ('r', 0.5, -0.3125),
            'r6': ('r', 0.5, -0.1875),
            'r5': ('r', 0.5, -0.0625),
            'r4': ('r', 0.5, 0.0625),
            'r3': ('r', 0.5, 0.1875),
            'r2': ('r', 0.5, 0.3125),
            'r1': ('r', 0.5, 0.4375),
            't1': ('t', -0.4375, 0.5),
            't2': ('t', -0.3125, 0.5),
            't3': ('t', -0.1875, 0.5),
            't4': ('t', -0.0625, 0.5),
            't5': ('t', 0.0625, 0.5),
            't6': ('t', 0.1875, 0.5),
            't7': ('t', 0.3125, 0.5),
            't8': ('t', 0.4375, 0.5)}


class Mux(Chip):
    """Multiplexer"""

    @property
    def path(self):
        return ((-0.25, -0.5), (-0.25, 0.5), (0.25, 0.25), (0.25, -0.25), (-0.25, -0.5))


class Umux21(Mux):
    """Multiplexer 2 to 1"""

    pins = {'l1': ('l', -0.25, 0.25),
            'l2': ('l', -0.25, -0.25),
            'b': ('b', 0, -0.375),
            't': ('t', 0, 0.375),
            'r': ('r', 0.25, 0.0)}


class Umux41(Mux):
    """Multiplexer 4 to 1"""

    pins = {'l1': ('l', -0.25, 0.375),
            'l2': ('l', -0.25, 0.125),
            'l3': ('l', -0.25, -0.125),
            'l4': ('l', -0.25, -0.375),
            'b1': ('b', -0.125, -0.4375),
            'b2': ('b', 0.125, -0.3125),
            't1': ('t', -0.125, 0.4375),
            't2': ('t', 0.125, 0.3125),
            'r': ('r', 0.25, 0)}


class Umux42(Mux):
    """Multiplexer 4 to 2"""

    pins = {'l1': ('l', -0.25, 0.375),
            'l2': ('l', -0.25, 0.125),
            'l3': ('l', -0.25, -0.125),
            'l4': ('l', -0.25, -0.375),
            'b1': ('b', -0.125, -0.4375),
            'b2': ('b', 0.125, -0.3125),
            't1': ('t', -0.125, 0.4375),
            't2': ('t', 0.125, 0.3125),
            'r1': ('r', 0.25, 0.125),
            'r2': ('r', 0.25, -0.125)}


class Uadc(Chip):
    """ADC"""

    pins = {'in': ('l', -0.5, 0),
            'in+': ('l', -0.4375, 0.125),
            'in-': ('l', -0.4375, -0.125),
            'vref-': ('l', -0.375, -0.25),
            'vref+': ('l', -0.375, 0.25),
            'avss': ('b', -0.1, -0.5),
            'vss': ('b', 0.1, -0.5),
            'dvss': ('b', 0.3, -0.5),
            'clk': ('r', 0.5, -0.25),
            'data': ('r', 0.5, 0),
            'fs': ('r', 0.5, 0.25),
            'dvdd': ('t', 0.3, 0.5),
            'vdd': ('t', 0.1, 0.5),
            'avdd': ('t', -0.1, 0.5)}

    pinlabels = {'vref-': 'VREF-', 'vref+': 'VREF+',
                 'vss': 'VSS', 'vdd': 'VDD',
                 'dvss': 'DVSS', 'dvdd': 'DVDD',
                 'avss': 'AVSS', 'avdd': 'AVDD',
                 'clk': '>', 'data': 'DATA', 'fs': 'FS'}

    @property
    def path(self):
        return ((-0.5, 0.0), (-0.25, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.25, 0.5))


class Uregulator(Chip):
    """Voltage regulator"""

    default_aspect = 4.0 / 3.0
    pins = {'in': ('l', -0.5, 0),
            'en': ('b', -0.25, -0.5),
            'gnd': ('b', 0, -0.5),
            'out': ('r', 0.5, 0)}

    pinlabels = {'en': 'E', 'gnd': 'GND'}


class Udac(Chip):
    """DAC"""

    pins = {'out': ('r', 0.5, 0),
            'out+': ('r', 0.4375, 0.125),
            'out-': ('r', 0.4375, -0.125),
            'vref-': ('r', 0.375, -0.25),
            'vref+': ('r', 0.375, 0.25),
            'avss': ('b', 0.1, -0.5),
            'vss': ('b', -0.1, -0.5),
            'dvss': ('b', -0.3, -0.5),
            'clk': ('l', -0.5, -0.25),
            'data': ('l', -0.5, 0),
            'fs': ('l', -0.5, 0.25),
            'dvdd': ('t', -0.3, 0.5),
            'vdd': ('t', -0.1, 0.5),
            'avdd': ('t', 0.1, 0.5)}

    pinlabels = {'vref-': 'VREF-', 'vref+': 'VREF+',
                 'vss': 'VSS', 'vdd': 'VDD',
                 'dvss': 'DVSS', 'dvdd': 'DVDD',
                 'avss': 'AVSS', 'avdd': 'AVDD',
                 'clk': '>', 'data': 'DATA', 'fs': 'FS'}

    @property
    def path(self):
        return ((-0.5, -0.5), (0.25, -0.5), (0.5, 0), (0.25, 0.5), (-0.5, 0.5))


class Udiffamp(Chip):
    """Differential amplifier.  It is not automatically annotated with + and - symbols for inputs.
    This can be achieved using pinlabels={in+, in-}."""

    default_width = 1.0

    pins = {'in+': ('l', -0.5, 0.25),
            'in-': ('l', -0.5, -0.25),
            'vss': ('b', 0, -0.25),
            'out': ('r', 0.5, 0),
            'vdd': ('t', 0, 0.25)}

    pinlabels = {'in+': '$+$', 'in-': '$-$', 'vss': 'VSS', 'vdd': 'VDD'}

    @property
    def path(self):
        return ((-0.5, 0.5), (-0.5, -0.5), (0.5, 0))


class Ubuffer(Chip):
    """Buffer with power supplies"""

    default_width = 1.0

    pins = {'in': ('l', -0.5, 0),
            'vss': ('b', 0, -0.25),
            'out': ('r', 0.5, 0),
            'vdd': ('t', 0, 0.25),
            'vdd1': ('t', -0.25, 0.375),
            'vdd2': ('t', 0.25, 0.125),
            'en': ('b', -0.25, -0.375)}

    pinlabels = {'vss': 'VSS', 'vdd': 'VDD', 'en': 'E'}

    @property
    def path(self):
        return ((-0.5, 0.5), (0.5, 0), (-0.5, -0.5))


class Uinverter(Chip):
    """Inverter with power supplies"""

    default_width = 1.0

    pins = {'in': ('l', -0.5, 0),
            'vss': ('b', 0, -0.22),
            'out': ('r', 0.5, 0),
            'vdd': ('t', 0, 0.22),
            'en': ('b', -0.25, -0.37)}

    pinlabels = {'vss': 'VSS', 'vdd': 'VDD', 'en': 'E'}

    @property
    def path(self):
        w = 0.05
        return ((-0.5, 0.5), (0.5 - 2 * w, 0), (-0.5, -0.5))

    def draw(self, **kwargs):

        s = super(Uinverter, self).draw(**kwargs)

        # Append inverting circle.
        centre = self.node('mid')
        q = self.tf(centre.pos, ((0.45, 0)))
        s += r'  \draw[thick] (%s) node[ocirc, scale=%s, %s] {};''\n' % (
            q, 1.8 * self.size * self.scale, self.draw_args_str(**kwargs))
        return s


class Udiffdriver(Chip):
    """Differential driver with power supplies"""

    default_width = 1.0

    pins = {'in': ('l', -0.5, 0),
            'vss': ('b', -0.15, -0.3),
            'out+': ('r', -0.02, 0.25),
            'out-': ('r', 0.1, -0.25),
            'vdd': ('t', 0, 0.22),
            'en': ('b', -0.3, -0.4)}

    pinlabels = {'vss': 'VSS', 'vdd': 'VDD', 'en': 'E'}

    @property
    def path(self):
        w = 0.05
        return ((-0.5, 0.5), (0.5 - 2 * w, 0), (-0.5, -0.5))

    def draw(self, **kwargs):

        s = super(Udiffdriver, self).draw(**kwargs)

        # Append inverting circle.
        centre = self.node('mid')
        q = self.tf(centre.pos, ((0.05, -0.25)))
        s += r'  \draw[thick] (%s) node[ocirc, scale=%s, %s] {};''\n' % (
            q, 1.8 * self.size * self.scale, self.draw_args_str(**kwargs))
        return s


class Flipflop(Chip):

    default_width = 1.0


class Udff(Flipflop):
    """D flip-flop"""

    default_pins = ('d', 'clk', 'q')

    pins = {'d': ('l', -0.5, 0.25),
            'clk': ('l', -0.5, 0),
            'vss': ('b', 0.0, -0.5),
            '/q': ('r', 0.5, -0.25),
            'q': ('r', 0.5, 0.25),
            'vdd': ('t', 0.0, 0.5)}

    pinlabels = {'vss': 'VSS', 'vdd': 'VDD',
                 'd': 'D', 'q': 'Q', '/q': '$\overline{\mathrm{Q}}$', 'clk': '>'}


class Ujkff(Flipflop):
    """JK flip-flop"""

    default_pins = ('j', 'k', 'clk', 'q')

    pins = {'j': ('l', -0.5, 0.25),
            'clk': ('l', -0.5, 0),
            'k': ('l', -0.5, -0.25),
            'vss': ('b', 0.0, -0.5),
            '/q': ('r', 0.5, -0.25),
            'q': ('r', 0.5, 0.25),
            'vdd': ('t', 0.0, 0.5)}

    pinlabels = {'vss': 'VSS', 'vdd': 'VDD',
                 'j': 'J', 'k': 'K',
                 'q': 'Q', '/q': '$\overline{\mathrm{Q}}$', 'clk': '>'}


class Urslatch(Flipflop):
    """RS latch"""

    default_pins = ('r', 's', 'q')

    pins = {'r': ('l', -0.5, 0.25),
            's': ('l', -0.5, -0.25),
            'vss': ('b', 0.0, -0.5),
            '/q': ('r', 0.5, -0.25),
            'q': ('r', 0.5, 0.25),
            'vdd': ('t', 0.0, 0.5)}

    pinlabels = {'vss': 'VSS', 'vdd': 'VDD',
                 'r': 'R', 's': 'S',
                 'q': 'Q', '/q': '$\overline{\mathrm{Q}}$'}


class Gate2(Chip):

    can_mirror = False
    can_scale = True
    can_invert = False

    pins = {'in1': ('l', -0.6, 0.15),
            'in2': ('l', -0.6, -0.15),
            'out': ('r', 0.55, 0)}

    def draw(self, **kwargs):

        centre = self.node('mid')
        q = self.tf(centre.pos, ((-0.3, 0.075)))
        args = self.draw_args(self.opts, **kwargs)
        args.append('anchor=in 1')
        s = self.draw_cptnode(q, cpt='ieeestd ' + self.kind + ' port',
                              args=args)
        s += self.draw_label(centre.s, **kwargs)
        return s


class Uand(Gate2):

    kind = 'and'


class Uor(Gate2):

    kind = 'or'


class Unand(Gate2):

    kind = 'nand'


class Unor(Gate2):

    kind = 'nor'


class Uxor(Gate2):

    kind = 'xor'


class Uxnor(Gate2):

    kind = 'xnor'


class Eopamp(Chip):
    """This is for an opamp created with the E netlist type as used for
    simulation."""

    can_scale = True
    can_mirror = True
    do_transpose = False
    default_width = 1.0

    # The Nm node is not used (ground).
    node_pinnames = ('out', '', 'in+', 'in-')

    normal_pins = {'out': ('rx', 1.25, 0.0),
                   'in+': ('lx', -1.25, 0.5),
                   'in-': ('lx', -1.25, -0.5),
                   'vdd': ('t', 0, 0.5),
                   'vdd2': ('t', -0.45, 0.755),
                   'vss2': ('b', -0.45, -0.755),
                   'vss': ('b', 0, -0.5),
                   'ref': ('b', 0.45, -0.245),
                   'r+': ('l', -0.85, 0.25),
                   'r-': ('l', -0.85, -0.25)}

    mirror_pins = {'out': ('rx', 1.25, 0.0),
                   'in-': ('lx', -1.25, 0.5),
                   'in+': ('lx', -1.25, -0.5),
                   'vdd': ('t', 0, 0.5),
                   'vdd2': ('t', -0.45, 0.755),
                   'vss2': ('b', -0.45, -0.755),
                   'vss': ('b', 0, -0.5),
                   'ref': ('b', 0.45, -0.245),
                   'r-': ('l', -0.85, 0.25),
                   'r+': ('l', -0.85, -0.25)}

    pinlabels = {'vdd': 'VDD', 'vss': 'VSS'}

    @property
    def pins(self):
        return self.mirror_pins if (self.mirrorinputs ^ self.mirror) else self.normal_pins

    def draw(self, **kwargs):

        if not self.check():
            return ''

        yscale = 2 * 0.95 * self.scale
        if not (self.mirror ^ self.mirrorinputs):
            yscale = -yscale

        centre = self.node('mid')

        # The circuitikz opamp has short wires on input and output.

        # Note, scale scales by area, xscale and yscale scale by length.
        s = r'  \draw (%s) node[op amp, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            centre.s,
            self.draw_args_str(**kwargs), 2 * self.scale * 0.95, yscale,
            -self.angle, self.s)
        if not self.nowires:
            s += r'  \draw (%s.out) |- (%s);''\n' % (self.s,
                                                     self.node('out').s)
            s += r'  \draw (%s.+) |- (%s);''\n' % (self.s, self.node('in+').s)
            s += r'  \draw (%s.-) |- (%s);''\n' % (self.s, self.node('in-').s)
        s += self.draw_label(centre.s, **kwargs)
        return s


class Efdopamp(Chip):
    """This is for a fully differential opamp created with the E netlist
    type.  See also Ufdopamp for a fully differential opamp created
    with the U netlist type."""

    can_scale = True
    can_mirror = True
    do_transpose = False
    default_width = 1.0

    node_pinnames = ('out+', 'out-', 'in+', 'in-', 'ocm')

    normal_pins = {'out+': ('r', 0.85, -0.5),
                   'out-': ('r', 0.85, 0.5),
                   'in+': ('l', -1.25, 0.5),
                   'ocm': ('l', -0.85, 0),
                   'in-': ('l', -1.25, -0.5),
                   'vdd': ('t', -0.25, 0.645),
                   'vss': ('b', -0.25, -0.645),
                   'r+': ('l', -0.85, 0.25),
                   'r-': ('l', -0.85, -0.25)}

    mirror_pins = {'out-': ('r', 0.85, -0.5),
                   'out+': ('r', 0.85, 0.5),
                   'in-': ('l', -1.25, 0.5),
                   'ocm': ('l', -0.85, 0),
                   'in+': ('l', -1.25, -0.5),
                   'vdd': ('t', -0.25, 0.645),
                   'vss': ('b', -0.25, -0.645),
                   'r-': ('l', -0.85, 0.25),
                   'r+': ('l', -0.85, -0.25)}

    pinlabels = {'vdd': 'VDD', 'vss': 'VSS'}

    @property
    def pins(self):
        return self.mirror_pins if (self.mirrorinputs ^ self.mirror) else self.normal_pins

    def draw(self, **kwargs):

        if not self.check():
            return ''

        centre = self.node('mid')

        yscale = 2 * 0.952 * self.scale
        if not (self.mirror ^ self.mirrorinputs):
            yscale = -yscale

        s = r'  \draw (%s) node[fd op amp, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            centre.s, self.draw_args_str(**kwargs), 2 * self.scale * 0.95, yscale,
            -self.angle, self.s)
        s += r'  \draw (%s.out +) |- (%s);''\n' % (self.s, self.node('out+').s)
        s += r'  \draw (%s.out -) |- (%s);''\n' % (self.s, self.node('out-').s)
        s += r'  \draw (%s.+) |- (%s);''\n' % (self.s, self.node('in+').s)
        s += r'  \draw (%s.-) |- (%s);''\n' % (self.s, self.node('in-').s)

        s += self.draw_label(centre.s, **kwargs)
        return s


class Einamp(Eopamp):
    """Instrumentation amplifier created with E netlist type.
       See also Uinamp for an instrumentation amplifier created
       with the U netlist type."""

    can_scale = True
    can_mirror = True
    do_transpose = False
    default_width = 1.0

    node_pinnames = ('out', 'ref', 'in+', 'in-', 'r+', 'r-')

    auxiliary = {'lin+': ('c', -0.375, 0.3),
                 'lin-': ('c', -0.375, -0.3)}
    auxiliary.update(Chip.auxiliary)


class Eamp(Chip):
    """Amplifier."""

    can_scale = True
    do_transpose = False
    default_width = 1.0

    # The Nm and Ncm nodes are not used (ground).
    node_pinnames = ('out', '', 'in', '')

    pins = {'out': ('rx', 1.25, 0.0),
            'in': ('lx', -1.25, 0.0),
            'vdd': ('t', 0, 0.5),
            'vdd2': ('t', -0.45, 0.755),
            'vss2': ('b', -0.45, -0.755),
            'vss': ('b', 0, -0.5),
            'ref': ('b', 0.45, -0.245),
            'r+': ('l', -0.85, 0.25),
            'r-': ('l', -0.85, -0.25)}

    pinlabels = {'vdd': 'VDD', 'vss': 'VSS'}

    def draw(self, **kwargs):

        if not self.check():
            return ''

        centre = self.node('mid')
        yscale = 2 * 0.952 * self.scale

        # The circuitikz buffer has short wires on input and output.

        # Note, scale scales by area, xscale and yscale scale by length.
        s = r'  \draw (%s) node[buffer, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            centre.s,
            self.draw_args_str(**kwargs), 2 * self.scale * 0.95, yscale,
            -self.angle, self.s)
        if not self.nowires:
            s += r'  \draw (%s.out) |- (%s);''\n' % (self.s,
                                                     self.node('out').s)
            s += r'  \draw (%s.in) |- (%s);''\n' % (self.s, self.node('in').s)
        s += self.draw_label(centre.s, **kwargs)
        return s


class Uopamp(Chip):
    """This is for an opamp created with the U netlist type.  It has no wires.
    See also Opamp for an opamp created with the E netlist type."""

    can_mirror = True

    auxiliary = {'lin+': ('c', -0.375, 0.25),
                 'lin-': ('c', -0.375, -0.25)}
    auxiliary.update(Chip.auxiliary)
    required_auxiliary = ('lin+', 'lin-', 'mid')

    normal_pins = {'out': ('r', 0.5, 0.0),
                   'in+': ('l', -0.5, 0.25),
                   'in-': ('l', -0.5, -0.25),
                   'vdd': ('t', 0, 0.25),
                   'vdd2': ('t', -0.225, 0.365),
                   'vss2': ('b', -0.225, -0.365),
                   'vss': ('b', 0, -0.25),
                   'ref': ('b', 0.225, -0.135),
                   'r+': ('l', -0.5, 0.125),
                   'r-': ('l', -0.5, -0.125)}

    mirror_pins = {'out': ('r', 0.5, 0.0),
                   'in-': ('l', -0.5, 0.25),
                   'in+': ('l', -0.5, -0.25),
                   'vdd': ('t', 0, 0.25),
                   'vdd2': ('t', -0.225, 0.365),
                   'vss2': ('b', -0.225, -0.365),
                   'vss': ('b', 0, -0.25),
                   'ref': ('b', 0.225, -0.135),
                   'r-': ('l', -0.5, 0.125),
                   'r+': ('l', -0.5, -0.125)}

    pinlabels = {'vdd': 'VDD', 'vss': 'VSS'}

    @property
    def pins(self):
        return self.mirror_pins if self.mirrorinputs else self.normal_pins

    @property
    def path(self):
        return ((-0.5, -0.5), (-0.5, 0.5), (0.5, 0))

    def draw(self, **kwargs):

        s = super(Uopamp, self).draw(**kwargs)

        if not self.nolabels:
            if self.mirrorinputs:
                s += self.annotate(self.node('lin+').s, '$-$', bold=True)
                s += self.annotate(self.node('lin-').s, '$+$', bold=True)
            else:
                s += self.annotate(self.node('lin+').s, '$+$', bold=True)
                s += self.annotate(self.node('lin-').s, '$-$', bold=True)
        return s


class Ufdopamp(Chip):
    """This is for a fully differential opamp created with the U netlist
    type.  It has no wires.  See also Efdopamp for a fully differential
    opamp created with the E netlist type.

    """
    can_mirror = True

    auxiliary = {'lin+': ('c', -0.375, 0.2),
                 'lin-': ('c', -0.375, -0.2),
                 'lout+': ('c', 0, -0.17),
                 'lout-': ('c', 0, 0.17)}
    auxiliary.update(Chip.auxiliary)
    required_auxiliary = ('lin+', 'lin-', 'lout+', 'lout-', 'mid')

    normal_pins = {'out-': ('r', 0.1, 0.2),
                   'out+': ('r', 0.1, -0.2),
                   'in+': ('l', -0.5, 0.2),
                   'in-': ('l', -0.5, -0.2),
                   'vdd': ('t', -0.1, 0.3),
                   'vss': ('b', -0.1, -0.3),
                   'ocm': ('l', -0.5, 0)}

    mirror_pins = {'out+': ('r', 0.1, +0.2),
                   'out-': ('r', 0.1, -0.2),
                   'in-': ('l', -0.5, 0.2),
                   'in+': ('l', -0.5, -0.2),
                   'vdd': ('t', -0.1, 0.3),
                   'vss': ('b', -0.1, -0.3),
                   'ocm': ('l', -0.5, 0)}

    pinlabels = {'vdd': 'VDD', 'vss': 'VSS'}

    @property
    def pins(self):
        return self.mirror_pins if self.mirrorinputs else self.normal_pins

    @property
    def path(self):
        return ((-0.5, -0.5), (-0.5, 0.5), (0.5, 0))

    def draw(self, **kwargs):

        s = super(Ufdopamp, self).draw(**kwargs)

        if not self.nolabels:
            if self.mirrorinputs:
                s += self.annotate(self.node('lin+').s, '$-$', bold=True)
                s += self.annotate(self.node('lin-').s, '$+$', bold=True)
                s += self.annotate(self.node('lout+').s, '$-$', bold=True)
                s += self.annotate(self.node('lout-').s, '$+$', bold=True)
            else:
                s += self.annotate(self.node('lin+').s, '$+$', bold=True)
                s += self.annotate(self.node('lin-').s, '$-$', bold=True)
                s += self.annotate(self.node('lout+').s, '$+$', bold=True)
                s += self.annotate(self.node('lout-').s, '$-$', bold=True)
        return s


class Uinamp(Uopamp):
    """Instrumentation amplifier created with U netlist type."""

    can_mirror = True

    auxiliary = {'lin+': ('c', -0.375, 0.3),
                 'lin-': ('c', -0.375, -0.3)}
    auxiliary.update(Chip.auxiliary)

    normal_pins = {'out': ('r', 0.5, 0.0),
                   'in+': ('l', -0.5, 0.3),
                   'in-': ('l', -0.5, -0.3),
                   'vdd': ('t', 0, 0.25),
                   'vdd2': ('t', -0.225, 0.365),
                   'vss2': ('b', -0.225, -0.365),
                   'vss': ('b', 0, -0.25),
                   'ref': ('b', 0.225, -0.135),
                   'r+': ('l', -0.5, 0.2),
                   'r-': ('l', -0.5, -0.2)}

    mirror_pins = {'out': ('r', 0.5, 0.0),
                   'in-': ('l', -0.5, 0.3),
                   'in+': ('l', -0.5, -0.3),
                   'vdd': ('t', 0, 0.25),
                   'vdd2': ('t', -0.225, 0.365),
                   'vss2': ('b', -0.225, -0.365),
                   'vss': ('b', 0, -0.25),
                   'ref': ('b', 0.225, -0.135),
                   'r-': ('l', -0.5, 0.2),
                   'r+': ('l', -0.5, -0.2)}


class Uisoamp(Ufdopamp):
    """Isolated amplifier created with U netlist type."""

    auxiliary = {'lin+': ('c', -0.4, 0.2),
                 'lin-': ('c', -0.4, -0.2),
                 'lout+': ('c', 0.1, 0.12),
                 'lout-': ('c', 0.1, -0.12)}
    auxiliary.update(Chip.auxiliary)

    pins = {'out-': ('r', 0.2, -0.15),
            'out+': ('r', 0.2, 0.15),
            'out': ('r', 0.5, 0),
            'in+': ('l', -0.5, 0.2),
            'in-': ('l', -0.5, -0.2),
            'in': ('l', -0.5, 0.0),
            'vdd1': ('t', -0.3, 0.4),
            'vss1': ('b', -0.3, -0.4),
            'vdd2': ('t', 0.0, 0.25),
            'vss2': ('b', 0.0, -0.25)}

    @property
    def path(self):
        return [((-0.5, -0.5), (-0.5, 0.5), (-0.2, 0.35), (-0.2, -0.35)),
                ((-0.1, -0.3), (-0.1, 0.3), (0.5, 0), (-0.1, -0.3))]



class FB(Bipole):
    """Ferrite bead"""

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.nodes
        centre = (n1.pos + n2.pos) * 0.5
        w = 0.125
        h = 0.2

        q1 = self.tf(centre, ((-0.5 * w, -h), (-0.5 * w, h),
                              (0.5 * w, h), (0.5 * w, -h)), -30)
        q = self.tf(centre, ((-0.53 * w, 0), (0.53 * w, 0), (-w, -2 * h)))

        s = self.draw_path(q1, closed=True, style='thick')
        s += self.draw_path((n1.s, q[0]))
        s += self.draw_path((q[1], n2.s))
        s += self.draw_label(q[2], **kwargs)
        return s


class CPE(Bipole):
    """Constant phase element"""

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.nodes
        centre = (n1.pos + n2.pos) * 0.5
        w = 0.125
        h = 0.2

        q1 = self.tf(centre, ((0, -h), (-w, 0), (0, h)), 0)
        q2 = self.tf(centre, ((w, -h), (0, 0), (w, h)), 0)
        q = self.tf(centre, ((-w, 0), (0, 0), (-w, -2 * h), (-w, 2 * h)))

        s = self.draw_path(q1, closed=False, style='thick')
        s += self.draw_path(q2, closed=False, style='thick')
        s += self.draw_path((n1.s, q[0]))
        s += self.draw_path((q[1], n2.s))
        # FOO
        s += self.draw_label(q[2], ('l', 'l_'), default=True, **kwargs)
        s += self.draw_label(q[3], ('l^', ), default=False, **kwargs)
        return s


class XX(Cpt):
    directive = True

    def draw(self, **kwargs):

        if self.string.startswith(';;'):
            return ' ' + self.string[2:] + '\n'
        return ''

    def process_opts_nodes(self):
        pass


classes = {}


def defcpt(name, base, docstring, cpt=None):

    if isinstance(base, str):
        base = classes[base]

    newclass = type(name, (base, ), {'__doc__': docstring})

    if cpt is not None:
        newclass.tikz_cpt = cpt
    classes[name] = newclass


def make(classname, parent, namespace, name, cpt_type, cpt_id,
         string, opts_string, node_names, *args):

    # Create instance of component object
    try:
        newclass = getattr(module, classname)
    except:
        newclass = classes[classname]

    cpt = newclass(parent, namespace, name, cpt_type, cpt_id,
                   string, opts_string, node_names, *args)
    # Add named attributes for the args?   Lname1, etc.

    return cpt


# Dynamically create classes.
defcpt('ADC', Bipole, 'ADC', 'adc')
defcpt('AM', Bipole, 'Ammeter', 'ammeter')

defcpt('DAC', Bipole, 'DAC', 'dac')
# TODO, rationalise D 1 2 led; with D 1 2; kind=led, etc.
defcpt('Dled', D, 'LED', 'leD')
defcpt('Dphoto', D, 'Photo diode', 'pD')
defcpt('Dschottky', D, 'Schottky diode', 'zD')
defcpt('Dtunnel', D, 'Tunnel diode', 'tD')
defcpt('Dzener', D, 'Zener diode', 'zD')

defcpt('E', VCS, 'VCVS', 'american controlled voltage source')
defcpt('F', VCS, 'CCCS', 'american controlled current source')
defcpt('G', CCS, 'VCCS', 'american controlled current source')
defcpt('H', CCS, 'CCVS', 'american controlled voltage source')

defcpt('FS', Bipole, 'Fuse', 'fuse')

defcpt('GY', Gyrator, 'Gyrator', 'gyrator')

defcpt('TVtriode', Triode, 'Triode', 'triode')

defcpt('sI', Bipole, 'Current source', 'I')
defcpt('Isin', I, 'Sinusoidal current source', 'sI')
defcpt('Idc', I, 'DC current source', 'I')
defcpt('Istep', I, 'Step current source', 'I')
defcpt('Iac', I, 'AC current source', 'sI')
defcpt('Inoise', I, 'Noise current source', 'sI')

defcpt('J', JFET, 'N JFET transistor', 'njfet')
defcpt('Jnjf', 'J', 'N JFET transistor', 'njfet')
defcpt('Jpjf', 'J', 'P JFET transistor', 'pjfet')

defcpt('M', MOSFET, 'N MOSFET transistor', 'nmos')
defcpt('Mnmos', 'M', 'N channel MOSFET transistor', 'nmos')
defcpt('Mpmos', 'M', 'P channel MOSFET transistor', 'pmos')
defcpt('MISC', Bipole, 'Generic circuitikz component', '')

defcpt('NR', Bipole, 'Noiseless resistor', 'R')

defcpt('O', Bipole, 'Open circuit', 'open')
defcpt('P', Bipole, 'Port', 'open')

defcpt('Q', BJT, 'NPN transistor', 'npn')
defcpt('Qpnp', 'Q', 'PNP transistor', 'pnp')
defcpt('Qnpn', 'Q', 'NPN transistor', 'npn')

defcpt('RV', Potentiometer, 'Potentiometer', 'pR')

defcpt('Sbox', Box, 'Box shape')
defcpt('Scircle', Circle, 'Circle shape')
defcpt('Sellipse', Ellipse, 'Ellipse shape')
defcpt('Striangle', Triangle, 'Triangle shape')

defcpt('SW', Bipole, 'Switch', 'closing switch')
defcpt('SWno', 'SW', 'Normally open switch', 'closing switch')
defcpt('SWnc', 'SW', 'Normally closed switch', 'opening switch')
defcpt('SWpush', 'SW', 'Pushbutton switch', 'push button')
defcpt('SWspdt', SPDT, 'SPDT switch', 'spdt')

defcpt('TF', Transformer, 'Transformer', 'ideal transformer')
defcpt('TFcore', Transformer, 'Transformer with core', 'transformer core')
defcpt('TFtapcore', TFtap, 'Tapped transformer with core', 'transformer core')
defcpt('TP', TwoPort, 'Two port', '')
defcpt('TPA', TwoPort, 'A-parameter two port', '')
defcpt('TPB', TwoPort, 'B-parameter two port', '')
defcpt('TPG', TwoPort, 'G-parameter two port', '')
defcpt('TPH', TwoPort, 'H-parameter two port', '')
defcpt('TPY', TwoPort, 'Y-parameter two port', '')
defcpt('TPZ', TwoPort, 'Z-parameter two port', '')

defcpt('Ubox', Box2, 'Box')
defcpt('Ucircle', Circle2, 'Circle')
defcpt('Ubox4', Box4, 'Box')
defcpt('Ubox12', Box12, 'Box')
defcpt('Ucircle4', Circle4, 'Circle')
defcpt('sV', Bipole, 'Voltage source', 'V')
defcpt('Vsin', V, 'Sinusoidal voltage source', 'sV')
defcpt('Vdc', V, 'DC voltage source', 'V')
defcpt('Vstep', V, 'Step voltage source', 'V')
defcpt('Vac', V, 'AC voltage source', 'sV')
defcpt('Vnoise', V, 'Noise voltage source', 'sV')

defcpt('VCVS', VCS, 'VCVS', 'american controlled voltage source')
defcpt('CCCS', VCS, 'CCCS', 'american controlled current source')
defcpt('VCCS', CCS, 'VCCS', 'american controlled current source')
defcpt('CCVS', CCS, 'CCVS', 'american controlled voltage source')

defcpt('TLlossless', TL, 'Lossless transmission line', '')

defcpt('VM', Bipole, 'Voltmeter', 'voltmeter')

defcpt('W', Wire, 'Wire', 'short')

defcpt('XT', Bipole, 'Crystal', 'piezoelectric')

defcpt('m', Bipole, 'Mass', 'mass')
defcpt('k', Bipole, 'Spring', 'spring')
defcpt('r', Bipole, 'Damper', 'damper')

# Perhaps AM for ammeter, VM for voltmeter, VR for variable resistor?
# Currently, a variable resistor is supported with the variable
# option.
