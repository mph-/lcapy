from .chip import Chip

# The mirror_pins have the inputs (and outputs for FDA) flipped, not
# the other pins


class Eopamp(Chip):
    """This is for an opamp created with the E netlist type as used for
    simulation."""

    can_scale = True
    can_mirror = True
    do_transpose = False
    default_width = 1.0

    # The Nm node is electrically connected but is not drawn.
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
