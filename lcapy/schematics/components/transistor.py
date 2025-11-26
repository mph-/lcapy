from warnings import warn
from .fixedcpt import FixedCpt


class Transistor(FixedCpt):
    """Transistor"""

    can_mirror = True
    can_scale = True
    can_invert = True
    auxiliary = {'mid': ('c', 0.0, 0.0)}

    @property
    def tpins(self):

        if (self.kind is not None
            and (self.kind.startswith('pigfet')
                 or self.kind.startswith('nigfet'))):
            pins = self.pins2
        else:
            pins = self.pins
        if (self.classname in ('Qpnp', 'Mpmos', 'Jpjf')
            or self.kind in ('pmos', 'pmosd', 'pfetd', 'pfet',
                             'pfetd-bodydiode', 'pfet-bodydiode',
                             'pigfetd', 'pigfete', 'pigfetebulk')):
            tpins = self.tweak_pins(pins, True)
        else:
            tpins = self.tweak_pins(pins)

        return tpins

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3 = self.nodes[0:3]
        centre = (n1.pos + n3.pos) * 0.5
        # centre = self.node('mid').pos

        xscale = self.scale
        yscale = self.scale
        if self.mirror:
            yscale = -yscale
        if self.invert:
            xscale = -xscale

        cpt = self.tikz_cpt
        if self.kind is not None:
            # It may be better to have classes for each transistor kind
            # but this is a lot of extra work.

            if self.kind not in self.kinds:
                warn('Kind %s not in known kinds: %s' %
                     (self.kind, ', '.join(self.kinds)))
            # Handle bodydiode
            cpt = self.kind.replace('-', ', ')

        label = self.label(**kwargs)

        label = self.label_tweak(label, xscale, yscale, self.angle)

        s = r'  \draw (%s) node[%s, %s, xscale=%s, yscale=%s, rotate=%d] (%s) {%s};''\n' % (
            centre, cpt, self.cpt_args_str(**kwargs), xscale, yscale,
            self.angle, self.s, label)

        # Add additional wires.  These help to compensate for the
        # slight differences in sizes of the different transistors.
        if self.tikz_cpt in ('pnp', 'npn'):
            s += r'  \draw[%s] (%s.C) -- (%s) (%s.B) -- (%s) (%s.E) -- (%s);''\n' % (
                self.cpt_args_str(**kwargs), self.s, n1.s, self.s, n2.s, self.s, n3.s)
        else:
            s += r'  \draw (%s.D) -- (%s) (%s.G) -- (%s) (%s.S) -- (%s);''\n' % (
                self.s, n1.s, self.s, n2.s, self.s, n3.s)

        return s
