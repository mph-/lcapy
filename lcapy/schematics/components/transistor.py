from warnings import warn
from .fixedcpt import FixedCpt


class Transistor(FixedCpt):
    """Transistor"""

    can_mirror = True
    can_scale = True
    can_invert = True

    @property
    def pins(self):

        if (self.kind is not None
            and (self.kind.startswith('pigfet')
                 or self.kind.startswith('nigfet'))):
            xpins = [[self.normal_pins2, self.invert_pins2],
                     [self.mirror_pins2, self.mirror_invert_pins2]]
        else:
            xpins = [[self.normal_pins, self.invert_pins],
                     [self.mirror_pins, self.mirror_invert_pins]]
        if (self.classname in ('Qpnp', 'Mpmos', 'Jpjf')
            or self.kind in ('pmos', 'pmosd', 'pfetd', 'pfet',
                             'pfetd-bodydiode', 'pfet-bodydiode',
                             'pigfetd', 'pigfete', 'pigfetebulk')):
            pins = xpins[not self.mirror][self.invert]
        else:
            pins = xpins[self.mirror][self.invert]

        if self.size != 1 or self.scale != 1:
            if 'g' in pins:
                # Apply hack to draw gate in correct place when
                # size is not 1.  Only required if pos != 0.5.
                pins = pins.copy()
                gpin = pins['g']
                y = ((1 - self.scale) / 2 +
                     gpin[2] * self.scale + (self.size - 1) / 2) / self.size
                pins['g'] = (gpin[0], gpin[1], y)
        return pins

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3 = self.nodes
        centre = (n1.pos + n3.pos) * 0.5

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
