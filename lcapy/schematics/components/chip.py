from .shape import Shape


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