from .shape import Shape


class Cable(Shape):
    """Cable"""

    default_aspect = 4
    a = 0.3
    pins = {'ignd': ('l', -0.5, -0.5),
            'ognd': ('r', 0.46, -0.5),
            't': ('c', 0, 0.5),
            'b': ('c', 0, -0.5)}

    auxiliary = {'in+': ('l', -0.5, a),
                 'in': ('l', -0.5, 0),
                 'in-': ('l', -0.5, -a),
                 'out+': ('r', 0.5, a),
                 'out': ('r', 0.5, 0),
                 'out-': ('r', 0.5, -a)}
    auxiliary.update(Shape.auxiliary)

    def draw(self, **kwargs):

        if not self.check():
            return ''

        centre = self.centre
        # Tikz uses height for length
        length = self.width
        width = self.height

        kind = self.kind
        if kind is None:
            kind = 'coax'

        if kind not in ('coax', 'twinax', 'twistedpair', 'shieldedtwistedpair',
                        'tline'):
            raise ValueError('Unknown cable kind %s' % kind)

        s = ''

        if kind in ('coax', 'twinax', 'shieldedtwistedpair', 'tline'):
            xscale = -1.025

            q = self.tf(centre, ((0.0125, 0)))

            s += r'  \draw[%s] (%s) node[cylinder, draw, rotate=%s, minimum width=%scm, minimum height=%scm, xscale=%s] {};''\n' % (
                self.draw_args_str(**kwargs), q, self.angle, width, length, xscale)

        if kind == 'tline':
            s += self.draw_label(centre, **kwargs)
        else:
            q = self.tf(centre, ((0, 0.9)))
            s += self.draw_label(q, **kwargs)

        if self.kind in ('twistedpair', 'shieldedtwistedpair'):
            # Needs to be even...
            twists = 4
            sections = twists * 4
            deltax = length / sections
            startx = centre.x - length / 2
            w = self.a * width

            x = startx
            y = self.node('mid').pos.y
            s += r'  \draw (%.2f,%.2f)' % (x, y + w)
            # 0 -1 0 1
            for m in range(sections):
                x += deltax
                a = [0, -1, 0, 1][m % 4] * w
                s += r' %s (%.2f,%.2f)' % (('cos', 'sin')[m % 2], x, y + a)
            s += ';\n'

            x = startx
            s += r'  \draw (%.2f,%.2f)' % (x, y - w)
            # 0 1 0 -1
            for m in range(sections):
                x += deltax
                a = [0, -1, 0, 1][m % 4] * w
                s += r' %s (%.2f,%.2f)' % (('cos', 'sin')[m % 2], x, y - a)
            s += ';\n'

        elif self.kind == 'coax':
            s += r'  \draw[-] (%s) to (%s);''\n' % (self.node('in').s,
                                                    self.node('out').s)
        elif self.kind == 'twinax':
            s += r'  \draw[-] (%s) to (%s);''\n' % (self.node('in-').s,
                                                    self.node('out-').s)
            s += r'  \draw[-] (%s) to (%s);''\n' % (self.node('in+').s,
                                                    self.node('out+').s)
        return s
