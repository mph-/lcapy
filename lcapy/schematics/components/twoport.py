from .shape import Shape


class TwoPort(Shape):
    """Two-port"""

    shape_scale = 4.0 / 3
    default_width = 1
    default_aspect = 1
    shape = 'rectangle'
    x = 0.5
    p = 0.375
    n = 0.75
    pins = {'w': ('l', -x, 0),
            'ssw': ('b', -p, -x),
            's': ('b', 0, -x),
            'sse': ('b', p, -x),
            'e': ('r', x, 0),
            'nne': ('t', p, x),
            'n': ('t', 0, x),
            'nnw': ('t', -p, x),
            'in+': ('lx', -n, p),
            'in-': ('lx', -n, -p),
            'out+': ('rx', n, p),
            'out-': ('rx', n, -p)}

    auxiliary = {'wnw': ('l', -x, p),
                 'wsw': ('l', -x, -p),
                 'ene': ('l', x, p),
                 'ese': ('l', x, -p)}
    auxiliary.update(Shape.auxiliary)
    required_auxiliary = ('wnw', 'wsw', 'ene', 'ese', 'mid')

    node_pinnames = ('out+', 'out-', 'in+', 'in-')

    def draw(self, **kwargs):

        self.shape = self.opts.get('shape', self.shape)
        scale = self.scale
        # TODO, change scaling if cloud puffs is specified
        # (the default is 10) or if aspect != 1.
        if self.shape == 'cloud':
            scale *= 1.07

        s = super(TwoPort, self).draw(**kwargs, scale=scale)
        s += r'  \draw (%s) -- (%s);''\n' % (self.node('in+').s,
                                             self.node('wnw').s)
        s += r'  \draw (%s) -- (%s);''\n' % (self.node('in-').s,
                                             self.node('wsw').s)
        s += r'  \draw (%s) -- (%s);''\n' % (self.node('ene').s,
                                             self.node('out+').s)
        s += r'  \draw (%s) -- (%s);''\n' % (self.node('ese').s,
                                             self.node('out-').s)

        return s
