import numpy as np

class Pos(object):

    def __init__(self, x, y=0):

        if isinstance(x, tuple):
            x, y = x
        elif isinstance(x, np.ndarray):
            x, y = x[0], x[1]

        self.x = x
        self.y = y

    def __mul__(self, scale):

        return Pos(self.x * scale, self.y * scale)

    def __add__(self, arg):

        if not isinstance(arg, Pos):
            arg = Pos(arg)

        return Pos(self.x + arg.x, self.y + arg.y)

    def __str__(self):

        xstr = ('%.2f' % self.x).rstrip('0').rstrip('.')
        ystr = ('%.2f' % self.y).rstrip('0').rstrip('.')

        return "%s,%s" % (xstr, ystr)

    def __repr__(self):

        return 'Pos(%s)' % self

    @property
    def xy(self):

        return np.array((self.x, self.y))
