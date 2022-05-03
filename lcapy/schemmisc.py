"""This module provides miscellaneous support for schematic drawing.

Copyright 2014--2022 Michael Hayes, UCECE

"""


class Pos(object):

    def __init__(self, x, y=0):

        from numpy import ndarray

        if isinstance(x, tuple):
            x, y = x
        elif isinstance(x, ndarray):
            x, y = x[0], x[1]

        self.x = x
        self.y = y

    def __mul__(self, scale):

        return Pos(self.x * scale, self.y * scale)

    def __add__(self, arg):

        if not isinstance(arg, Pos):
            arg = Pos(arg)

        return Pos(self.x + arg.x, self.y + arg.y)

    def __sub__(self, arg):

        if not isinstance(arg, Pos):
            arg = Pos(arg)

        return Pos(self.x - arg.x, self.y - arg.y)

    def __str__(self):

        xstr = ('%.3f' % self.x).rstrip('0').rstrip('.')
        ystr = ('%.3f' % self.y).rstrip('0').rstrip('.')

        return "%s,%s" % (xstr, ystr)

    def __repr__(self):

        return 'Pos(%s)' % self

    @property
    def xy(self):

        from numpy import array

        return array((self.x, self.y))


class Steps(object):

    def __init__(self, string, pos1, pos2):

        self.string = string

        self.numv = 0
        self.numh = 0
        self.steps = []
        self.pos1 = pos1
        self.pos2 = pos2
        self.pos = pos1

        self.dx = (pos2.x - pos1.x)
        self.dy = (pos2.y - pos1.y)

        if string == '':
            if abs(self.dx) > abs(self.dy):
                string = '-|'
            else:
                string = '|-'

        s = string
        m = 0
        while m < len(s):
            if s[m].isdigit():
                n = m
                while s[n].isdigit():
                    n += 1
                num = int(s[m:n])
                m = n
            else:
                num = 1
            if s[m] == '|':
                self.numv += num
                self.steps.append((s[m], num))
            elif s[m] == '-':
                self.numh += num
                self.steps.append((s[m], num))
            else:
                raise ValueError('Expecting - or | in %s' % string)
            m += 1

        self.dh = self.dx / self.numh
        self.dv = self.dy / self.numv

    def __iter__(self):
        self.n = 0
        self.pos = self.pos1
        return self

    def __next__(self):
        if self.n >= len(self.steps):
            raise StopIteration

        c, n = self.steps[self.n]
        if c == '-':
            self.pos += Pos(self.dh * n, 0)
        elif c == '|':
            self.pos += Pos(0, self.dv * n)

        self.n += 1
        return self.pos
