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

        xstr = ('%.3f' % self.x).rstrip('0').rstrip('.')
        ystr = ('%.3f' % self.y).rstrip('0').rstrip('.')

        return "%s,%s" % (xstr, ystr)

    def __repr__(self):

        return 'Pos(%s)' % self

    @property
    def xy(self):

        return np.array((self.x, self.y))


class Opts(dict):

    def __init__(self, arg=None):

        if arg is None:
            return

        if isinstance(arg, str):
            self.add(arg)
            return

        for key, val in arg.items():
            self[key] = val

    def add(self, string):

        def split(s):
            """Split a string by , except if in braces"""
            parts = []
            bracket_level = 0
            current = []
            for c in (s + ','):
                if c == ',' and bracket_level == 0:
                    parts.append(''.join(current))
                    current = []
                else:
                    if c == '{':
                        bracket_level += 1
                    elif c == '}':
                        bracket_level -= 1
                    current.append(c)
            if bracket_level != 0:
                raise ValueError('Mismatched braces for ' + s)
            return parts

        if string == '':
            return
        if string[0] == ';':
            self['append'] += string[1:] + '\n'
            return

        for part in split(string):
            part = part.strip()
            if part == '':
                continue

            fields = part.split('=')
            key = fields[0].strip()
            arg = '='.join(fields[1:]).strip() if len(fields) > 1 else ''
            if arg in ('true', 'True'):
                arg = True
            elif arg in ('false', 'False'):
                arg = False

            self[key] = arg

    def __str__(self):
        return self.format()

    @property
    def size(self):

        return float(self.get('size', 1))

    def format(self):

        def fmt(key, val):
            if val == '':
                return key
            return '%s=%s' % (key, val)

        return ', '.join([fmt(key, val) for key, val in self.items()])

    def copy(self):
        
        return self.__class__(super(Opts, self).copy())

    def strip(self, *args):

        stripped = Opts()
        for opt in args:
            if opt in self:
                stripped[opt] = self.pop(opt)        
        return stripped

    def strip_voltage_labels(self):

        return self.strip('v', 'vr', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<')

    def strip_current_labels(self):

        return self.strip('i', 'ir', 'i_', 'i^', 'i_>', 'i_<', 'i^>', 'i^<',
                          'i>_', 'i<_', 'i>^', 'i<^')

    def strip_labels(self):

        return self.strip('l', 'l^', 'l_')

    def strip_all_labels(self):

        self.strip_voltage_labels()
        self.strip_current_labels()
        self.strip_labels()
