"""
This module converts values into engineering format.

Copyright 2021 Michael Hayes, UCECE
"""
from math import floor, log10


class EngFormatter(object):

    def __init__(self, trim=True, hundreds=False, num_digits=3):
        """If `hundreds` True format like 100 pF rather than 0.1 nF"""

        self.trim = trim
        self.hundreds = hundreds
        self.num_digits = num_digits

    def _fmt(self, valstr, unit='', prefix='', aslatex=True):

        if unit == '' and prefix == '':
            return valstr

        if not aslatex:
            return valstr + ' ' + prefix + unit

        s = valstr + '\\,'
        if prefix == 'u':
            if unit.startswith('$'):
                s += '\\mu' + unit[1:-1]
            else:
                s += '\\mu\\mathrm{' + unit + '}'
        else:
            if unit.startswith('$'):
                s += '\\mathrm{' + prefix + '}' + unit[1:-1]
            else:
                s += '\\mathrm{' + prefix + unit + '}'
        return s


    def _do(self, value, unit, aslatex):

        prefixes = ('f', 'p', 'n', 'u', 'm', '', 'k', 'M', 'G', 'T')

        value = value
        if value == 0:
            return self._fmt('0', unit, '', aslatex)

        m = log10(abs(value))

        if m < -1 or m >= 3.0:
            if not self.hundreds:
                m += 1
            n = int(floor(m / 3))
            k = int(floor(m)) - n * 3
        else:
            n = 0
            k = m - 1

        dp = self.num_digits - k

        idx = n + 5
        if idx < 0:
            idx = 0
            return self._fmt('%e' % value, unit, '', aslatex)
        elif idx >= len(prefixes):
            idx = len(prefixes) - 1
            return self._fmt('%e' % value, unit, '', aslatex)

        if dp < 0:
            dp = 0
        fmt = '%%.%df' % dp

        n = idx - 5
        value = value * 10**(-3 * n)

        valstr = fmt % value

        if self.trim:
            # Remove trailing zeroes after decimal point
            valstr = valstr.rstrip('0').rstrip('.')

        return self._fmt(valstr, unit, prefixes[idx], aslatex)

    def latex_math(self, value, unit):
        """This is equivalent to the `latex()` method but encloses in $ $."""

        return '$' + self.latex(value, unit) + '$'

    def latex(self, value, unit):
        """Make latex string assuming math mode."""

        return self._do(value, unit, aslatex=True)

    def str(self, value, unit):
        """Make string."""

        return self._do(value, unit, aslatex=False)
