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

    def _fmt(self, valstr, unit='', prefix='', space='', mbox_prefix='',
             mbox_suffix=''):

        if unit == '' and prefix == '':
            return valstr

        return valstr + space + mbox_prefix + prefix + unit + mbox_suffix

    def _do(self, value, unit, prefixes, space='', mbox_prefix='',
            mbox_suffix=''):

        value = value
        if value == 0:
            return self._fmt('0', unit, '', space, mbox_prefix, mbox_suffix)

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
            return self._fmt('%e' % value, unit, '', space, mbox_prefix, mbox_suffix)
        elif idx >= len(prefixes):
            idx = len(prefixes) - 1
            return self._fmt('%e' % value, unit, '', space, mbox_prefix, mbox_suffix)

        if dp < 0:
            dp = 0
        fmt = '%%.%df' % dp

        n = idx - 5
        value = value * 10**(-3 * n)

        valstr = fmt % value

        if self.trim:
            # Remove trailing zeroes after decimal point
            valstr = valstr.rstrip('0').rstrip('.')

        return self._fmt(valstr, unit, prefixes[idx], space,
                         mbox_prefix, mbox_suffix)

    def latex_math(self, value, unit):
        """Make latex math-mode string."""

        return '$' + self.latex(value, unit) + '$'

    def latex(self, value, unit):
        """Make latex string."""

        return self._do(value, unit,
                        ('f', 'p', 'n', '$\mu$', 'm', '', 'k', 'M', 'G', 'T'),
                        '\,', r'\mbox{', r'}')

    def str(self, value, unit):
        """Make string."""

        return self._do(value, unit,
                        ('f', 'p', 'n', 'u', 'm', '', 'k', 'M', 'G', 'T'),
                        ' ', '', '')
