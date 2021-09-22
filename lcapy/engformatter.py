"""
This module converts values into engineering format.

Copyright 2021 Michael Hayes, UCECE
"""
from math import floor, log10

class EngFormatter(object):

    def __init__(self, trim=True, hundreds=False, sfmax=3):
        """If `hundreds` True format like 100 pF rather than 0.1 nF"""
        
        self.trim = trim
        self.hundreds = hundreds
        self.sfmax = sfmax
        
    def _do(self, value, unit, prefixes, space='', mbox_prefix='',
            mbox_suffix=''):

        value = value
        if value == 0:
            return '0' + space + mbox_prefix + unit + mbox_suffix

        m = log10(abs(value))

        if m < -1 or m >= 3.0:
            if not self.hundreds:
                m += 1
            n = int(floor(m / 3))
            k = int(floor(m)) - n * 3
        else:
            n = 0
            k = m - 1

        dp = self.sfmax - k

        idx = n + 5
        if idx < 0:
            idx = 0
            return '%e\,' % value + unit
        elif idx >= len(prefixes):
            idx = len(prefixes) - 1
            return '%e\,' % value + unit

        fmt = '%%.%df' % dp

        n = idx - 5
        value = value * 10**(-3 * n)

        string = fmt % value

        if self.trim:
            # Remove trailing zeroes after decimal point
            string = string.rstrip('0').rstrip('.')

        return string + space + mbox_prefix + prefixes[idx] + unit + mbox_suffix

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
    
