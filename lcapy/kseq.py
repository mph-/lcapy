"""This module provides the DiscreteFourierDomainSequence class to
represent discrete-Fourier domain sequences.

Copyright 2021 Michael Hayes, UCECE

"""

from .domains import DiscreteFourierDomain
from .sequence import Sequence
from .sym import ksym

__all__ = ('kseq', )


class DiscreteFourierDomainSequence(DiscreteFourierDomain, Sequence):
    """Discrete-Fourier domain sequence."""

    var = ksym
    domain = 'discrete fourier sequence'

    def IDFT(self):
        """Calculate IDFT and return as sequence."""

        from sympy import exp
        from .sym import j, pi
        from .nexpr import n

        results = []
        vals = self.vals
        N = len(vals)
        for ni in range(N):
            result = 0
            for ki in range(N):
                result += vals[ki] * exp(2 * j * pi * ni * self.n[ki] / N)
            result = result.change(result, domain='discrete time')
            results.append(result / N)

        return self.change(results, domain='discrete time sequence')


def kseq(arg, ni=None, origin=None):
    """Create a discrete-Fourier domain Sequence from a tuple, list, ndarray, or str.

    >>> a = kseq((1, 2, 3))

    The sequence indices are specified with the optional `ni` argument.
    For example:

    >>> a = kseq((1, 2, 3, 4), (-1, 0, 1, 2))

    If the `ni` argument is not specified, the sequence indices
    are enumerated from 0.

    With a string argument, an underscore indicates the zero sequence
    index:

    >>> a = kseq('{1, _2, 3, 4}')

    The sequence indices can be found using the `n` attribute.
    This returns a list.

    >>> a = kseq('{1, _2, 3, 4}').n
    [-1, 0, 1, 2]
    """

    return DiscreteFourierDomainSequence(arg, ni, origin)


from .expressionclasses import expressionclasses  # nopep8

expressionclasses.register('discrete fourier sequence',
                           DiscreteFourierDomainSequence)
