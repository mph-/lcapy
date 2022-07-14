"""This module provides the ZDomainSequence class to represent
z-domain sequences.

Copyright 2021 Michael Hayes, UCECE

"""

from .domains import ZDomain
from .sequence import Sequence
from .sym import zsym

__all__ = ('zseq', )


class ZDomainSequence(ZDomain, Sequence):
    """z-domain sequence."""

    var = zsym
    domain = 'Z sequence'

    def IZT(self):
        """Calculate inverse z-transform and return as sequence."""

        from .zexpr import z

        results = []
        vals = self.vals
        N = len(vals)
        for ni in range(N):
            result = vals[ni] * z**ni
            result = result.change(result, domain='discrete time')
            results.append(vals[ni] * z**ni)

        return self.change(results, domain='discrete time sequence')


def zseq(arg, ni=None, origin=None):
    """Create a discrete-time domain Sequence from a tuple, list, ndarray, or str.

    >>> a = zseq((1, 2, 3))

    The sequence indices are specified with the optional `ni` argument.
    For example:

    >>> a = zseq((1, 2, 3, 4), (-1, 0, 1, 2))

    If the `ni` argument is not specified, the sequence indices
    are enumerated from 0.

    With a string argument, an underscore indicates the zero sequence
    index:

    >>> a = zseq('{1, _2, 3, 4}')

    The sequence indices can be found using the `n` attribute.
    This returns a list.

    >>> a = zseq('{1, _2, 3, 4}').n
    [-1, 0, 1, 2]
    """

    return ZDomainSequence(arg, ni, origin)


from .expressionclasses import expressionclasses  # nopep8

expressionclasses.register('Z sequence', ZDomainSequence)
