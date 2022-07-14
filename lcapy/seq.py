"""This module creates sequences.

Copyright 2020--2022 Michael Hayes, UCECE

"""

from .expr import expr
from .sequence import Sequence
from .nexpr import n
from .utils import isiterable


def seq(arg, ni=None, origin=None, domain=None):
    """Create a Sequence from a tuple, list, ndarray, or str.

    >>> a = seq((1, 2, 3))

    The sequence indices are specified with the optional `ni` argument.
    For example:

    >>> a = seq((1, 2, 3, 4), (-1, 0, 1, 2))

    If the `ni` argument is not specified, the sequence indices
    are enumerated from 0.

    With a string argument, an underscore indicates the zero sequence
    index:

    >>> a = seq('{1, _2, 3, 4}')

    The sequence indices can be found using the `n` attribute.
    This returns a list.

    >>> a = seq('{1, _2, 3, 4}').n
    [-1, 0, 1, 2]

    This function creates a discrete-time domain sequence.  See also
    `nseq` to create a discrete-time domain sequence, `kseq` to create
    a discrete-Fourier domain sequence, and `zseq` to create a
    z-domain sequence.

    """

    if domain is None:
        # For backward-compatibility
        # TODO: If find a z in the arg, then create a zseq.
        return nseq(arg, ni, origin)
    elif domain is n:
        return nseq(arg, ni, origin)
    elif domain is k:
        return kseq(arg, ni, origin)
    elif domain is z:
        return zseq(arg, ni, origin)
    else:
        raise ValueError('Unknown domain %s' % domain)


from .nseq import nseq  # nopep8
from .kseq import kseq  # nopep8
from .zseq import zseq  # nopep8
from .symbols import n, k, z  # nopep8
