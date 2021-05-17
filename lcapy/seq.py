"""This module creates sequences.

Copyright 2020-2021 Michael Hayes, UCECE

"""

from numpy import ndarray
from .expr import expr
from .sequence import Sequence
from .nexpr import n
from .utils import isiterable

def seq(arg, ni=None):
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
    """

    if not isiterable(arg):
        arg = (arg, )
    
    if isinstance(arg, (tuple, list, ndarray)):
        return Sequence(arg, ni, var=n)        

    if not isinstance(arg, str):
        raise ValueError('Argument not scalar, string, tuple, list, or ndarray.')
    
    s = arg
    if s.startswith('{'):
        if not s.endswith('}'):
            raise ValueError('Mismatched braces for %s' % s)
        s = s[1:-1]
    
    parts = s.split(',')
    N = len(parts)

    vals = []
    m0 = None
    for m, item in enumerate(parts):
        item = item.strip()
        if item.startswith('_'):
            if m0 is not None:
                raise ValueError('Cannot have multiple zero index indicators')
            m0 = m
            item = item[1:]

        val = expr(item)
        vals.append(val)

    if m0 is None:
        m0 = 0
        
    nv = list(range(-m0, N - m0))
    return Sequence(vals, nv, var=n)
