from numpy import ndarray
from .expr import expr
from .sequence import Sequence
from .nexpr import n

def seq(arg):

    if isinstance(arg, (tuple, list, ndarray)):
        nv = list(range(len(arg)))
        return Sequence(arg, nv, var=n)        

    s = arg
    if s.startswith('{'):
        if not s.endswith('}'):
            raise ValueError('Mismatched braces for %s' % s)
    
    parts = s.split(',')
    N = len(parts)

    vals = []
    m0 = None
    for m, item in enumerate(parts):
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
