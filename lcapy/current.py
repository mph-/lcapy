"""This module provides current support.

Copyright 2020 Michael Hayes, UCECE

"""
from .expr import expr
from .units import u as uu
from .classmap import domain_kind_to_symbol, domain_kind_quantity_to_class

def Iname(name, kind, cache=False):

    # Not caching is a hack to avoid conflicts of Vn1 with Vn1(s) etc.
    # when using subnetlists.  The alternative is a proper context
    # switch.  This would require every method to set the context.

    if kind == 't':
        name = name.lower()
    
    undef = domain_kind_to_symbol(kind, name)
    cls = domain_kind_quantity_to_class(kind, 'current')
    return cls(undef, cache=cache)


def Itype(kind):

    return domain_kind_quantity_to_class(kind, 'current')        


def current(arg, **assumptions):

    expr1 = expr(arg, **assumptions)

    try:
        expr1 = expr1.as_current()
    except:        
        raise ValueError('Cannot represent %s(%s) as current' % (expr1.__class__.__name__, expr1))

    return expr1.apply_unit(uu.amperes)
