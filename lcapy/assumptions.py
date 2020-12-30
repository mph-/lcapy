"""This module provides the Assumptions class.  It manages
assumptions for expressions needed for inverse Laplace transforms
and expression simplification.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from copy import copy
from .acdc import is_dc, is_ac, is_causal


class Assumptions(dict):
    """Manage expression assumptions.  These include 'ac', 'dc', 'causal'
    and 'unknown' to specify expression behavior for t < 0, and 'positive',
    'complex', etc., to help with expression simplification.

    'ac', 'dc', 'causal' and 'unknown' refer to the time domain
    expression.  They can be inferred for a time domain expression but
    can be overridden.  The inference is performed lazily.  In
    particular, it is performed when converting from the time domain
    to the Laplace or Phasor domains.  The latter is required to
    handle the case of conversions from the Phasor to the Laplace
    domain.  Conversions from the Fourier domain to the Laplace domain
    go via the time domain and so the inference is performed while in
    the time domain.

    """

    def __init__(self, *args, **kwargs):

        super (Assumptions, self).__init__(*args)

        for assumption, value in kwargs.items():
            self.set(assumption, value)
    
    def set(self, assumption, value):

        if assumption in ('dc', 'ac', 'causal', 'unknown'):
            if value:
                self.pop('dc', None)
                self.pop('ac', None)
                self.pop('causal', None)
                self.pop('unknown', None)                
                self[assumption] = value
            else:
                self.pop(assumption, None)
        else:
            self[assumption] = value

    def get(self, assumption):

        return self[assumption]

    def merge(self, **assumptions):

        new = self.copy()
        for assumption, value in assumptions.items():
            new.set(assumption, value)
        return new

    def copy(self):

        return copy(self)
    
    def sympy_assumptions(self):
        """Return dict of the SymPy assumptions such as complex, positive, etc."""

        assumptions = {}
        for assumption, value in self.items():
            if assumption not in ('nid', 'ac', 'dc', 'causal', 'unknown'):
                assumptions[assumption] = value
        return assumptions

    def infer_from_expr(self, expr):

        if expr.is_transform_domain:
            return
        
        var = expr.var
        if is_dc(expr, var):
            self.set('dc', True)
            return

        if is_ac(expr, var):
            self.set('ac', True)            
            return

        if is_causal(expr, var):
            self.set('causal', True)            
            return

        self.set('unknown', True)                    
