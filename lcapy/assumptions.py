"""This module provides the Assumptions class.  It manages
assumptions for expressions needed for inverse Laplace transforms
and expression simplification.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from copy import copy
from .acdc import is_dc, is_ac, is_causal


class Assumptions(dict):
    """Manage expression assumptions.  These include 'ac', 'dc', 'causal'
    and 'unknown' to specify expression behavior for t < 0, 'positive',
    'complex', etc., to help with expression simplification, 'nid' as a noise
    identifier, and `omega` for the phasor angular frequency.

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

    def get(self, assumption, default=False):
        
        return super(Assumptions, self).get(assumption, default)

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
            self.set('unknown', True)
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

    @property
    def has_unspecified(self):

        if 'ac' in self:
            return False
        if 'dc' in self:
            return False
        if 'causal' in self:
            return False
        if 'unknown' in self:
            return False        
        return True

    def infer_from_expr_and_merge(self, expr, **assumptions):
        """If the user has not specified ac, dc, or causal, infer them.
         Perhaps should check what the user says and issue warning?
        Then override with `assumptions`."""

        if self.has_unspecified:
            self.infer_from_expr(expr)
        return self.merge(**assumptions)

    @property
    def is_ac(self):

        return self.get('ac')

    @property
    def is_dc(self):

        return self.get('dc')

    @property
    def is_causal(self):

        return self.get('causal')

    @property
    def is_unknown(self):

        return self.get('unknown')    
