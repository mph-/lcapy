"""This modules provides the dExpr class to provide common methods for
the discrete-time representations.

Copyright 2020 Michael Hayes, UCECE

"""

import numpy as np
from .expr import Expr

class dExpr(Expr):

    """Superclass of discrete-time, discrete-frequency, and z-domain
    expressions."""
    

    def __call__(self, arg, **assumptions):
        """Transform domain or substitute arg for variable. 
        
        Substitution is performed if arg is a tuple, list, numpy
        array, or constant.  If arg is a tuple or list return a list.
        If arg is an numpy array, return numpy array.

        Domain transformation is performed if arg is a domain variable
        or an expression of a domain variable.

        See also evaluate.

        """

        if isinstance(arg, (tuple, list)):
            return [self._subs1(self.var, arg1) for arg1 in arg]

        if isinstance(arg, np.ndarray):
            return np.array([self._subs1(self.var, arg1) for arg1 in arg])

        from .discretetime import call        
        return call(self, arg, **assumptions)
