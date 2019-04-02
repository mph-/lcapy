from __future__ import division
from .sym import j
from .functions import sin, cos, exp

__all__ = ('Phasor', 'Vphasor', 'Iphasor')

from .omegaexpr import omegaExpr


class Phasor(omegaExpr):

    # Could convert Vphasor + Vconst -> VSuper but that is not really
    # the scope for types such as Vphasor and Vconst.

    def __init__(self, val, **assumptions):

        assumptions['positive'] = True  # ????
        assumptions['ac'] = True        
        super (Phasor, self).__init__(val, **assumptions)

    def __compat_add__(self, x, op):

        cls = self.__class__
        xcls = x.__class__

        # Special case for zero.
        if isinstance(x, int) and x == 0:
            return cls, self, cls(x), self.assumptions

        if not isinstance(x, Phasor):
            raise TypeError('Incompatible arguments %s and %s for %s' %
                            (repr(self), repr(x), op))

        if self.omega != x.omega:
            raise ValueError('Cannot combine %s(%s, omega=%s)'
                             ' with %s(%s, omega=%s)' %
                             (cls.__name__, self, self.omega,
                              xcls.__name__, x, x.omega))
        return cls, self, x, self.assumptions

    def __compat_mul__(self, x, op):

        cls = self.__class__
        xcls = x.__class__

        # Perhaps check explicitly for int, float?
        if not isinstance(x, Expr):
            return cls, self, cls(x), self.assumptions

        if isinstance(x, omegaExpr):
            return cls, self, x, self.assumptions

        if not isinstance(x, Phasor):
            raise TypeError('Incompatible arguments %s and %s for %s' %
                            (repr(self), repr(x), op))

        if self.omega != x.omega:
            raise ValueError('Cannot combine %s(%s, omega=%s)'
                             ' with %s(%s, omega=%s)' %
                             (cls.__name__, self, self.omega,
                              xcls.__name__, x, x.omega))
        return cls, self, x, self.assumptions

    def time(self, **assumptions):
        """Convert to time domain representation."""
        from .symbols import t
        
        omega = self.omega
        if isinstance(omega, Expr):
            # TODO: Fix inconsistency.  Sometimes omega is a symbol.
            omega = omega.expr
            
        if self.is_complex:
            result = self.expr * exp(j * omega * t)
        else:
            result = self.real.expr * cos(omega * t) - self.imag.expr * sin(omega * t)

        if hasattr(self, '_laplace_conjugate_class'):
            return self._laplace_conjugate_class(result)
        return tExpr(result)

    def fourier(self, **assumptions):
        """Attempt Fourier transform."""

        # TODO, could optimise.
        return self.time().fourier()        

    def laplace(self, **assumptions):
        """Convert to Laplace domain representation."""

        return self.time().laplace()

    def phasor(self):
        """Convert to phasor representation."""
        return self.__class__(self, **self.assumptions)

    def rms(self):
        return {Vphasor: Vt, Iphasor : It}[self.__class__](0.5 * self)

    # def plot(self, fvector=None, **kwargs):

    #     if self.omega != omegasym:
    #         self.fourier.plot(fvector, **kwargs)
    #     return omegaExpr(self).plot(fvector, **kwargs)


class Vphasor(Phasor):

    def __init__(self, val, **assumptions):

        super(Vphasor, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Vt

    def cpt(self):
        from .oneport import Vac
        return Vac(self, 0, self.omega)            

    
class Iphasor(Phasor):

    def __init__(self, val, **assumptions):
        super(Iphasor, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = It

    def cpt(self):
        from .oneport import Iac
        return Iac(self, 0, self.omega)

from .texpr import It, Vt, tExpr
from .expr import Expr
from .phasor import Phasor
