from __future__ import division
from .sym import j, omegasym, symsimplify
from .context import context
from .omegaexpr import omegaExpr
import sympy as sym

class noiseExpr(omegaExpr):
    """Frequency domain (one-sided) noise spectrum expression (amplitude
    spectral density).

    This characterises a zero-mean Gaussian noise process.

    When performing arithmetic on two noiseExpr expressions it is
    assumed that they are uncorrelated unless they have the same nid
    (noise indentifier).  If the nid is not specified, a new one is
    created.

    Uncorrelated noise expressions are added in quadrature (on a power
    basis).  Thus (Vn(3) + Vn(4)).expr = 5 since 5 = sqrt(3**2 + 4**2)

    Vn(3) != Vn(3) since they are different noise realisations albeit
    with the same properties.  However, Vn(3).expr == Vn(3).expr.
    Similarly, Vn(3, nid='n1') == Vn(3, nid='n1') since they have the
    same noise identifier and thus have the same realisation.

    Caution: The sum of two noise expressions generates a noise
    expression with a new nid.  This can lead to unexpected results
    since noise expressions with different nids are assumed to be
    uncorrelated.  For example, consider:
    a = Vn(3); b = Vn(4)
    a + b - b gives sqrt(41) and  a + b - a gives sqrt(34).

    This case is correctly handled by the Super class since each noise
    component is stored and considered separately.

    (Vsuper(a) + Vsuper(b) - Vsuper(b)).n gives 3 as expected.

    """
    one_sided = True

    def _new_nid(self):
        context.nid += 1
        return 'n%d' % context.nid

    def __init__(self, val, **assumptions):
        if 'nid' not in assumptions or assumptions['nid'] is None:
            if val == 0:
                assumptions['nid'] = 'n0'
            else:
                assumptions['nid'] = self._new_nid()
        super(noiseExpr, self).__init__(val, **assumptions)

    @property
    def nid(self):
        return self.assumptions['nid']
        
    def __add__(self, x):
        """Add noise spectra (on power basis if uncorrelated)."""

        if not isinstance(x, noiseExpr):
            raise ValueError('Cannot add %s and %s' % (self, x))
        
        if x == 0:
            return self.__class__(self, nid=self.nid)
        
        if self.nid == x.nid:
            return self.__class__(self.expr + x.expr, nid=self.nid)
        
        value1 = self.expr
        value2 = x.expr
        value1sq = symsimplify(value1 * sym.conjugate(value1))
        value2sq = symsimplify(value2 * sym.conjugate(value2))                  
        result = symsimplify(sym.sqrt((value1sq + value2sq)))
        return self.__class__(result)

    def __radd__(self, x):
        raise ValueError('Cannot add %s and %s' % (self, x))        

    def __sub__(self, x):
        if not isinstance(x, noiseExpr):
            raise ValueError('Cannot subtract %s and %s' % (self, x))

        if x == 0:
            return self.__class__(self, nid=self.nid)        

        if self.nid == x.nid:
            return self.__class__(self.expr - x.expr, nid=self.nid)        
        return self + x

    def __rsub__(self, x):
        raise ValueError('Cannot subtract %s and %s' % (self, x))                

    def __mul__(self, x):
        if isinstance(x, noiseExpr) and self.nid != x.nid:
            raise ValueError('Cannot multiply %s and %s' % (self, x))
        return self.__class__(self.expr * x, nid=self.nid)

    def __rmul__(self, x):
        return self.__class__(self.expr * x, nid=self.nid)    

    def __div__(self, x):
        if isinstance(x, noiseExpr) and self.nid != x.nid:
            raise ValueError('Cannot divide %s and %s' % (self, x))
        return self.__class__(self.expr / x, nid=self.nid)

    def __rdiv__(self, x):
        return self.__class__(x / self.expr, nid=self.nid)        

    def __eq__(self, x):
        try:
            if self.nid != x.nid:
                return False
        except:
            pass

        try:
            if self.expr == x.expr:
                return True
        except:
            pass
        return self.expr == x

    def __ne__(self, x):
        return not (self == x)

    def rms(self):
        """Calculate rms value."""

        P = sym.integrate(self.expr**2, (self.var, 0, sym.oo)) / (2 * sym.pi)
        rms = sym.sqrt(P)
        # TODO: Use rms class?
        return self._fourier_conjugate_class(rms)

    def sample(self, t):
        """Return a sample function (realisation) of the noise process
        evaluated at time values specified by vector t."""

        N = len(t)
        if N < 3:
            raise ValueError('Require at least 3 samples')            
        
        td = np.diff(t)
        if not np.allclose(np.diff(td), 0):
            raise ValueError('Require uniform sampling')

        fs = 1 / td[0]
        f = np.arange(N // 2 + 1) * fs / N
        Sn = self.evaluate(f * 2 * np.pi)
    
        x = np.random.randn(N)
        X = np.fft.rfft(x)
        
        Y = X * np.sqrt(Sn * fs / 2)
        y = np.fft.irfft(Y)
        return y

    def time(self):
        print('Warning: no time representation for noise expression'
              ', assumed zero: use rms()')
        return 0

    def laplace(self):
        print('Warning: no Laplace representation for noise expression'
              ', assumed zero')
        return 0    

    def autocorrelation(self):
        # Convert to two-sided spectrum
        S = self.subs(self.var, abs(self.var)) / sqrt(2)
        return S.inverse_fourier()

    def plot(self, fvector=None, **kwargs):

        """Plot frequency response at values specified by fvector.

        There are many plotting options, see matplotlib.pyplot.plot.

        For example:
            V.plot(fvector, log_frequency=True)
            V.real.plot(fvector, color='black')
            V.phase.plot(fvector, color='black', linestyle='--')

        By default complex data is plotted as separate plots of magnitude (dB)
        and phase.
        """

        from .plot import plot_frequency
        # Hack so show as linear frequency. 
        obj = self.subs(omega, 2 * pi * omega)
        obj.domain_name = 'Frequency'
        obj.domain_units = 'Hz'
        if hasattr(self, 'part'):
            obj.part = self.part
        return plot_frequency(obj, fvector, **kwargs)    

class Vn(noiseExpr):
    """Voltage noise amplitude spectral density (units V/rtHz).
    This can be a function of angular frequency, omega.  For example,
    to model an opamp voltage noise:

    v = Vn(1e-8 / sqrt(omega) + 8e-9)
    
    """

    quantity = 'Voltage noise spectral density'
    units = 'V/rtHz'

    def __init__(self, val, **assumptions):

        assumptions['positive'] = True
        super(Vn, self).__init__(val, **assumptions)
        # FIXME
        self._fourier_conjugate_class = Vt


class In(noiseExpr):
    """Current noise amplitude spectral density (units A/rtHz).

    This can be a function of angular frequency, omega.  For example,
    to model an opamp current noise:

    i = In(3e-12 / sqrt(omega) + 200e-15)
    """

    quantity = 'Current noise spectral density'
    units = 'A/rtHz'

    def __init__(self, val, **assumptions):

        assumptions['positive'] = True
        super(In, self).__init__(val, **assumptions)
        # FIXME
        self._fourier_conjugate_class = It

        
from .texpr import Ht, It, Vt, Yt, Zt, tExpr        
        
