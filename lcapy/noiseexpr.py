from __future__ import division
from .sym import symsimplify
from .functions import sqrt
from .sym import pi, omegasym, fsym
from .state import state
from .expr import Expr
import sympy as sym
import numpy as np

class NoiseExpression(Expr):
    """Frequency domain (one-sided) noise spectrum expression (amplitude
    spectral density).

    This characterises a zero-mean Gaussian noise process.

    When performing arithmetic on two noiseExpr expressions it is
    assumed that they are uncorrelated unless they have the same nid
    (noise indentifier).  If the nid is not specified, a new one is
    created.

    Uncorrelated noise expressions are added in quadrature (on a power
    basis).  Thus (NoiseExpression(3) + NoiseExpression(4)).expr = 5
    since 5 = sqrt(3**2 + 4**2)

    NoiseExpression(3) != NoiseExpression(3) since they are different
    noise realisations albeit with the same properties.  However,
    NoiseExpression(3).expr == NoiseExpression(3).expr.  Similarly,
    NoiseExpression(3, nid='n1') == NoiseExpression(3, nid='n1') since
    they have the same noise identifier and thus have the same
    realisation.

    Caution: The sum of two noise expressions generates a noise
    expression with a new nid.  This can lead to unexpected results
    since noise expressions with different nids are assumed to be
    uncorrelated.  For example, consider: 
    a = NoiseExpression(3); 
    b = NoiseExpression(4)
    a + b - b gives sqrt(41) and a + b - a gives sqrt(34).

    This case is correctly handled by the Super class since each noise
    component is stored and considered separately.

    (Voltage(a) + Voltage(b) - Voltage(b)).n gives 3 as expected.

    """
    is_real = True
    is_positive = True

    def _new_nid(self):
        state.context.nid += 1
        return 'n%d' % state.context.nid

    def __init__(self, val, **assumptions):
        if 'nid' not in assumptions or assumptions['nid'] is None:
            if isinstance(val, NoiseExpression):
                nid = val.nid
            elif val == 0:
                assumptions['nid'] = 'n0'
            else:
                assumptions['nid'] = self._new_nid()

        super(NoiseExpression, self).__init__(val, **assumptions)

    @property
    def nid(self):
        return self.assumptions['nid']

    def subs(self, *args, **kwargs):
        nid = self.nid
        result = super(NoiseExpression, self).subs(*args, **kwargs)
        result.assumptions['nid'] = nid
        return result

    def __compat__(self, x):

        if not isinstance(x, NoiseExpression):
            return x
        if self.var == x.var:
            return x
        return x(self.var)
    
    def __add__(self, x):
        """Add noise spectra (on power basis if uncorrelated)."""

        # Perhaps allow constant?
        
        if not isinstance(x, NoiseExpression):
            raise ValueError('Cannot add %s and %s' % (self, x))
        
        if x == 0:
            return self.__class__(self, nid=self.nid)

        x = self.__compat__(x)
        
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

        # Perhaps allow constant?
        
        if not isinstance(x, NoiseExpression):
            raise ValueError('Cannot subtract %s and %s' % (self, x))

        if x == 0:
            return self.__class__(self, nid=self.nid)

        x = self.__compat__(x)        

        if self.nid == x.nid:
            return self.__class__(self.expr - x.expr, nid=self.nid)        
        return self + x

    def __rsub__(self, x):
        raise ValueError('Cannot subtract %s and %s' % (self, x))          

    def __mul__(self, x):
        if isinstance(x, NoiseExpression) and self.nid != x.nid:
            raise ValueError('Cannot multiply %s and %s' % (self, x))

        x = self.__compat__(x)
        return self.__class__(self.expr * x, nid=self.nid)

    def __rmul__(self, x):
        return self.__class__(self.expr * x, nid=self.nid)    

    def __div__(self, x):
        if isinstance(x, NoiseExpression) and self.nid != x.nid:
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

        x = self.__compat__(x)        
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

        P = sym.integrate(self.expr**2, (self.var, 0, sym.oo))
        if self.var == omegasym:
            P /= 2 * sym.pi
        rms = sym.sqrt(P)
        # TODO: Use rms class?
        cls = self._class_by_quantity(self.quantity, 'time')
        return cls(rms)

    def sample(self, t, seed=None):
        """Return a sample function (realisation) of the noise process
        evaluated at time values specified by vector t.   If t is an integer,
        this returns `t` samples."""

        if isinstance(t, int):
            t = np.arange(t)
        
        N = len(t)
        if N < 3:
            raise ValueError('Require at least 3 samples')
        if N & 1:
            raise ValueError('Require even number of samples')

        if seed is not None:
            np.random.seed(seed)
        
        td = np.diff(t)
        if not np.allclose(np.diff(td), 0):
            raise ValueError('Require uniform sampling')

        fs = 1 / td[0]
        vf = np.arange(N // 2 + 1) * fs / N
        Sn = self(f).evaluate(vf)
    
        x = np.random.randn(N)
        X = np.fft.rfft(x)
        
        Y = X * np.sqrt(Sn * fs / 2)
        y = np.fft.irfft(Y)
        return y

    def time(self):
        print('Warning: no time representation for noise expression, '
              'assumed zero: use rms()')
        return 0

    def fourier(self):
        print('Warning: no Fourier representation for noise expression, '
              'assumed zero: use asd() or psd()')
        return 0    

    def laplace(self):
        print('Warning: no Laplace representation for noise expression, '
              'assumed zero')
        return 0

    def asd(self):
        return self

    def psd(self):
        # FIXME, the units are wrong
        return self ** 2

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
        from .plot import plot_frequency, plot_angular_frequency

        if self.var == fsym:
            return plot_frequency(self, fvector, **kwargs)
        elif self.var == omegasym:
            return plot_angular_frequency(self, fvector, **kwargs)
        raise ValueError('Invalid var for NoiseExpression')

    
from .fexpr import f
from .omegaexpr import omega
