from __future__ import division
from .sym import symsimplify
from .functions import sqrt
from .sym import pi, omegasym, fsym
from .state import state
from .noiseexpr import noiseExpr
import sympy as sym
import numpy as np

class noisefExpr(noiseExpr):
    """Frequency domain (one-sided) noise spectrum expression (amplitude
    spectral density).

    This characterises a zero-mean Gaussian noise process.

    When performing arithmetic on two noiseExpr expressions it is
    assumed that they are uncorrelated unless they have the same nid
    (noise indentifier).  If the nid is not specified, a new one is
    created.

    Uncorrelated noise expressions are added in quadrature (on a power
    basis).  Thus (Vnoisy(3) + Vnoisy(4)).expr = 5 since 5 = sqrt(3**2 + 4**2)

    Vnoisy(3) != Vnoisy(3) since they are different noise realisations albeit
    with the same properties.  However, Vnoisy(3).expr == Vnoisy(3).expr.
    Similarly, Vnoisy(3, nid='n1') == Vnoisy(3, nid='n1') since they have the
    same noise identifier and thus have the same realisation.

    Caution: The sum of two noise expressions generates a noise
    expression with a new nid.  This can lead to unexpected results
    since noise expressions with different nids are assumed to be
    uncorrelated.  For example, consider:
    a = Vnoisy(3); b = Vnoisy(4)
    a + b - b gives sqrt(41) and  a + b - a gives sqrt(34).

    This case is correctly handled by the Super class since each noise
    component is stored and considered separately.

    (Voltage(a) + Voltage(b) - Voltage(b)).n gives 3 as expected.

    """
    one_sided = True
    var = fsym

    domain_name = 'Frequency'
    domain_units = 'Hz'    

    def __call__(self, arg, **assumptions):
        """Transform domain or substitute arg for variable. 
        
        Substitution is performed if arg is a tuple, list, numpy
        array, or constant.  If arg is a tuple or list return a list.
        If arg is an numpy array, return numpy array.

        Domain transformation is performed if arg is a domain variable
        or an expression of a domain variable.

        See also evaluate.

        """
        from .noiseomegaexpr import noiseomegaExpr        

        if id(arg) == id(self.var):
            return self.copy()

        if arg == omegasym:
            # Convert fsym to omegasym
            return noiseomegaExpr(self.subs(self.var, arg / (2 * pi)),
                                  nid=self.nid, **assumptions)
        return super(noisefExpr, self).__call__(arg, nid=self.nid, **assumptions)

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

        return plot_frequency(self, fvector, **kwargs)

    
from .texpr import It, Vt
from .fexpr import f
from .omegaexpr import omega
