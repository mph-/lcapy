from __future__ import division
from .sym import pi, omegasym, fsym
from .noiseomegaexpr import noiseomegaExpr


class Vnoisy(noiseomegaExpr):
    """Voltage noise amplitude spectral density (units V/rtHz).
    This can be a function of angular frequency, omega.  For example,
    to model an opamp voltage noise:

    v = Vnoisy(1e-8 / sqrt(omega) + 8e-9)
    
    """

    quantity = 'Voltage noise spectral density'
    units = 'V/rtHz'

    def __init__(self, val, **assumptions):

        assumptions['positive'] = True
        super(Vnoisy, self).__init__(val, **assumptions)
        # FIXME
        self._fourier_conjugate_class = Vt


class Inoisy(noiseomegaExpr):
    """Current noise amplitude spectral density (units A/rtHz).

    This can be a function of angular frequency, omega.  For example,
    to model an opamp current noise:

    i = Inoisy(3e-12 / sqrt(omega) + 200e-15)
    """

    quantity = 'Current noise spectral density'
    units = 'A/rtHz'

    def __init__(self, val, **assumptions):

        assumptions['positive'] = True
        super(Inoisy, self).__init__(val, **assumptions)
        # FIXME
        self._fourier_conjugate_class = It

from .texpr import It, Vt
        
