"""
This module supports simple linear two-port networks.

Copyright 2014--2022 Michael Hayes, UCECE
"""

from __future__ import division
from warnings import warn
import sympy as sym
from .symbols import s
from .sexpr import LaplaceDomainVoltage, LaplaceDomainCurrent, LaplaceDomainImpedance
from .sexpr import LaplaceDomainAdmittance, LaplaceDomainTransferFunction
from .sexpr import LaplaceDomainExpression
from .smatrix import LaplaceDomainVoltageMatrix, LaplaceDomainCurrentMatrix
from .smatrix import LaplaceDomainImpedanceMatrix, LaplaceDomainAdmittanceMatrix
from .cexpr import ConstantDomainExpression
from .expr import expr
from .vector import Vector
from .matrix import Matrix
from .oneport import OnePort, I, V, Y, Z
from .network import Network
from .functions import exp, sqrt, Eq, MatMul, MatAdd

# TODO:
# 1. Defer the choice of the two-port model.  For example, a T section
# would store the three sub-networks rather than generating a B matrix.
# The appropriate model would be generated when desired.  This would
# avoid the inversion of singular matrices. The downside is that each
# object would require methods to generate each type of two-port model.
#
# Some multiport networks, such as a shunt R (or L, C) have a singular
# Y matrix.  Thus switching to the Z matrix and back to the Y matrix
# produces a bogus result.  The same thing occurs for a series R (or
# L, C); this has a singular Z matrix.
#
# 2. Fix handling of buffered two ports (amplifier / delay).
#
# 3. Inherit a subset of the Network methods

# Consider chaining a resistor shunt to a two-port network described by the
# A matrix A1.  The shunt has an A matrix A2.  The result has an A matrix A3.
#
# A1 = [A11 A12]
#      [A21 A22]
#
# A2 = [1     0]
#      [1/R   1]
#
# A3 = A1 * A2 = [A11 + A12/R  A12]
#                [A21 + A22/R  A22]
#
# We should get the same result by adding the Y matrices.
#
# Y1 = [A22/A12  -A11 A22 / A12 + A21]
#      [-1/A12                A11/A12]
#
# Y3 = [A22/A12  -A11 A22 / A12 + A21]
#      [-1/A12          A11/A12 + 1/R]
#
# but
#
# Y2 = [inf   inf]
#      [inf   inf]
#
# We can get the correct answer using:
#
# A2 = lim x->0  [1    x]
#                [1/R  1]
#
# Now det(A2) = lim x->0 (1 - x/R)
#
# and Y2 = lim x->0 [1/x   (1/x - 1/R)]
#                   [-1/x          1/x]
#
# Note when x=0, then
#
# Y2 = lim x->0 [1/0   1/0]
#               [-1/0  1/0]
#
# and we lose the information on R.
#
# The same problem occurs with a series R and the Z matrix.
#
# A2 = lim x->0  [1    R]
#                [x    1]
#
# Now det(A2) = lim x->0 (1 - R x)
#
# and Z2 = lim x->0 [1/x  1/x - R]
#                   [1/x      1/x]
#
# Thus it is advantageous to represent two-ports by the A (or B)
# matrix.  However, things will go wrong when we transform to the Y or
# Z matrix for specific cases.

__all__ = ('Chain', 'Par2', 'Ser2', 'Hybrid2', 'InverseHybrid2',
           'Series', 'Shunt', 'Transformer', 'IdealTransformer', 'IdealGyrator',
           'VoltageFollower', 'VoltageAmplifier',
           'IdealVoltageAmplifier', 'IdealDelay',
           'IdealVoltageDifferentiator', 'IdealVoltageIntegrator',
           'CurrentFollower', 'IdealCurrentAmplifier',
           'IdealCurrentDifferentiator', 'IdealCurrentIntegrator',
           'OpampInverter', 'OpampIntegrator', 'OpampDifferentiator',
           'TSection', 'TwinTSection', 'BridgedTSection', 'PiSection',
           'LSection', 'Ladder', 'GeneralTxLine', 'LosslessTxLine', 'TL',
           'TxLine', 'GeneralTransmissionLine', 'LosslessTransmissionLine',
           'TransmissionLine', 'AMatrix', 'BMatrix', 'GMatrix', 'HMatrix',
           'SMatrix', 'TMatrix', 'YMatrix', 'ZMatrix',
           'TwoPortAModel', 'TwoPortBModel', 'TwoPortYModel', 'TwoPortZModel',
           'TwoPortGModel', 'TwoPortHModel', 'TP',
           'TPA', 'TPB', 'TPG', 'TPH', 'TPY', 'TPZ')


def DeltaWye(Z1, Z2, Z3):

    ZZ = (Z1 * Z2 + Z2 * Z3 + Z3 * Z1)
    return (ZZ / Z1, ZZ / Z2, ZZ / Z3)


def WyeDelta(Z1, Z2, Z3):

    ZZ = Z1 + Z2 + Z3
    return (Z2 * Z3 / ZZ, Z1 * Z3 / ZZ, Z1 * Z2 / ZZ)


def _check_oneport_args(args):

    for arg1 in args:
        if not isinstance(arg1, OnePort):
            raise ValueError('%s not a OnePort' % arg1)


class TwoPortMixin(object):

    @property
    def is_buffered(self):
        """Return true if two-port is buffered, i.e., any load
        on the output has no affect on the input. """
        # return self._A12 == 0 and self._A22 == 0
        return self._B12 == 0 and self._B22 == 0

    @property
    def is_bilateral(self):
        """Return true if two-port is bilateral. """
        return self.Bparams.det() == 1

    @property
    def is_reciprocal(self):
        """Return true if two-port is reciprocal. """
        # This also applies to Y12 == Y21, S12 == S21, or Aparams.det() == 1.
        return self._Z12 == self._Z21

    @property
    def is_reversible(self):
        """Return true if two-port is reversible. """
        return self.params.det() != 0

    @property
    def is_symmetrical(self):
        """Return true if two-port is symmetrical. """
        return self._B11 == self._B22

    @property
    def is_series(self):
        """Return true if two-port is a series network. """
        # return (self._A11 == 1) and (self._A22 == 1) and (self._A21 == 0)
        return (self._B11 == 1) and (self._B22 == 1) and (self._B21 == 0)

    @property
    def is_shunt(self):
        """Return true if two-port is a shunt network. """
        # return (self._A11 == 1) and (self._A22 == 1) and (self._A12 == 0)
        return (self._B11 == 1) and (self._B22 == 1) and (self._B12 == 0)

    @property
    def _A11(self):
        """Open-circuit inverse voltage ratio"""
        return self.Aparams[0, 0]

    @property
    def _A12(self):
        """Negative short-circuit transfer impedance"""
        return self.Aparams[0, 1]

    @property
    def _A21(self):
        """Negative short-circuit inverse current ratio"""
        return self.Aparams[1, 0]

    @property
    def _A22(self):
        """Open circuit transfer admittance"""
        return self.Aparams[1, 1]

    @property
    def _B11(self):
        """Open-circuit voltage gain"""
        return self.Bparams[0, 0]

    @property
    def _B12(self):
        """Negative short-circuit transfer impedance"""
        return self.Bparams[0, 1]

    @property
    def _B21(self):
        """Negative short-circuit current gain"""
        return self.Bparams[1, 0]

    @property
    def _B22(self):
        """Open-circuit transfer admittance"""
        return self.Bparams[1, 1]

    @property
    def _G11(self):
        """Open-circuit input admittance"""
        return self.Gparams[0, 0]

    @property
    def _G12(self):
        """Short-circuit reverse current gain"""
        return self.Gparams[0, 1]

    @property
    def _G21(self):
        """Open-circuit forward voltage gain"""
        return self.Gparams[1, 0]

    @property
    def _G22(self):
        """Short-circuit output impedance"""
        return self.Gparams[1, 1]

    @property
    def _H11(self):
        """Short-circuit input impedance"""
        return self.Hparams[0, 0]

    @property
    def _H12(self):
        """Open-circuit reverse voltage gain"""
        return self.Hparams[0, 1]

    @property
    def _H21(self):
        """Short-circuit forward current gain"""
        return self.Hparams[1, 0]

    @property
    def _H22(self):
        """Open-circuit output admittance"""
        return self.Hparams[1, 1]

    @property
    def _S11(self):
        return self.Sparams[0, 0]

    @property
    def _S12(self):
        return self.Sparams[0, 1]

    @property
    def _S21(self):
        return self.Sparams[1, 0]

    @property
    def _S22(self):
        return self.Sparams[1, 1]

    @property
    def _T11(self):
        return self.Tparams[0, 0]

    @property
    def _T12(self):
        return self.Tparams[0, 1]

    @property
    def _T21(self):
        return self.Tparams[1, 0]

    @property
    def _T22(self):
        return self.Tparams[1, 1]

    @property
    def _Y11(self):
        """Short-circuit input admittance"""
        return self.Yparams[0, 0]

    @property
    def _Y12(self):
        """Short-circuit reverse transfer admittance"""
        return self.Yparams[0, 1]

    @property
    def _Y21(self):
        """Short-circuit transfer admittance"""
        return self.Yparams[1, 0]

    @property
    def _Y22(self):
        """Short-circuit output admittance"""
        return self.Yparams[1, 1]

    @property
    def _Z11(self):
        """Open-cicuit input impedance"""
        return self.Zparams[0, 0]

    @property
    def _Z12(self):
        """Open-cicuit transfer impedance"""
        return self.Zparams[0, 1]

    @property
    def _Z21(self):
        """Open-cicuit reverse transfer impedance"""
        return self.Zparams[1, 0]

    @property
    def _Z22(self):
        """Open-cicuit output impedance"""
        return self.Zparams[1, 1]

    @property
    def A11(self):
        """Open-circuit inverse voltage ratio"""
        return expr(self._A11)

    @property
    def A12(self):
        """Negative short-circuit transfer impedance"""
        return expr(self._A12)

    @property
    def A21(self):
        """Negative short-circuit inverse current ratio"""
        return expr(self._A21)

    @property
    def A22(self):
        """Open circuit transfer admittance"""
        return expr(self._A22)

    @property
    def B11(self):
        """Open-circuit voltage gain"""
        return expr(self._B11)

    @property
    def B12(self):
        """Negative short-circuit transfer impedance"""
        return expr(self._B12)

    @property
    def B21(self):
        """Negative short-circuit current gain"""
        return expr(self._B21)

    @property
    def B22(self):
        """Open-circuit transfer admittance"""
        return expr(self._B22)

    @property
    def G11(self):
        """Open-circuit input admittance"""
        return expr(self._G11)

    @property
    def G12(self):
        """Short-circuit reverse current gain"""
        return expr(self._G12)

    @property
    def G21(self):
        """Open-circuit forward voltage gain"""
        return expr(self._G21)

    @property
    def G22(self):
        """Short-circuit output impedance"""
        return expr(self._G22)

    @property
    def H11(self):
        """Short-circuit input impedance"""
        return expr(self._H11)

    @property
    def H12(self):
        """Open-circuit reverse voltage gain"""
        return expr(self._H12)

    @property
    def H21(self):
        """Short-circuit forward current gain"""
        return expr(self._H21)

    @property
    def H22(self):
        """Open-circuit output admittance"""
        return expr(self._H22)

    @property
    def S11(self):
        """S11"""
        return expr(self._S11)

    @property
    def S12(self):
        """S12"""
        return expr(self._S12)

    @property
    def S21(self):
        """S21"""
        return expr(self._S21)

    @property
    def S22(self):
        """S22"""
        return expr(self._S22)

    @property
    def T11(self):
        """T11"""
        return expr(self._T11)

    @property
    def T12(self):
        """T12"""
        return expr(self._T12)

    @property
    def T21(self):
        """T21"""
        return expr(self._T21)

    @property
    def T22(self):
        """T22"""
        return expr(self._T22)

    @property
    def Y11(self):
        """Short-circuit input admittance"""
        return expr(self._Y11)

    @property
    def Y12(self):
        """Short-circuit reverse transfer admittance"""
        return expr(self._Y12)

    @property
    def Y21(self):
        """Short-circuit transfer admittance"""
        return expr(self._Y21)

    @property
    def Y22(self):
        """Short-circuit output admittance"""
        return expr(self._Y22)

    @property
    def Z11(self):
        """Open-cicuit input impedance"""
        return expr(self._Z11)

    @property
    def Z12(self):
        """Open-cicuit transfer impedance"""
        return expr(self._Z12)

    @property
    def Z21(self):
        """Open-cicuit reverse transfer impedance"""
        return expr(self._Z21)

    @property
    def Z22(self):
        """Open-cicuit output impedance"""
        return expr(self._Z22)

    @property
    def voltage_gain(self):
        """Return V2 / V1 for I2 = 0 with internal sources zero.

        This is an alias for forward_voltage_gain."""

        return self.Vgain12

    @property
    def forward_voltage_gain(self):
        """Return V2 / V1 for I2 = 0 with internal sources zero."""

        return self.Vgain12

    @property
    def reverse_voltage_gain(self):
        """Return V1 / V2 for I1 = 0 with internal sources zero."""

        return self.Vgain21

    @property
    def current_gain(self):
        """Return I2 / I1 for V2 = 0 with internal sources zero.

        This is an alias for forward_current_gain."""

        return self.Igain12

    @property
    def forward_current_gain(self):
        """Return I2 / I1 for V2 = 0 with internal sources zero."""

        return self.Igain12

    @property
    def reverse_current_gain(self):
        """Return I1 / I2 for I2 = 0 with internal sources zero."""

        return self.Igain21

    @property
    def transadmittance(self):
        """Return I2 / V1 for V2 = 0 with internal sources zero.

        This is an alias for forward_transadmittance."""

        return LaplaceDomainAdmittance(self._Y21)

    @property
    def forward_transadmittance(self):
        """Return I2 / V1 for V2 = 0 with internal sources zero."""

        return LaplaceDomainAdmittance(self._Y21)

    @property
    def reverse_transadmittance(self):
        """Return I1 / V2 for V1 = 0 with internal sources zero."""

        return LaplaceDomainAdmittance(self._Y12)

    @property
    def transimpedance(self):
        """Return V2 / I1 for I2 = 0 with internal sources zero.

        This is an alias for forward_transadmittance."""

        return LaplaceDomainImpedance(self._Z21)

    @property
    def forward_transimpedance(self):
        """Return V2 / I1 for I2 = 0 with internal sources zero."""

        return LaplaceDomainImpedance(self._Z21)

    @property
    def reverse_transimpedance(self):
        """Return V1 / I2 for I1 = 0 with internal sources zero."""

        return LaplaceDomainImpedance(self._Z12)


class TwoPortMatrix(Matrix, TwoPortMixin):

    # The following default properties are fallbacks when other conversions have
    # not been defined.

    @property
    def Aparams(self):
        return AMatrix(self.Bparams.inv()).simplify()

    @property
    def Bparams(self):
        if not hasattr(self, '_Bparams'):
            self._Bparams = BMatrix(self.Aparams.inv()).simplify()
        return self._Bparams

    @property
    def Gparams(self):
        return GMatrix(self.Hparams.inv()).simplify()

    @property
    def Hparams(self):
        return HMatrix(self.Gparams.inv()).simplify()

    @property
    def Sparams(self):
        return self.Aparams.Sparams

    @property
    def Tparams(self):
        return self.Sparams.Tparams

    @property
    def Yparams(self):
        return YMatrix(self.Zparams.inv()).simplify()

    @property
    def Zparams(self):
        return ZMatrix(self.Yparams.inv()).simplify()

    def pdb(self):
        import pdb
        pdb.set_trace()
        return self

    @property
    def Z1oc(self):
        """Open-circuit input impedance"""
        return self.Aparams.Z1oc

    @property
    def Z1sc(self):
        """Short-circuit input impedance"""
        return self.Aparams.Z1sc

    @property
    def Z2oc(self):
        """Open-circuit output impedance"""
        return self.Aparams.Z2oc

    @property
    def Z2sc(self):
        """Short-circuit output impedance"""
        return self.Aparams.Z2sc

    @property
    def Vgain12(self):
        """Forward voltage gain"""
        return self.Aparams.Vgain12

    @property
    def Vgain21(self):
        """Reverse voltage gain"""
        return self.Aparams.Vgain21

    @property
    def Igain12(self):
        """Forward current gain"""
        return self.Aparams.Igain12

    @property
    def Igain21(self):
        """Reverse current gain"""
        return self.Aparams.Igain21


class AMatrix(TwoPortMatrix):
    """A-parameters (ABCD parameters, chain matrix)
    ::
       +-  -+     +-       -+   +-  -+
       | V1 |  =  | A11  A12|   | V2 |
       | I1 |     | A21  A22|   |-I2 |
       +-  -+     +-       -+   +-  -+

              +-         -+
       units  | 1     ohm |
              | 1/ohm   1 |
              +-         -+

    A buffered two-port has A12 = A22 = 0.

    A = inv(B)
    """

    @classmethod
    def generic(cls):
        return cls((('A_11', 'A_12'), ('A_21', 'A_22')))

    def equation(self):
        return Eq(Matrix(('V1', 'I1')), MatMul(self, Matrix(('V2', '-I2'))),
                  evaluate=False)

    @property
    def Aparams(self):
        # Perhaps we should make a copy?
        return self

    @property
    def Bparams(self):
        if not hasattr(self, '_Bparams'):
            # Inverse
            det = self.det().expr
            if det == 0:
                warn('Producing dodgy B matrix')
            self._Bparams = BMatrix(self.inv()).simplify()
        return self._Bparams

    @property
    def Hparams(self):

        if self._A22 == 0:
            warn('Producing dodgy H matrix')
        det = self.det().expr
        return HMatrix(((self._A12 / self._A22, det / self._A22),
                        (-1 / self._A22, self._A21 / self._A22))).simplify()

    @property
    def Sparams(self):
        Z0 = LaplaceDomainImpedance('Z_0').as_expr()
        d = self._A12 + Z0 * (self._A11 + self._A22) + Z0**2 * self._A21
        return SMatrix((((self._A12 + Z0 * (self._A11 - self._A22) - Z0**2 * self._A21) / d,
                         (2 * Z0 * (self._A11 * self._A22 - self._A12 * self._A21)) / d),
                        (2 * Z0 / d, (self._A12 - Z0 * (self._A11 - self._A22) - Z0**2 * self._A21) / d))).simplify()

    @property
    def Yparams(self):

        # This produces a bogus Y matrix when A12 is zero (say for a
        # shunt element).   Note, it doesn't use A21.
        if self._A12 == 0:
            warn('Producing dodgy Y matrix')
        det = self.det().expr
        return YMatrix(((self._A22 / self._A12, -det / self._A12),
                        (-1 / self._A12, self._A11 / self._A12))).simplify()

    @property
    def Zparams(self):

        # This produces a bogus Z matrix when A21 is zero (say for a
        # series element).   Note, it doesn't use A12.
        if self._A21 == 0:
            warn('Producing dodgy Z matrix')
        det = self.det().expr
        return ZMatrix(((self._A11 / self._A21, det / self._A21),
                        (1 / self._A21, self._A22 / self._A21))).simplify()

    @property
    def Z1oc(self):
        """Open-circuit input impedance"""
        return LaplaceDomainImpedance(self._A11 / self._A21)

    @property
    def Z1sc(self):
        """Short-circuit input impedance"""
        return LaplaceDomainImpedance(self._A12 / self._A22)

    @property
    def Z2oc(self):
        """Open-circuit output impedance"""
        return LaplaceDomainImpedance(self._A22 / self._A21)

    @property
    def Z2sc(self):
        """Short-circuit output impedance"""
        return LaplaceDomainImpedance(self._A12 / self._A11)

    @property
    def Vgain12(self):
        """Forward voltage gain"""
        return LaplaceDomainTransferFunction(1 / self._A11)

    @property
    def Vgain21(self):
        """Reverse voltage gain"""
        return LaplaceDomainTransferFunction((self._A11 * self._A22 - self._A12 * self._A21) / self.A22)

    @property
    def Igain12(self):
        """Forward current gain"""
        return LaplaceDomainTransferFunction(-1 / self._A22)

    @property
    def Igain21(self):
        """Reverse current gain"""
        return LaplaceDomainTransferFunction(-(self._A11 * self._A22 - self._A12 * self._A21) / self.A11)

    @property
    def forward_transadmittance(self):
        """Return I2 / V1 for V2 = 0 with internal sources zero."""

        return LaplaceDomainAdmittance(-1 / self._A12)

    @property
    def reverse_transadmittance(self):
        """Return I1 / V2 for V1 = 0 with internal sources zero."""

        return LaplaceDomainAdmittance(self._A21 - self._A11 * self._A22 / self._A12)

    @classmethod
    def Zseries(cls, Zval):

        if not isinstance(Zval, LaplaceDomainImpedance):
            raise ValueError('Zval not LaplaceDomainImpedance')

        return cls(((1, Zval),
                    (0, 1)))

    @classmethod
    def Yseries(cls, Yval):

        if not isinstance(Yval, LaplaceDomainAdmittance):
            raise ValueError('Yval not LaplaceDomainAdmittance')

        return cls(((1, 1 / Yval),
                    (0, 1)))

    @classmethod
    def Yshunt(cls, Yval):

        if not isinstance(Yval, LaplaceDomainAdmittance):
            raise ValueError('Yval not LaplaceDomainAdmittance')

        return cls(((1, 0),
                    (Yval, 1)))

    @classmethod
    def Zshunt(cls, Zval):

        if not isinstance(Zval, LaplaceDomainImpedance):
            raise ValueError('Zval not LaplaceDomainImpedance')

        return cls(((1, 0),
                    (1 / Zval, 1)))

    @classmethod
    def transformer(cls, alpha):
        """The voltage gain alpha = 1 / a, where a is the turns ratio."""

        alpha = ConstantDomainExpression(alpha)

        return cls(((1 / alpha, 0),
                    (0, alpha)))

    @classmethod
    def gyrator(cls, R):

        R = ConstantDomainExpression(R)

        return cls(((0, R),
                    (1 / R, 0)))

    @classmethod
    def Lsection(cls, Z1, Z2):

        return cls.Zseries(Z1).chain(cls.Zshunt(Z2))

    @classmethod
    def Tsection(cls, Z1, Z2, Z3):

        return cls.Lsection(Z1, Z2).chain(cls.Zseries(Z3))

    @classmethod
    def Pisection(cls, Z1, Z2, Z3):

        return cls.Zshunt(Z1).chain(cls.Lsection(Z2, Z3))

    def chain(self, TP):

        return self * TP

    def cascade(self, TP):

        return self.chain(TP)


class BMatrix(TwoPortMatrix):
    """B-parameters (inverse ABCD parameters)
    ::
       +-  -+     +-       -+   +-  -+
       | V2 |  =  | B11  B12|   | V1 |
       |-I2 |     | B21  B22|   | I1 |
       +-  -+     +-       -+   +-  -+

              +-         -+
       units  | 1     ohm |
              | 1/ohm   1 |
              +-         -+

    B = inv(A)
    """

    @classmethod
    def generic(cls):
        return cls((('B_11', 'B_12'), ('B_21', 'B_22')))

    def equation(self):
        return Eq(Matrix(('V2', '-I2')), MatMul(self, Matrix(('V1', 'I1'))),
                  evaluate=False)

    @property
    def Aparams(self):
        return AMatrix(self.inv()).simplify()

    @property
    def Bparams(self):
        # Perhaps we should make a copy?
        return self

    @property
    def Gparams(self):

        det = self.det().expr
        return GMatrix(((-self._B21 / self._B22, -1 / self._B22),
                        (det / self._B22, -self._B12 / self._B22))).simplify()

    @property
    def Hparams(self):

        return HMatrix(((-self._B12 / self._B11, 1 / self._B11),
                        (self._B21 * self._B12 / self._B11 - self._B22,
                         -self._B21 / self._B11))).simplify()

    @property
    def Yparams(self):

        det = self.det().expr
        return YMatrix(((-self._B11 / self._B12, 1 / self._B12),
                        (det / self._B12, -self._B22 / self._B12))).simplify()

    @property
    def Zparams(self):

        det = self.det().expr
        return ZMatrix(((-self._B22 / self._B21, -1 / self._B21),
                        (-det / self._B21, -self._B11 / self._B21))).simplify()

    @property
    def Z1oc(self):
        """Open-circuit input impedance"""
        return LaplaceDomainImpedance(-self._B22 / self._B21)

    @property
    def Z1sc(self):
        """Short-circuit input impedance"""
        return LaplaceDomainImpedance(-self._B12 / self._B11)

    @property
    def Z2oc(self):
        """Open-circuit output impedance"""
        return LaplaceDomainImpedance(-self._B11 / self._B21)

    @property
    def Z2sc(self):
        """Short-circuit output impedance"""
        return LaplaceDomainImpedance(-self._B12 / self._B22)

    @property
    def Vgain12(self):
        """Forward voltage gain"""
        return LaplaceDomainTransferFunction((self._B11 * self._B22 - self._B12 * self._B21) / self.B22)

    @property
    def Vgain21(self):
        """Reverse voltage gain"""
        return LaplaceDomainTransferFunction(1 / self._B11)

    @property
    def Igain12(self):
        """Forward current gain"""
        return LaplaceDomainTransferFunction(-(self._B11 * self._B22 - self._B12 * self._B21) / self.B11)

    @property
    def Igain21(self):
        """Reverse current gain"""
        return LaplaceDomainTransferFunction(-1 / self._B22)

    @property
    def forward_transadmittance(self):
        """Return I2 / V1 for V2 = 0 with internal sources zero."""

        return LaplaceDomainAdmittance(self._B11 * self._B22 / self._B12 - self._B21)

    @property
    def reverse_transadmittance(self):
        """Return I1 / V2 for V1 = 0 with internal sources zero."""

        return LaplaceDomainAdmittance(-1 / self._B12)

    @classmethod
    def Zseries(cls, Zval):

        if not isinstance(Zval, LaplaceDomainImpedance):
            raise ValueError('Zval not LaplaceDomainImpedance')

        return cls(((1, -Zval),
                    (0, 1)))

    @classmethod
    def Yseries(cls, Yval):

        if not isinstance(Yval, LaplaceDomainAdmittance):
            raise ValueError('Yval not LaplaceDomainAdmittance')

        return cls(((1, -1 / Yval),
                    (0, 1)))

    @classmethod
    def Yshunt(cls, Yval):

        if not isinstance(Yval, LaplaceDomainAdmittance):
            raise ValueError('Yval not LaplaceDomainAdmittance')

        return cls(((1, 0),
                    (-Yval, 1)))

    @classmethod
    def Zshunt(cls, Zval):

        if not isinstance(Zval, LaplaceDomainImpedance):
            raise ValueError('Zval not LaplaceDomainImpedance')

        return cls(((1, 0),
                    (-1 / Zval, 1)))

    @classmethod
    def voltage_amplifier(cls, Af, Ar=1e-9, Yin=1e-9, Zout=1e-9):
        """Voltage amplifier
        Af forward voltage gain
        Ar reverse voltage gain (ideally 0)
        Yin input admittance (ideally 0)
        Zout output impedance (ideally 0)
        """

        if Ar == 0 and Yin == 0 and Zout == 0:
            warn('Should use G matrix; tweaking B matrix to make invertible')
            Ar = 1e-9
            Yin = 1e-9
            Zout = 1e-9

        Af = LaplaceDomainExpression(Af)
        Ar = LaplaceDomainExpression(Ar)
        Yin = LaplaceDomainExpression(Yin)
        Zout = LaplaceDomainExpression(Zout)

        # This should be defined with a G matrix
        #
        # G = [0   0]
        #     [Af  0]
        #
        # With this model, the inverse voltage gain is 1 / Af
        #
        # G = lim x->0  [0   x]
        #               [Af  0]
        #
        # B = lim x->0  [Af    0/x]
        #               [0/x  -1/x]
        #
        # A = lim x->0  [1/Af 0/Af]
        #               [0/Af   -x]

        # Perhaps default Ar, Yin, and Zout to 1e-10 to get a reasonable
        # B matrix?

        return cls(((1 / Ar, -1 / (Ar * Yin)),
                    (-1 / (Ar * Zout), -1 / (Ar * Yin * Zout * (Af * Ar - 1)))))

    @classmethod
    def current_amplifier(cls, Af, Ar=1e-9, Zin=1e-9, Yout=1e-9):
        """Current amplifier
        Af forward current gain
        Ar reverse current gain (ideally 0)
        Yin input admittance (ideally 0)
        Yout output impedance (ideally 0)
        """

        if Ar == 0 and Zin == 0 and Yout == 0:
            warn('Should use G matrix; tweaking B matrix to make invertible')
            Ar = 1e-9
            Zin = 1e-9
            Yout = 1e-9

        Af = LaplaceDomainExpression(Af)
        Ar = LaplaceDomainExpression(Ar)
        Zin = LaplaceDomainExpression(Zin)
        Yout = LaplaceDomainExpression(Yout)

        # This should be defined with a H matrix
        #
        # H = [0   0]
        #     [Af  0]
        #
        # With this model, the inverse current gain is 1 / Af
        #
        # H = lim x->0  [0   x]
        #               [Af  0]
        #
        # B = lim x->0  [1/x    0]
        #               [0    -Af]
        #
        # A = lim x->0  [1/x  -0/x]
        #               [0/x   -Af]

        return cls(((1 / Ar, -1 / (Ar * Yout)),
                    (-1 / (Ar * Zin), -1 / (Ar * Yout * Zin * (Af * Ar - 1)))))

    @classmethod
    def voltage_differentiator(cls, Av=1):

        return cls.voltage_amplifier(LaplaceDomainExpression(Av).differentiate())

    @classmethod
    def voltage_integrator(cls, Av):

        return cls.voltage_amplifier(LaplaceDomainExpression(Av).integrate())

    @classmethod
    def current_differentiator(cls, Av):

        return cls.current_amplifier(LaplaceDomainExpression(Av).differentiate())

    @classmethod
    def current_integrator(cls, Av):

        return cls.current_amplifier(LaplaceDomainExpression(Av).integrate())

    @classmethod
    def transformer(cls, alpha):
        """The voltage gain alpha = 1 / a, where a is the turns ratio."""

        alpha = expr(alpha)

        return cls(((alpha, 0),
                    (0, 1 / alpha)))

    @classmethod
    def gyrator(cls, R):

        R = ConstantDomainExpression(R)

        return cls(((0, R),
                    (1 / R, 0)))

    @classmethod
    def Lsection(cls, Z1, Z2):

        Y = 1 / Z2
        return cls(((1 + Y * Z1, -Z1),
                    (-Y, 1)))
        # return cls.Zseries(Z1).chain(cls.Zshunt(Z2))

    @classmethod
    def Tsection(cls, Z1, Z2, Z3):

        Y = 1 / Z2
        return cls(((1 + Y * Z1, -Z1 - Z3 * (1 + Y * Z1)),
                    (-Y, 1 + Y * Z3)))
        # return cls.Lsection(Z1, Z2).chain(cls.Zseries(Z3))

    @classmethod
    def Pisection(cls, Z1, Z2, Z3):

        return cls.Zshunt(Z1).chain(cls.Lsection(Z2, Z3))

    def chain(self, TP):

        # Note reverse order compared to AMatrix.
        return TP * self

    def cascade(self, TP):

        return self.chain(TP)


class GMatrix(TwoPortMatrix):
    """G-parameters (inverse hybrid parameters)

    ::
       +-  -+     +-       -+   +-  -+
       | V2 |  =  | G11  G12|   | I2 |
       | I1 |     | G21  G22|   | V1 |
       +-  -+     +-       -+   +-  -+

              +-         -+
       units  | ohm     1 |
              | 1   1/ohm |
              +-         -+

    G = inv(H)
    """

    @classmethod
    def generic(cls):
        return cls((('G_11', 'G_12'), ('G_21', 'G_22')))

    def equation(self):
        return Eq(Matrix(('I1', 'V2')), MatMul(self, Matrix(('V1', 'I2'))),
                  evaluate=False)

    @property
    def Aparams(self):
        # return self.Hparams.Aparams
        det = self.det().expr
        return AMatrix(((1 / self._G21, self._G22 / self._G21),
                        (self._G11 / self._G21, det / self._G21))).simplify()

    @property
    def Bparams(self):
        if not hasattr(self, '_Bparams'):
            det = self.det().expr
            self._Bparams = BMatrix(((-det / self._G12, (self._G22 / self._G12)),
                                     (self._G11 / self._G12, -1 / self._G12))).simplify()
        return self._Bparams

    @property
    def Gparams(self):
        # Perhaps we should make a copy?
        return self

    @property
    def Hparams(self):
        return HMatrix(self.inv()).simplify()

    @property
    def Yparams(self):
        return self.Hparams.Yparams

    @property
    def Zparams(self):
        return self.Hparams.Zparams


class HMatrix(TwoPortMatrix):
    """H-parameters (hybrid parameters)
    ::
       +-  -+     +-       -+   +-  -+
       | V1 |  =  | H11  H12|   | I1 |
       | I2 |     | H21  H22|   | V2 |
       +-  -+     +-       -+   +-  -+

              +-         -+
       units  | ohm     1 |
              | 1   1/ohm |
              +-         -+

    H = inv(G)
    """

    @classmethod
    def generic(cls):
        return cls((('H_11', 'H_12'), ('H_21', 'H_22')))

    def equation(self):
        return Eq(Matrix(('V1', 'I2')), MatMul(self, Matrix(('I1', 'V2'))),
                  evaluate=False)

    @property
    def Aparams(self):
        det = self.det().expr
        return AMatrix(((-det / self._H21, -self._H11 / self._H21),
                        (-self._H22 / self._H21, -1 / self._H21))).simplify()

    @property
    def Bparams(self):
        if not hasattr(self, '_Bparams'):
            self._Bparams = BMatrix(((1 / self._H12, -self._H11 / self._H12),
                                     (-self._H22 / self._H12,
                                     self._H22 * self._H11 / self._H12 - self._H21))).simplify()

        return self._Bparams

    @property
    def Hparams(self):
        # Perhaps we should make a copy?
        return self

    @property
    def Yparams(self):
        det = self.det().expr
        return YMatrix(((1 / self._H11, -self._H12 / self._H11),
                        (self._H21 / self._H11, det / self._H11))).simplify()

    @property
    def Zparams(self):
        det = self.det().expr
        return ZMatrix(((det / self._H22, self._H12 / self._H22),
                        (-self._H21 / self._H22, 1 / self._H22))).simplify()


class SMatrix(TwoPortMatrix):
    """S-parameters (scattering parameters)
    ::
       +-  -+     +-       -+   +-  -+
       | b1 |  =  | S11  S12|   | a1 |
       | b2 |     | S21  S22|   | a2 |
       +-  -+     +-       -+   +-  -+

    Each element in the s-matrix has units of impedance.
    """

    @classmethod
    def generic(cls):
        return cls((('S_11', 'S_12'), ('S_21', 'S_22')))

    def equation(self):
        return Eq(Matrix(('b1', 'b2')), MatMul(self, Matrix(('a1', 'a2'))),
                  evaluate=False)

    @property
    def Aparams(self):
        det = self.det().expr
        Z0 = LaplaceDomainImpedance('Z_0').as_expr()
        A = AMatrix(((1 + self._S11 - self._S22 - det,
                      (1 + self._S11 + self._S22 + det) * Z0),
                     ((1 - self._S11 - self._S22 + det) / Z0,
                      1 - self._S11 + self._S22 - det))) / (2 * self._S21)
        return A.simplify()

    @property
    def Hparams(self):
        return self.Aparams.Hparams

    @property
    def Sparams(self):
        return self

    @property
    def Tparams(self):
        det = self.det().expr
        return TMatrix(((-det / self._S21, self._S11 / self._S21),
                        (-self._S22 / self._S21, 1 / self._S21))).simplify()

    @property
    def Zparams(self):
        return self.Aparams.Zparams


class TMatrix(TwoPortMatrix):
    """T-parameters (scattering transfer parameters)
    ::
       +-  -+     +-       -+   +-  -+
       | b1 |  =  | T11  T12|   | a2 |
       | a2 |     | T21  T22|   | b2 |
       +-  -+     +-       -+   +-  -+
    """

    @classmethod
    def generic(cls):
        return cls((('T_11', 'T_12'), ('T_21', 'T_22')))

    # Note, another convention uses a1, b1 in terms of b2, a2.
    def equation(self):
        return Eq(Matrix(('b1', 'a1')), MatMul(self, Matrix(('a2', 'b2'))),
                  evaluate=False)

    @property
    def Aparams(self):
        return self.Sparams.Aparams

    @property
    def Hparams(self):
        return self.Aparams.Hparams

    @property
    def Sparams(self):
        det = self.det().expr
        return SMatrix(((self._T12 / self._T22, det / self._T22),
                        (1 / self._T22, -self._T21 / self._T22))).simplify()

    @property
    def Tparams(self):
        return self

    @property
    def Zparams(self):
        return self.Aparams.Zparams


class YMatrix(TwoPortMatrix):
    """Y-parameters (admittance parameters)
    ::
       +-  -+     +-       -+   +-  -+
       | I1 |  =  | Y11  Y12|   | V1 |
       | I2 |     | Y21  Y22|   | V2 |
       +-  -+     +-       -+   +-  -+

              +-           -+
       units  | 1/ohm 1/ohm |
              | 1/ohm 1/ohm |
              +-           -+

    Y = inv(Z)
    """

    @classmethod
    def generic(cls):
        return cls((('Y_11', 'Y_12'), ('Y_21', 'Y_22')))

    def equation(self):
        return Eq(Matrix(('I1', 'I2')), MatMul(self, Matrix(('V1', 'V2'))),
                  evaluate=False)

    @property
    def Ysc(self):
        return LaplaceDomainAdmittanceMatrix((self._Y11, self._Y22))

    @property
    def Aparams(self):
        det = self.det().expr
        return AMatrix(((-self._Y22 / self._Y21, -1 / self._Y21),
                        (-det / self._Y21, -self._Y11 / self._Y21))).simplify()

    @property
    def Bparams(self):
        if not hasattr(self, '_Bparams'):
            det = self.det().expr
            self._Bparams = BMatrix(((-self._Y11 / self._Y12, 1 / self._Y12),
                                     (det / self._Y12, -self._Y22 / self._Y12))).simplify()
        return self._Bparams

    @property
    def Hparams(self):
        det = self.det().expr
        return HMatrix(((1 / self._Y11, -self._Y12 / self._Y11),
                        (self._Y21 / self._Y11, det / self._Y11))).simplify()

    @property
    def Yparams(self):
        # Perhaps we should make a copy?
        return self

    @property
    def Zparams(self):
        # Inverse
        det = self.det().expr
        return ZMatrix(((self._Y22 / det, -self._Y12 / det),
                        (-self._Y21 / det, self._Y11 / det))).simplify()

    @property
    def Z1sc(self):
        """Short-circuit input impedance"""
        return LaplaceDomainImpedance(1 / self._Y11)

    @property
    def Z2sc(self):
        """Short-circuit input impedance"""
        return LaplaceDomainImpedance(1 / self._Y22)

    @property
    def Vgain12(self):
        """Forward voltage gain"""
        return LaplaceDomainTransferFunction(-self._Y21 / self._Y22)

    @property
    def Vgain21(self):
        """Reverse voltage gain"""
        return LaplaceDomainTransferFunction(-self._Y12 / self._Y11)

    @property
    def Igain12(self):
        """Forward current gain"""
        return LaplaceDomainTransferFunction(self._Y21 / self._Y11)

    @property
    def Igain21(self):
        """Reverse current gain"""
        return LaplaceDomainTransferFunction(self._Y12 / self._Y22)


class ZMatrix(TwoPortMatrix):
    """Z-parameters (impedance parameters)
    ::
       +-  -+     +-       -+   +-  -+
       | V1 |  =  | Z11  Z12|   | I1 |
       | V2 |     | Z21  Z22|   | I2 |
       +-  -+     +-       -+   +-  -+

              +-         -+
       units  | ohm   ohm |
              | ohm   ohm |
              +-         -+

    Z = inv(Y)
    """

    @classmethod
    def generic(cls):
        return cls((('Z_11', 'Z_12'), ('Z_21', 'Z_22')))

    def equation(self):
        return Eq(Matrix(('V1', 'V2')), MatMul(self, Matrix(('I1', 'I2'))),
                  evaluate=False)

    @property
    def Aparams(self):
        det = self.det().expr
        return AMatrix(((self._Z11 / self._Z21, det / self._Z21),
                        (1 / self._Z21, self._Z22 / self._Z21))).simplify()

    @property
    def Bparams(self):
        if not hasattr(self, '_Bparams'):
            det = self.det().expr
            self._Bparams = BMatrix(((self._Z22 / self._Z12, -det / self._Z12),
                                     (-1 / self._Z12, self._Z11 / self._Z12))).simplify()
        return self._Bparams

    @property
    def Hparams(self):
        det = self.det().expr
        return HMatrix(((det / self._Z22, self._Z12 / self._Z22),
                        (-self._Z21 / self._Z22, 1 / self._Z22))).simplify()

    @property
    def Yparams(self):

        det = self.det().expr
        return YMatrix(((self._Z22 / det, -self._Z12 / det),
                        (-self._Z21 / det, self._Z11 / det))).simplify()

    @property
    def Zparams(self):
        # Perhaps we should make a copy?
        return self

    @property
    def Z1oc(self):
        """Open-circuit input impedance"""
        return LaplaceDomainImpedance(self._Z11)

    @property
    def Z2oc(self):
        """Open-circuit input impedance"""
        return LaplaceDomainImpedance(self._Z22)

    @property
    def Vgain12(self):
        """Forward voltage gain"""
        return LaplaceDomainTransferFunction(self._Z21 / self._Z11)

    @property
    def Vgain21(self):
        """Reverse voltage gain"""
        return LaplaceDomainTransferFunction(self._Z12 / self._Z22)

    @property
    def Igain12(self):
        """Forward current gain"""
        return LaplaceDomainTransferFunction(-self._Z21 / self._Z22)

    @property
    def Igain21(self):
        """Reverse current gain"""
        return LaplaceDomainTransferFunction(-self._Z21 / self._Z11)

    @classmethod
    def Lsection(cls, Z1, Z2):
        return cls.Tsection(Z1 + Z2, Z2, Z2, Z2)

    @classmethod
    def Tsection(cls, Z1, Z2, Z3):
        # Note if Z3 is infinity then all elements of Z are infinite.
        # Thus we cannot model a single series R with a Z matrix.
        # A single shunt R works though.
        return cls(((Z1 + Z2, Z2),
                    (Z2, Z2 + Z3)))

    @classmethod
    def Pisection(cls, Z1, Z2, Z3):

        Za, Zb, Zc = DeltaWye(Z1, Z2, Z3)
        return cls.Tsection(Za, Zb, Zc)


# Probably should only inherit a subset of Network
class TwoPort(Network, TwoPortMixin):
    """
    General class for two-port networks.  Two-port networks are
    constrained to have the same current at each port (but flowing in
    opposite directions).  This is called the port condition.
    """

    def __init__(self, *args, **kwargs):

        self.args = args
        self.kwargs = kwargs

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        if hasattr(self, 'tp'):
            return self.tp._net_make(netlist, n1, n2, n3, n4, dir)

        n2, n1, n4, n3 = netlist._make_nodes(n2, n1, n4, n3)

        return 'TP? %s %s %s %s B %s %s %s %s; %s' % (n3, n4, n1, n2,
                                                      netlist._netarg(
                                                          self.B11), netlist._netarg(self.B12),
                                                      netlist._netarg(
                                                          self.B21),
                                                      netlist._netarg(
                                                          self.B22),
                                                      self._opts_str(l=''))

    def _TP_make(self, netlist, n1, n2, n3, n4, kind, *args):

        n2, n1, n4, n3 = netlist._make_nodes(n2, n1, n4, n3)

        args = ' '.join([netlist._netarg(arg) for arg in args])

        s = 'TP? %s %s %s %s %s %s; right, %s' % (n3, n4, n1, n2,
                                                  kind, args,
                                                  self._opts_str(l=''))

        return s

    def _add_elements(self):
        raise ValueError('Cannot generate netlist for two-port objects')

    def netlist(self, layout='horizontal', evalf=None):
        """Create a netlist.

        `layout` can be 'horizontal'.

        `evalf` can be False or an integer specifying the number of
        decimal places used to evaluate floats.
        """

        from .netlistmaker import NetlistMaker
        return NetlistMaker(self, layout=layout, evalf=evalf)()

    def _check_twoport_args(self, args):

        # This is an interim measure until Par2, Ser2, etc. generalised.
        if len(args) != 2:
            raise ValueError('Only two args supported for %s' %
                             self.__class__.__name__)
        for arg1 in args:
            if not isinstance(arg1, TwoPort):
                raise ValueError('%s not a TwoPort' % arg1)

    @property
    def params(self):
        return self._params

    @property
    def sources(self):
        return self._sources

    @property
    def Aparams(self):
        """Return chain parameters"""
        return self._params.Aparams

    @property
    def Bparams(self):
        """Return inverse chain parameters"""
        if not hasattr(self, '_Bparams'):
            self._Bparams = self._params.Bparams
        return self._Bparams

    @property
    def Gparams(self):
        """Return inverse hybrid parameters"""
        return self._params.Gparams

    @property
    def Hparams(self):
        """Return hybrid parameters"""
        return self._params.Hparams

    @property
    def Sparams(self):
        """Return scattering parameters"""
        return self._params.Sparams

    @property
    def Tparams(self):
        """Return scattering transfer parameters"""
        return self._params.Tparams

    @property
    def Yparams(self):
        """Return admittance parameters"""
        return self._params.Yparams

    @property
    def Zparams(self):
        """Return impedance parameters"""
        return self._params.Zparams

    @property
    def I1a(self):
        return -LaplaceDomainCurrent(self.A21 * self.V2b) - self.A22 * self.I2b

    @property
    def V1a(self):
        return -self.A11 * self.V2b - LaplaceDomainVoltage(self.A12 * self.I2b)

    @property
    def I1g(self):
        return LaplaceDomainCurrent(-self.I2b / self._B22)

    @property
    def V2g(self):
        return self.V2b - LaplaceDomainVoltage(self._B21 / self._B22 * self.I2b)

    @property
    def V1h(self):
        return LaplaceDomainVoltage(-self.V2b / self._B11)

    @property
    def I2h(self):
        return LaplaceDomainCurrent(-self.V2b * self._B21 / self._B11) - self.I2b

    @property
    def I1y(self):
        return LaplaceDomainCurrent(-self.V2b / self._B12)

    @property
    def I2y(self):
        return LaplaceDomainCurrent(self.V2b * self._B22 / self._B12) - self.I2b

    @property
    def V1z(self):
        return LaplaceDomainVoltage(-self.I2b / self._B21)

    @property
    def V2z(self):
        return self.V2b - LaplaceDomainVoltage(self.I2b * self._B11 / self._B21)

    @property
    def Yoc(self):
        """Return admittance vector with ports open circuit"""
        return LaplaceDomainAdmittanceMatrix((LaplaceDomainAdmittance(1 / self.Z1oc),
                                              LaplaceDomainAdmittance(1 / self.Z2oc)))

    @property
    def Y1oc(self):
        """Return input admittance with the output port open circuit"""
        return LaplaceDomainAdmittance(1 / self.Z1oc)

    @property
    def Y2oc(self):
        """Return output admittance with the input port open circuit"""
        return LaplaceDomainAdmittance(1 / self.Z2oc)

    @property
    def Ysc(self):
        """Return admittance vector with ports short circuit"""
        return self.Yparams.Ysc

    @property
    def Y1sc(self):
        """Return input admittance with output port short circuit"""
        # Y11, A22 / A12
        return LaplaceDomainAdmittance(self._Y11)

    @property
    def Y2sc(self):
        """Return output admittance with output port short circuit"""
        # Y22, A11 / A12
        return LaplaceDomainAdmittance(self._Y22)

    @property
    def Zoc(self):
        """Return impedance vector with ports open circuit"""
        return LaplaceDomainImpedanceMatrix((self.Z1oc, self.Z2oc))

    @property
    def Z1oc(self):
        """Return input impedance with the output port open circuit"""
        return self.params.Z1oc

    @property
    def Z2oc(self):
        """Return output impedance with the input port open circuit"""
        return self.params.Z2oc

    @property
    def Zsc(self):
        """Return impedance vector with ports short circuit"""
        return LaplaceDomainImpedanceMatrix((self.Z1sc, self.Z2sc))

    @property
    def Z1sc(self):
        """Return input impedance with the output port short circuit"""
        return self.params.Z1sc

    @property
    def Z2sc(self):
        """Return output impedance with the input port short circuit"""
        return self.params.Z2sc

    def Vgain(self, inport=1, outport=2):
        """Return voltage gain for specified ports with internal
        sources zero"""

        # Av  = G21 = 1 / A11 = -det(B) / B22 = Z21 / Z11 =  Y21 / Y22
        # Av' = H12 = 1 / B11 =  |A| / A22 = Z12 / Z22 = -Y12 / Y11

        if inport == outport:
            return LaplaceDomainTransferFunction(1)
        if inport == 1 and outport == 2:
            return self.Vgain12
        if inport == 2 and outport == 1:
            return self.Vgain21
        raise ValueError('bad port values')

    def Igain(self, inport=1, outport=2):
        """Return current gain for specified ports with internal
         sources zero"""

        # Ai  = H21 = -1 / A22 = -det(B) / B11 = -Z21 / Z22 = Y21 / Y11
        # Ai' = G12 =  1 / B22 =  |A| / A11 = -Z12 / Z11 = Y12 / Y22

        if inport == outport:
            return LaplaceDomainTransferFunction(1)
        if inport == 1 and outport == 2:
            return self.Igain12
        if inport == 2 and outport == 1:
            return self.Igain21
        raise ValueError('bad port values')

    @property
    def Vgain12(self):
        """Return V2 / V1 for I2 = 0 (forward voltage gain) with
        internal sources zero

        Av = G21 = 1 / A11 = -det(B) / B22 = Z21 / Z11 =  Y21 / Y22
        """

        return self.params.Vgain12

    @property
    def Vgain21(self):
        """Return V1 / V2 for I1 = 0 (reverse voltage gain) with
        internal sources zero

        """

        return self.params.Vgain21

    @property
    def Vtransfer(self):
        """Return V2 / V1 for I2 = 0 (forward voltage gain) with
        internal sources zero  (see Vgain12)"""

        return self.params.Vgain12

    @property
    def Igain12(self):
        """Return I2 / I1 for V2 = 0 (forward current gain) with
        internal sources zero

        Ai = H21 = -1 / A22 = -det(B) / B11 = -Z21 / Z22 = Y21 / Y11
        """

        return self.params.Igain12

    @property
    def Igain21(self):
        """Return I1 / I2 for V1 = 0 (reverse current gain) with
        internal sources zero

        """

        return self.params.Igain21

    @property
    def Itransfer(self):
        """Return I2 / I1 for V2 = 0 (forward current gain) with
        internal sources zero  (sett Igain12)"""

        return self.params.Igain12

    @property
    def forward_transadmittance(self):
        """Return I2 / V1 for V2 = 0 with internal sources zero."""

        return self.params.forward_transadmittance

    @property
    def reverse_transadmittance(self):
        """Return I1 / V2 for V1 = 0 with internal sources zero."""

        return self.params.reverse_transadmittance

    def Vresponse(self, V, inport=1, outport=2):
        """Return voltage response for specified applied voltage and
        specified ports"""

        if issubclass(V.__class__, OnePort):
            V = V.Voc.laplace()

        p1 = inport - 1
        p2 = outport - 1

        H = self.Zparams[p2, p1] / self.Zparams[p1, p1]
        return LaplaceDomainVoltage(self.Voc[p2] + (V - self.Voc[p1]) * H)

    def Iresponse(self, I, inport=1, outport=2):
        """Return current response for specified applied current and
        specified ports"""

        if issubclass(I.__class__, OnePort):
            I = I.Isc.laplace()

        p1 = inport - 1
        p2 = outport - 1

        Y = self.Yparams
        Isc = self.Isc

        return LaplaceDomainCurrent(Isc[p2] + Y[p2, p1] / Y[p1, p1] * (I - Isc[p1]))

    def Ytrans(self, inport=1, outport=2):
        """Return transadmittance for specified ports with internal
        sources zero"""

        return LaplaceDomainAdmittance(self.Yparams[outport - 1, inport - 1])

    @property
    def Ytrans12(self):
        """Return I2 / V1 for V2 = 0 (forward transadmittance) with
         internal sources zero

         Y21 = -1 / A12 = det(B) / B12
         """

        return LaplaceDomainAdmittance(self._Y21)

    @property
    def Ytransfer(self):
        """Return I2 / V1 for V2 = 0 (forward transadmittance) with
         internal sources zero.  This is an alias for Ytrans12.

         Y21 = -1 / A12 = det(B) / B12
         """

        return self.Ytrans12

    def Ztrans(self, inport=1, outport=2):
        """Return transimpedance for specified ports with internal
        sources zero"""

        return LaplaceDomainImpedance(self.Zparams[outport - 1, inport - 1])

    def Ztrans12(self):
        """Return V2 / I1 for I2 = 0 (forward transimpedance) with
        internal sources zero

        Z21 = 1 / A21 = -det(B) / B21
        """

        return LaplaceDomainImpedance(self._Z21)

    @property
    def Ztransfer(self):
        """Return V2 / I1 for I2 = 0 (forward transimpedance) with
        internal sources zero.  This is an alias for Ztrans12.

        Z21 = 1 / A21 = -det(B) / B21
        """

        return self.Ztrans12

    @property
    def V1oc(self):
        """Return V1 with all ports open-circuited (i.e., I1 = I2 = 0)"""
        return LaplaceDomainVoltage(self.Voc[0])

    @property
    def V2oc(self):
        """Return V2 with all ports open-circuited (i.e., I1 = I2 = 0)"""
        return LaplaceDomainVoltage(self.Voc[1])

    @property
    def I1sc(self):
        """Return I1 with all ports short-circuited, i.e, V1 = V2 = 0"""
        return LaplaceDomainCurrent(self.Isc[0])

    @property
    def I2sc(self):
        """Return I2 with all ports short-circuited, i.e, V1 = V2 = 0"""
        return LaplaceDomainCurrent(self.Isc[1])

    @property
    def Voc(self):
        """Return voltage vector with all ports open-circuited
        (i.e., In = 0)"""
        return LaplaceDomainVoltageMatrix((self.V1z, self.V2z))

    @property
    def Isc(self):
        """Return current vector with all ports short-circuited
        (i.e., V1 = V2 = 0)"""
        return LaplaceDomainCurrentMatrix((self.I1y, self.I2y))

    @property
    def Amodel(self):

        return TwoPortAModel(self.Aparams, V1a=self.V1a, I1a=self.I1a)

    @property
    def Bmodel(self):

        return TwoPortBModel(self.Bparams, V2b=self.V2b, I2b=self.I2b)

    @property
    def Gmodel(self):

        return TwoPortGModel(self.Gparams, I1g=self.I1g, V2g=self.V2g)

    @property
    def Hmodel(self):

        return TwoPortHModel(self.Hparams, V1h=self.V1h, I2h=self.I2h)

    @property
    def Ymodel(self):

        if self.is_shunt:
            warn('Converting a shunt two-port to a Y model is dodgy...')
        return TwoPortYModel(self.Yparams, I1y=self.I1y, I2y=self.I2y)

    @property
    def Zmodel(self):

        if self.is_series:
            warn('Converting a series two-port to a Z model is dodgy...')
        return TwoPortZModel(self.Zparams, V1z=self.V1z, V2z=self.V2z)

    def chain(self, TP):
        """Return the model with, TP, appended (cascade or
        chain connection)"""

        if not issubclass(TP.__class__, TwoPort):
            raise TypeError('Argument not', TwoPort)

        return Chain(self, TP)

    def append(self, TP):
        """Return the model with, TP, appended (this is equivalent to chain)"""

        return self.chain(TP)

    def prepend(self, TP):
        """Return the model with, TP, prepended"""

        return TP.chain(self)

    def cascade(self, TP):
        """Return the model with, TP, appended (this is equivalent to chain)"""

        return self.chain(TP)

    def series(self, TP, port=None):
        """Return the model with, TP, in series.

         In general, this is tricky to ensure that the port condition
         is valid.  The common ground connection of the first two-port
         shorts out the top of the T of the second two-port.
         """

        if issubclass(TP.__class__, OnePort):
            raise NotImplementedError('TODO')

        warn('Do you mean chain?  The result of a series combination'
             ' of two two-ports may be dodgy')

        return Ser2(self, TP)

    def terminate(self, OP, port=2):
        """Connect one-port in parallel to specified port and return
        a Thevenin (one-port) object"""

        if port == 1:
            return self.source(OP)
        if port == 2:
            return self.load(OP)
        raise ValueError('Invalid port ' + port)

    def parallel(self, TP, port=None):
        """Return the model with, TP, in parallel"""

        if issubclass(TP.__class__, OnePort):
            raise NotImplementedError('TODO')

        return Par2(self, TP)

    def hybrid(self, TP, port=None):
        """Return the model with, TP, in hybrid connection (series
        input, parallel output)"""

        if issubclass(TP.__class__, OnePort):
            raise NotImplementedError('TODO')

        return Hybrid2(self, TP)

    def inverse_hybrid(self, TP, port=None):
        """Return the model with, TP, in inverse hybrid connection
        (parallel input, series output)"""

        if issubclass(TP.__class__, OnePort):
            raise NotImplementedError('TODO')

        return InverseHybrid2(self, TP)

    # Other operations: swapping the input terminals negates the A matrix.
    # switching ports.

    def bridge(self, OP):
        """Bridge the ports with a one-port element"""

        if not issubclass(OP.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        # FIXME
        return self.parallel(Series(OP))

    def load(self, OP):
        """Apply a one-port load and return a Thevenin (one-port) object"""

        if not issubclass(OP.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        foo = self.chain(Shunt(OP))
        return V(foo.V1oc) + Z(foo.Z1oc)

    def source(self, OP):
        """Apply a one-port source and return a Thevenin (one-port) object"""

        if not issubclass(OP.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        foo = Shunt(OP).chain(self)
        return V(foo.V2oc) + Z(foo.Z2oc)

    def short_circuit(self, port=2):
        """Apply a short-circuit to specified port and return a
        one-port object"""

        p = port - 1
        Yval = self.Yparams[1 - p, 1 - p]
        Ival = self.Isc[1 - p]

        return (I(Ival) | Y(Yval)).simplify()

    def open_circuit(self, port=2):
        """Apply a open-circuit to specified port and return a
        one-port object"""

        p = port - 1
        Zval = self.Zparams[1 - p, 1 - p]
        Vval = self.Voc[1 - p]

        return (V(Vval) + Z(Zval)).simplify()

    def simplify(self):

        if self.Bparams == sym.eye(2):
            # Have a pair of wires... perhaps could simplify
            # to an LSection comprised of a V and I but
            # may have a weird voltage expression.
            pass
        return self

    def __add__(self, OP):
        """Series combination"""

        return self.series(OP)

    def __or__(self, OP):
        """Parallel combination"""

        return self.parallel(OP)

    def __mul__(self, OP):
        """Chained combination"""

        return self.chain(OP)

    def __eq__(self, x):

        if x.__class__ != self.__class__:
            return False

        return self.params == x.params and self.sources == x.sources

    def equation(self):
        """Return equation describing model."""

        input = Matrix(self.input).expr
        output = Matrix(self.output).expr
        params = self.params.expr
        sources = self.sources.expr

        if sources[0] == 0 and sources[1] == 0:
            return expr(sym.Eq(output, sym.MatMul(params, input),
                               evaluate=False))

        return expr(sym.Eq(output,
                           sym.MatAdd(sym.MatMul(params, input),
                                      sources), evaluate=False))


class TwoPortBModel(TwoPort):
    """
    ::
            +-------------------+    +------+
     I1     |                   | I2'| -  + |          I2
    -->-----+                   +-<--+  V2b +----+-----<--
            |    two-port       |    |      |    |
    +       |    network        | +  +------+    |       +
    V1      |    without        | V2'        +---+---+  V2
    -       |    sources        | -          |   |   |   -
            |    represented    |               I2b  |
            |    by B matrix    |            |   v   |
            |                   |            +---+---+
            |                   |                |
    --------+                   +----------------+--------
            |                   |
            +-------------------+

    +-   +     +-        -+   +-  -+     +-   -+
    | V2 |  =  | B11  B12 |   | V1 |  +  | V2b |
    |-I2 |     | B21  B22 |   | I1 |     | I2b |
    +-  -+     +-        -+   +-  -+     +-   -+


    +-    +     +-        -+   +-  -+
    | V2' |  =  | B11  B12 |   | V1 |
    |-I2' |     | B21  B22 |   | I1 |
    +-   -+     +-        -+   +-  -+

    +-    +     +-  -+    +-   -+
    | V2' |  =  | V2 | -  | V2b |
    | I2' |     | I2'|    | I2b |
    +-  - +     +-  -+    +-   -+

    +-         +     +-        -+   +-  -+
    | V2 - V2b |  =  | B11  B12 |   | V1 |
    |-I2 - I2b |     | B21  B22 |   | I1 |
    +-        -+     +-        -+   +-  -+

    +-   +     +-        -+   +-    +
    | V1 |  =  | A11  A12 |   | V2' |
    | I1 |     | A21  A22 |   |-I2' |
    +-  -+     +-        -+   +-   -+

    +-   +     +-        -+   +-       -+
    | V1 |  =  | A11  A12 |   | V2 - V2b |
    | I1 |     | A21  A22 |   |-I2 + I2b |
    +-  -+     +-        -+   +-        -+

    """

    model = 'B'
    input = ('V1', '-I1')
    output = ('V2', 'I2')
    offset = ('V2b', 'I2b')

    def __init__(self, B11=None, B12=None, B21=None, B22=None,
                 V2b=None, I2b=None, **kwargs):

        if B11 is not None and B12 is None and B21 is None and B22 is None:
            B = B11
        else:
            B11 = 'B11' if B11 is None else B11
            B12 = 'B12' if B12 is None else B12
            B21 = 'B21' if B21 is None else B21
            B22 = 'B22' if B22 is None else B22
            B = BMatrix(((B11, B12), (B21, B22)))

        if V2b is None:
            V2b = LaplaceDomainVoltage(0)
        if I2b is None:
            I2b = LaplaceDomainCurrent(0)

        if issubclass(B.__class__, TwoPortBModel):
            B, V2b, I2b = B._params, B._V2b, B._I2b

        if not isinstance(B, BMatrix):
            raise ValueError('B not BMatrix')

        V2b = LaplaceDomainVoltage(V2b)
        I2b = LaplaceDomainCurrent(I2b)

        super(TwoPortBModel, self).__init__(
            B[0, 0], B[0, 1], B[1, 0], B[1, 1], V2b, I2b, **kwargs)
        self._params = B
        self._V2b = V2b
        self._I2b = I2b
        self._sources = Vector(V2b, I2b)

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        return self._TP_make(netlist, n1, n2, n3, n4, 'B',
                             self.B11, self.B12, self.B21, self.B22,
                             self.V2b, self.I2b)

    @property
    def Bparams(self):
        """Return chain matrix"""
        return self._params

    @property
    def I2b(self):
        return self._I2b

    @property
    def V2b(self):
        return self._V2b

    @property
    def V1h(self):
        return LaplaceDomainVoltage(-self.V2b / self._B11)

    @property
    def I2h(self):
        return LaplaceDomainCurrent(-self.V2b * self._B21 / self._B11) - self.I2b

    @property
    def I1y(self):
        return LaplaceDomainCurrent(-self.V2b / self._B12)

    @property
    def I2y(self):
        return LaplaceDomainCurrent(self.V2b * self._B22 / self._B12) - self.I2b

    @property
    def V1z(self):
        return LaplaceDomainVoltage(-self.I2b / self._B21)

    @property
    def V2z(self):
        return self.V2b - LaplaceDomainVoltage(self.I2b * self._B11 / self._B21)


class TwoPortAModel(TwoPort):
    """
    """

    model = 'A'
    input = ('V2', '-I2')
    output = ('V1', 'I1')
    offset = ('V1a', 'I1a')

    def __init__(self, A11=None, A12=None, A21=None, A22=None,
                 V1a=None, I1a=None, **kwargs):

        if A11 is not None and A12 is None and A21 is None and A22 is None:
            A = A11
        else:
            A11 = 'A11' if A11 is None else A11
            A12 = 'A12' if A12 is None else A12
            A21 = 'A21' if A21 is None else A21
            A22 = 'A22' if A22 is None else A22
            A = AMatrix(((A11, A12), (A21, A22)))

        if V1a is None:
            V1a = LaplaceDomainVoltage(0)
        if I1a is None:
            I1a = LaplaceDomainCurrent(0)

        if issubclass(A.__class__, TwoPortAModel):
            A, V1a, I1a = A._params, A._V1a, A._I1a

        if not isinstance(A, AMatrix):
            raise ValueError('A not AMatrix')

        V1a = LaplaceDomainVoltage(V1a)
        I1a = LaplaceDomainCurrent(I1a)

        super(TwoPortAModel, self).__init__(
            A[0, 0], A[0, 1], A[1, 0], A[1, 1], V1a, I1a, **kwargs)
        self._params = A
        self._V1a = V1a
        self._I1a = I1a
        self._sources = Vector(V1a, I1a)

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        return self._TP_make(netlist, n1, n2, n3, n4, 'A',
                             self.A11, self.A12, self.A21, self.A22,
                             self.V1a, self.I1a)

    @property
    def Aparams(self):
        """Return chain matrix"""
        return self._params

    @property
    def I1a(self):
        return self._I1a

    @property
    def V1a(self):
        return self._V1a

    @property
    def I2b(self):
        # Avoid matrix inverse calculating Bparams
        if self.V1a == 0 and self.I1a == 0:
            return self.I1a

        return -LaplaceDomainCurrent(self.B21 * self.V1a) - self.B22 * self.I1a

    @property
    def V2b(self):
        # Avoid matrix inverse calculating Bparams
        if self.V1a == 0 and self.I1a == 0:
            return self.V1a

        return -self.B11 * self.V1a - LaplaceDomainVoltage(self.B12 * self.I1a)


class TwoPortGModel(TwoPort):
    """
    """

    model = 'G'
    input = ('V1', 'I2')
    output = ('I1', 'V2')
    offset = ('I1g', 'V2g')

    def __init__(self, G11=None, G12=None, G21=None, G22=None,
                 I1g=None, V2g=None, **kwargs):

        if G11 is not None and G12 is None and G21 is None and G22 is None:
            G = G11
        else:
            G11 = 'G11' if G11 is None else G11
            G12 = 'G12' if G12 is None else G12
            G21 = 'G21' if G21 is None else G21
            G22 = 'G22' if G22 is None else G22
            G = GMatrix(((G11, G12), (G21, G22)))

        if I1g is None:
            I1g = LaplaceDomainCurrent(0)
        if V2g is None:
            V2g = LaplaceDomainVoltage(0)

        if issubclass(G.__class__, TwoPortGModel):
            G, I1g, V2g = G._params, G._I1g, G._V2g

        if not isinstance(G, GMatrix):
            raise ValueError('G not GMatrix')

        I1g = LaplaceDomainCurrent(I1g)
        V2g = LaplaceDomainVoltage(V2g)

        super(TwoPortGModel, self).__init__(
            G[0, 0], G[0, 1], G[1, 0], G[1, 1], I1g, V2g, **kwargs)
        self._params = G
        self._I1g = I1g
        self._V2g = V2g
        self._sources = Vector(I1g, V2g)

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        return self._TP_make(netlist, n1, n2, n3, n4, 'G',
                             self.G11, self.G12, self.G21, self.G22,
                             self.I1g, self.V2g)

    @property
    def Gparams(self):
        """Return hybrid matrix"""
        return self._params

    @property
    def V2b(self):
        """Return V2b"""

        # return self._V2g - LaplaceDomainVoltage(self._I1g / self.Gparams._G12)
        return self._V2g - LaplaceDomainVoltage(self._I1g * self.B21)

    @property
    def I2b(self):
        """Return I2b"""

        return -self.B22 * self._I1g

    @property
    def I1g(self):
        return self._I1g

    @property
    def V2g(self):
        return self._V2g


class TwoPortHModel(TwoPort):
    """
    ::
         +------+   +-------------------+
     I1  | +  - |   |                   | I2'          I2
    -->--+  V1h +---+                   +-<-------+-----<--
         |      |   |    two-port       |         |
    +    +------+   |    network        | +       |       +
    V1              |    without        | V2' ---+---+  V2
    -               |    sources        | -  |   |   |   -
                    |    represented    |    |  I2h  |
                    |    by H matrix    |    |   v   |
                    |                   |    +---+---+
                    |                   |        |
    ----------------+                   +--------+--------
                    |                   |
                    +-------------------+


    +-   +     +-        -+   +-  -+     +-   -+
    | V1 |  =  | H11  H12 |   | I1 |  +  | V1h |
    | I2 |     | H21  H22 |   | V2 |     | I2h |
    +-  -+     +-        -+   +-  -+     +-   -+
    """

    model = 'H'
    input = ('I1', 'V2')
    output = ('V1', 'I2')
    offset = ('V1h', 'I2h')

    def __init__(self, H11=None, H12=None, H21=None, H22=None,
                 V1h=None, I2h=None, **kwargs):

        if H11 is not None and H12 is None and H21 is None and H22 is None:
            H = H11
        else:
            H11 = 'H11' if H11 is None else H11
            H12 = 'H12' if H12 is None else H12
            H21 = 'H21' if H21 is None else H21
            H22 = 'H22' if H22 is None else H22
            H = HMatrix(((H11, H12), (H21, H22)))

        if V1h is None:
            V1h = LaplaceDomainVoltage(0)
        if I2h is None:
            I2h = LaplaceDomainCurrent(0)

        if issubclass(H.__class__, TwoPortHModel):
            H, V1h, I2h = H._params, H._V1h, H._I2h

        if not isinstance(H, HMatrix):
            raise ValueError('H not HMatrix')

        V1h = LaplaceDomainVoltage(V1h)
        I2h = LaplaceDomainCurrent(I2h)

        super(TwoPortHModel, self).__init__(
            H[0, 0], H[0, 1], H[1, 0], H[1, 1], V1h, I2h, **kwargs)
        self._params = H
        self._V1h = V1h
        self._I2h = I2h
        self._sources = Vector(V1h, I2h)

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        return self._TP_make(netlist, n1, n2, n3, n4, 'H',
                             self.H11, self.H12, self.H21, self.H22,
                             self.V1h, self.I2h)

    @property
    def Hparams(self):
        """Return hybrid matrix"""
        return self._params

    @property
    def V2b(self):
        """Return V2b"""

        return LaplaceDomainVoltage(-self.V1h / self.Hparams._H12)

    @property
    def I2b(self):
        """Return I2b"""

        return LaplaceDomainCurrent(-self.H22 / self.H12 * self.V1h) - self.I2h

    @property
    def V1h(self):
        return self._V1h

    @property
    def I2h(self):
        return self._I2h


class TwoPortYModel(TwoPort):
    """
    ::
                     +-------------------+
     I1              |                   | I2'           I2
    -->----+---------+                   +-<-------+-----<--
           |         |    two-port       |         |
    +      |       + |    network        | +       +       +
    V1 +---+---+  V1'|    without        | V2' +---+---+  V2
    -  |   |   |   - |    sources        | -   |   |   |   -
       |  I1y  |     |    represented    |     |  I2y  |
       |   v   |     |    by Y matrix    |     |   v   |
       +---+---+     |                   |     +---+---+
           |         |                   |         |
    -------+---------+                   +---------+--------
                     |                   |
                     +-------------------+

    +-   +     +-        -+   +-  -+     +-   -+
    | I1 |  =  | Y11  Y12 |   | V1 |  +  | I1y |
    | I2 |     | Y21  Y22 |   | V2 |     | I2y |
    +-  -+     +-        -+   +-  -+     +-   -+

    Ymn = Im / Vn for Vm = 0
    """

    model = 'Y'
    input = ('V1', 'V2')
    output = ('I1', 'I2')
    offset = ('I1y', 'I2y')

    def __init__(self, Y11=None, Y12=None, Y21=None, Y22=None,
                 I1y=None, I2y=None, **kwargs):

        if Y11 is not None and Y12 is None and Y21 is None and Y22 is None:
            Y = Y11
        else:
            Y11 = 'Y11' if Y11 is None else Y11
            Y12 = 'Y12' if Y12 is None else Y12
            Y21 = 'Y21' if Y21 is None else Y21
            Y22 = 'Y22' if Y22 is None else Y22
            Y = YMatrix(((Y11, Y12), (Y21, Y22)))

        if I1y is None:
            I1y = LaplaceDomainCurrent(0)
        if I2y is None:
            I2y = LaplaceDomainCurrent(0)

        if issubclass(Y.__class__, TwoPortYModel):
            Y, I1y, I2y = Y._params, Y._I1y, Y._I2y

        if not isinstance(Y, YMatrix):
            raise ValueError('Y not YMatrix')

        I1y = LaplaceDomainCurrent(I1y)
        I2y = LaplaceDomainCurrent(I2y)

        super(TwoPortYModel, self).__init__(
            Y[0, 0], Y[0, 1], Y[1, 0], Y[1, 1], I1y, I2y, **kwargs)
        self._params = Y
        self._I1y = I1y
        self._I2y = I2y
        self._sources = Vector(I1y, I2y)

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        return self._TP_make(netlist, n1, n2, n3, n4, 'Y',
                             self.Y11, self.Y12, self.Y21, self.Y22,
                             self.I1y, self.I2y)

    @property
    def Yparams(self):
        """Return admittance matrix"""
        return self._params

    @property
    def I2b(self):
        return LaplaceDomainCurrent(self.I1y * self._Y22 / self._Y12) - self.I2y

    @property
    def V2b(self):
        return LaplaceDomainVoltage(-self.I1y / self._Y12)

    @property
    def I1y(self):
        return self._I1y

    @property
    def I2y(self):
        return self._I2y


class TwoPortZModel(TwoPort):
    """
    ::
         +------+    +-------------------+    +------+
    I1   | +  - | I1'|                   | I2'| -  + |  I2
    -->--+  V1z +-->-+                   +-<--+  V2z +--<--
         |      |    |    two-port       |    |      |
    +    +------+  + |    network        | +  +------+    +
    V1            V1'|    without        | V2'           V2
    -              - |    sources        | -              -
                     |    represented    |
                     |    by Z matrix    |
                     |                   |
                     |                   |
    -----------------+                   +-----------------
                     |                   |
                     +-------------------+

    +-   +     +-        -+   +-  -+     +-   -+
    | V1 |  =  | Z11  Z12 |   | I1 |  +  | V1z |
    | V2 |     | Z21  Z22 |   | I2 |     | V2z |
    +-  -+     +-        -+   +-  -+     +-   -+

    """

    model = 'Z'
    input = ('I1', 'I2')
    output = ('V1', 'V2')
    offset = ('V1z', 'V2z')

    def __init__(self, Z11=None, Z12=None, Z21=None, Z22=None,
                 V1z=None, V2z=None, **kwargs):

        if Z11 is not None and Z12 is None and Z21 is None and Z22 is None:
            Z = Z11
        else:
            Z11 = 'Z11' if Z11 is None else Z11
            Z12 = 'Z12' if Z12 is None else Z12
            Z21 = 'Z21' if Z21 is None else Z21
            Z22 = 'Z22' if Z22 is None else Z22
            Z = ZMatrix(((Z11, Z12), (Z21, Z22)))

        if V1z is None:
            V1z = LaplaceDomainVoltage(0)
        if V2z is None:
            V2z = LaplaceDomainVoltage(0)

        if issubclass(Z.__class__, TwoPortZModel):
            Z, V1z, V2z = Z._params, Z._V1z, Z._V2z

        if not isinstance(Z, ZMatrix):
            raise ValueError('Z not ZMatrix')

        V1z = LaplaceDomainVoltage(V1z)
        V2z = LaplaceDomainVoltage(V2z)

        super(TwoPortZModel, self).__init__(
            Z[0, 0], Z[0, 1], Z[1, 0], Z[1, 1], V1z, V2z, **kwargs)
        self._params = Z
        self._V1z = V1z
        self._V2z = V2z
        self._sources = Vector(V1z, V2z)

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        return self._TP_make(netlist, n1, n2, n3, n4, 'Z',
                             self.Z11, self.Z12, self.Z21, self.Z22,
                             self.V1z, self.V2z)

    @property
    def Zparams(self):
        """Return impedance matrix"""
        return self._params

    @property
    def I2b(self):
        return LaplaceDomainCurrent(self.V1z / self._Z12)

    @property
    def V2b(self):
        return self.V2z - LaplaceDomainVoltage(self.V1z * self._Z22 / self._Z12)

    @property
    def I1y(self):

        Zdet = self.Zparams.det().expr
        return LaplaceDomainCurrent(-self.V1z * self._Z22 / Zdet + self.V2z * self._Z12 / Zdet)

    @property
    def I2y(self):

        Zdet = self.Zparams.det().expr
        return LaplaceDomainCurrent(self.V1z * self._Z21 / Zdet - self.V2z * self._Z11 / Zdet)

    @property
    def V1z(self):
        return self._V1z

    @property
    def V2z(self):
        return self._V2z


class TPA(TwoPortAModel):
    """A-parameter two-port network."""

    pass


class TPB(TwoPortBModel):
    """B-parameter two-port network."""

    pass


class TPG(TwoPortGModel):
    """G-parameter two-port network."""

    pass


class TPH(TwoPortHModel):
    """H-parameter two-port network."""

    pass


class TPY(TwoPortYModel):
    """Y-parameter two-port network."""

    pass


class TPZ(TwoPortZModel):
    """Z-parameter two-port network."""

    pass


class TP(TPB):
    """A generic two-port network."""

    pass


class Chain(TwoPortBModel):
    """Connect two-port networks in a chain (aka cascade)"""

    def __init__(self, *args):

        self._check_twoport_args(args)

        # FIXME for non-invertible Aparams.

        arg1 = args[-1]
        B = arg1.Bparams

        foo = Vector(arg1.V2b, arg1.I2b)

        for arg in reversed(args[0:-1]):

            foo += B * Vector(arg.V2b, arg.I2b)
            B = B * arg.Bparams

        super(Chain, self).__init__(B, V2b=LaplaceDomainVoltage(foo[0, 0]),
                                    I2b=LaplaceDomainCurrent(foo[1, 0]))
        self.args = args

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        n2, n1, n4, n3, n6, n5 = netlist._make_nodes(
            n2, n1, n4, n3, None, None)

        nets = []
        nets.append(self.args[0]._net_make(netlist, n1, n2, n5, n6))
        nets.append(self.args[1]._net_make(netlist, n5, n6, n3, n4))
        return '\n'.join(nets)

    def simplify(self):

        if isinstance(self.args[0], Shunt) and isinstance(self.args[1], Shunt):
            return Shunt(
                (self.args[0].args[0] | self.args[1].args[0]).simplify())

        if isinstance(self.args[0], Series) and isinstance(
                self.args[1], Series):
            return Series(
                (self.args[0].args[0] + self.args[1].args[0]).simplify())

        return self


class Par2(TwoPortYModel):
    """Connect two-port networks in parallel"""

    def __init__(self, *args):

        self._check_twoport_args(args)

        # This will fail with a Shunt as an argument since it does
        # not have a valid Y model.
        # We can special case this.
        if isinstance(args[0], Shunt) or isinstance(args[1], Shunt):
            warn('A parallel Shunt not properly handled')

        arg = args[0]
        I1y = arg.I1y
        I2y = arg.I2y
        Y = arg.Yparams

        for arg in args[1:]:
            I1y += arg.I1y
            I2y += arg.I2y
            Y += arg.Yparams

        super(Par2, self).__init__(Y, I1y=I1y, I2y=I2y)
        self.args = args

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        n2, n1, n4, n3, n6, n5, n8, n7 = netlist._make_nodes(
            n2, n1, n4, n3, *([None] * 4))
        n10, n9, n12, n11, n14, n13, n16, n15, n18, n17, n20, n19 = netlist._make_nodes(
            *([None] * 12))

        nets = []
        nets.append(self.args[0]._net_make(netlist, n5, n6, n7, n8))
        nets.append(self.args[1]._net_make(netlist, n9, n10, n11, n12))
        nets.append('W %s %s; right=0.75' % (n1, n13))
        nets.append('W %s %s; right=0.25' % (n13, n5))
        nets.append('W %s %s; right=0.75' % (n7, n14))
        nets.append('W %s %s; right=0.25' % (n14, n3))
        nets.append('W %s %s; right=0.75' % (n15, n6))
        nets.append('W %s %s; right=0.25' % (n8, n16))
        nets.append('W %s %s; right=0.25' % (n17, n9))
        nets.append('W %s %s; right=0.75' % (n11, n18))
        nets.append('W %s %s; right=0.25' % (n2, n19))
        nets.append('W %s %s; right=0.75' % (n19, n10))
        nets.append('W %s %s; right=0.25' % (n12, n20))
        nets.append('W %s %s; right=0.75' % (n20, n4))
        nets.append('W %s %s; down=1.75' % (n15, n19))
        nets.append('W %s %s; down=1.75' % (n13, n17))
        nets.append('W %s %s; down=1.75' % (n16, n20))
        nets.append('W %s %s; down=1.75' % (n14, n18))
        nets.append('O %s %s; down' % (n1, n2))
        nets.append('O %s %s; down' % (n3, n4))
        return '\n'.join(nets)

    def simplify(self):

        if isinstance(self.args[0], Shunt) and isinstance(self.args[1], Shunt):
            return Shunt(
                (self.args[0].args[0] | self.args[1].args[0]).simplify())

        if isinstance(self.args[0], Series) and isinstance(
                self.args[1], Series):
            return Series(
                (self.args[0].args[0] | self.args[1].args[0]).simplify())

        return self


class Ser2(TwoPortZModel):
    """Connect two-port networks in series (note this is unusual and can
    break the port condition)"""

    def __init__(self, *args):

        self._check_twoport_args(args)

        # Need to be more rigorous.
        if isinstance(args[1], (Series, LSection, TSection)):
            warn('This can violate the port condition')

        arg = args[0]
        V1z = arg.V1z
        V2z = arg.V2z
        Z = arg.Zparams

        for arg in args[1:]:
            V1z += arg.V1z
            V2z += arg.V2z
            Z += arg.Zparams

        super(Ser2, self).__init__(Z, V1z=V1z, V2z=V2z)
        self.args = args

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        n2, n1, n4, n3 = netlist._make_nodes(n2, n1, n4, n3)
        n6, n5, n8, n7, n10, n9, n12, n11 = netlist._make_nodes(*([None] * 8))

        nets = []
        nets.append(self.args[0]._net_make(netlist, n5, n6, n7, n8))
        nets.append(self.args[1]._net_make(netlist, n9, n10, n11, n12))
        nets.append('W %s %s; right=0.75' % (n1, n5))
        nets.append('W %s %s; right=0.75' % (n7, n3))
        nets.append('W %s %s; right=0.75' % (n2, n10))
        nets.append('W %s %s; right=0.75' % (n12, n4))
        nets.append('W %s %s; down=0.75' % (n6, n9))
        nets.append('W %s %s; down=0.75' % (n8, n11))
        nets.append('O %s %s; down' % (n1, n2))
        nets.append('O %s %s; down' % (n3, n4))
        return '\n'.join(nets)

    def simplify(self):

        if isinstance(self.args[0], Shunt) and isinstance(self.args[1], Shunt):
            return Shunt(
                (self.args[0].args[0] + self.args[1].args[0]).simplify())

        return self


class Hybrid2(TwoPortHModel):
    """Connect two-port networks in hybrid configuration (inputs in
    series, outputs in parallel)"""

    def __init__(self, *args):

        self._check_twoport_args(args)

        arg = args[0]
        V1h = arg.V1h
        I2h = arg.I2h
        H = arg.Hparams

        for arg in args[1:]:
            V1h += arg.V1h
            I2h += arg.I2h
            H += arg.Hparams

        super(Hybrid2, self).__init__(H, V1h=V1h, I2h=I2h)
        self.args = args

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        n2, n1, n4, n3, n6, n5, n8, n7 = netlist._make_nodes(
            n2, n1, n4, n3, *([None] * 4))
        n10, n9, n12, n11, n14, n13, n16, n15, n18, n17 = netlist._make_nodes(
            *([None] * 10))

        nets = []
        nets.append(self.args[0]._net_make(netlist, n5, n6, n7, n8))
        nets.append(self.args[1]._net_make(netlist, n9, n10, n11, n12))
        nets.append('W %s %s; right=0.75' % (n1, n5))
        nets.append('W %s %s; right=0.75' % (n7, n13))
        nets.append('W %s %s; right=0.25' % (n13, n3))
        nets.append('W %s %s; right=0.25' % (n14, n6))
        nets.append('W %s %s; right=0.25' % (n8, n15))
        nets.append('W %s %s; right=0.25' % (n16, n9))
        nets.append('W %s %s; right=0.75' % (n11, n17))
        nets.append('W %s %s; right=0.75' % (n2, n10))
        nets.append('W %s %s; right=0.25' % (n12, n18))
        nets.append('W %s %s; right=0.75' % (n18, n4))
        nets.append('W %s %s; down=0.75' % (n14, n16))
        nets.append('W %s %s; down=1' % (n15, n18))
        nets.append('W %s %s; down=1' % (n13, n17))
        nets.append('O %s %s; down' % (n1, n2))
        nets.append('O %s %s; down' % (n3, n4))
        return '\n'.join(nets)


class InverseHybrid2(TwoPortGModel):
    """Connect two-port networks in inverse hybrid configuration (outputs in
    series, inputs in parallel)"""

    def __init__(self, *args):

        self._check_twoport_args(args)

        arg = args[0]
        I1g = arg.I1g
        V2g = arg.V2g
        G = arg.Gparams

        for arg in args[1:]:
            I1g += arg.I1g
            V2g += arg.V2g
            G += arg.Gparams

        super(InverseHybrid2, self).__init__(G, I1g=I1g, V2g=V2g)
        self.args = args

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        n2, n1, n4, n3, n6, n5, n8, n7 = netlist._make_nodes(
            n2, n1, n4, n3, *([None] * 4))
        n10, n9, n12, n11, n14, n13, n16, n15, n18, n17 = netlist._make_nodes(
            *([None] * 10))

        nets = []
        nets.append(self.args[0]._net_make(netlist, n5, n6, n7, n8))
        nets.append(self.args[1]._net_make(netlist, n9, n10, n11, n12))
        nets.append('W %s %s; right=0.75' % (n1, n13))
        nets.append('W %s %s; right=0.75' % (n7, n3))
        nets.append('W %s %s; right=0.25' % (n13, n5))
        nets.append('W %s %s; right=0.75' % (n14, n6))
        nets.append('W %s %s; right=0.25' % (n8, n15))
        nets.append('W %s %s; right=0.25' % (n16, n9))
        nets.append('W %s %s; right=0.25' % (n11, n17))
        nets.append('W %s %s; right=0.25' % (n2, n18))
        nets.append('W %s %s; right=0.75' % (n18, n10))
        nets.append('W %s %s; right=0.75' % (n12, n4))
        nets.append('W %s %s; down=1' % (n13, n16))
        nets.append('W %s %s; down=1' % (n14, n18))
        nets.append('W %s %s; down=0.75' % (n15, n17))
        nets.append('O %s %s; down' % (n1, n2))
        nets.append('O %s %s; down' % (n3, n4))
        return '\n'.join(nets)


class Series(TwoPortBModel):
    """
    Two-port comprising a single one-port in series configuration
    ::

           +---------+
         --+   OP    +---
           +---------+

         ----------------

    Note, this has a singular Y matrix.
    """

    def __init__(self, OP):

        _check_oneport_args((OP, ))

        super(Series, self).__init__(BMatrix.Zseries(OP.Z.laplace()),
                                     V2b=LaplaceDomainVoltage(
                                         OP.Voc.laplace()),
                                     I2b=LaplaceDomainCurrent(0))
        self.OP = OP
        self.args = (OP, )

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        n2, n1, n4, n3 = netlist._make_nodes(n2, n1, n4, n3)

        nets = []
        nets.append(self.args[0]._net_make(netlist, n1, n3, dir='right'))
        nets.append('W %s %s; right' % (n2, n4))
        nets.append('O %s %s; down' % (n1, n2))
        nets.append('O %s %s; down' % (n3, n4))
        return '\n'.join(nets)


class Shunt(TwoPortBModel):
    """
    Two-port comprising a single one-port in shunt configuration
    ::

         -----+----
              |
            +-+-+
            |   |
            |OP |
            |   |
            +-+-+
              |
         -----+----

    Note, this has a singular Z matrix.
    """

    def __init__(self, OP):

        _check_oneport_args((OP, ))
        super(Shunt, self).__init__(BMatrix.Yshunt(OP.Y.laplace()),
                                    V2b=LaplaceDomainVoltage(0),
                                    I2b=LaplaceDomainCurrent(OP.Isc.laplace()))
        self.OP = OP
        self.args = (OP, )

    def _net_make(self, netlist, n1=None, n2=None, n3=None, n4=None,
                  dir='right'):

        n2, n1, n4, n3, n6, n5 = netlist._make_nodes(
            n2, n1, n4, n3, None, None)

        nets = []
        nets.append(self.args[0]._net_make(netlist, n5, n6, dir='down'))
        nets.append('W %s %s; right=0.5' % (n1, n5))
        nets.append('W %s %s; right=0.5' % (n5, n3))
        nets.append('W %s %s; right=0.5' % (n2, n6))
        nets.append('W %s %s; right=0.5' % (n6, n4))
        return '\n'.join(nets)


class Transformer(TwoPortBModel):
    """Transformer voltage gain alpha, current gain 1 / alpha.
    Note, alpha = 1 / a where a is the turns ratio defined as the
    number of primary turns to the number of secondary turns, a = N_2 / N_1.

    Unlike with IdealTransformer, the parameter alpha can be a function of s.
    """

    def __init__(self, alpha=1):

        super(Transformer, self).__init__(BMatrix.transformer(alpha))
        self.alpha = expr(alpha)
        self.args = (alpha, )


class IdealTransformer(TwoPortBModel):
    """Ideal transformer voltage gain alpha, current gain 1 / alpha.
    Note, alpha = 1 / a where a is the turns ratio defined as the
    number of primary turns to the number of secondary turns, a = N_2 / N_1.

    alpha must be a constant, otherwise use Transformer if the
    parameter alpha can be a function of s.
    """

    def __init__(self, alpha=1):

        super(IdealTransformer, self).__init__(BMatrix.transformer(alpha))
        self.alpha = ConstantDomainExpression(alpha)
        self.args = (alpha, )


class TF(Transformer):
    pass


class IdealGyrator(TwoPortBModel):
    """Ideal gyrator with gyration resistance R.

    A gyrator converts a voltage to current and a current to voltage.
    Cascaded gyrators act like a transformer"""

    def __init__(self, R=1):

        super(IdealGyrator, self).__init__(BMatrix.gyrator(R))
        self.R = ConstantDomainExpression(R)
        self.args = (R, )


class VoltageFollower(TwoPortBModel):
    """Voltage follower"""

    def __init__(self):

        super(VoltageFollower, self).__init__(BMatrix.voltage_amplifier(1))
        self.args = ()


class VoltageAmplifier(TwoPortBModel):
    """Voltage amplifier"""

    def __init__(self, Av=1, Af=0, Yin=0, Zout=0):

        Av = LaplaceDomainExpression(Av)
        Af = LaplaceDomainExpression(Af)
        Yin = LaplaceDomainExpression(Yin)
        Zout = LaplaceDomainExpression(Zout)

        super(VoltageAmplifier, self).__init__(
            BMatrix.voltage_amplifier(Av, Af, Yin, Zout))
        self.args = (Av, Af, Yin, Zout)


class IdealVoltageAmplifier(TwoPortBModel):
    """Ideal voltage amplifier"""

    def __init__(self, Av=1):

        Av = LaplaceDomainExpression(Av)
        super(IdealVoltageAmplifier, self).__init__(
            BMatrix.voltage_amplifier(Av))
        self.args = (Av, )


class IdealDelay(TwoPortBModel):
    """Ideal buffered delay"""

    def __init__(self, delay=0):

        delay = ConstantDomainExpression(delay)
        super(IdealDelay, self).__init__(
            BMatrix.voltage_amplifier(exp(-s * delay)))
        self.args = (delay, )


class IdealVoltageDifferentiator(TwoPortBModel):
    """Voltage differentiator"""

    def __init__(self, Av=1):

        Av = LaplaceDomainExpression(Av)
        super(IdealVoltageDifferentiator, self).__init__(
            BMatrix.voltage_differentiator(Av))
        self.args = (Av, )


class IdealVoltageIntegrator(TwoPortBModel):
    """Ideal voltage integrator"""

    def __init__(self, Av=1):

        Av = LaplaceDomainExpression(Av)
        super(IdealVoltageIntegrator, self).__init__(
            BMatrix.voltage_integrator(Av))
        self.args = (Av, )


class CurrentFollower(TwoPortBModel):
    """Current follower"""

    def __init__(self):

        super(CurrentFollower, self).__init__(BMatrix.current_amplifier(1))
        self.args = ()


class IdealCurrentAmplifier(TwoPortBModel):
    """Ideal current amplifier"""

    def __init__(self, Ai=1):

        Ai = LaplaceDomainExpression(Ai)
        super(IdealCurrentAmplifier, self).__init__(
            BMatrix.current_amplifier(Ai))
        self.args = (Ai, )


class IdealCurrentDifferentiator(TwoPortBModel):
    """Ideal current differentiator"""

    def __init__(self, Ai=1):

        Ai = LaplaceDomainExpression(Ai)
        super(IdealCurrentDifferentiator, self).__init__(
            BMatrix.current_differentiator(Ai))
        self.args = (Ai, )


class IdealCurrentIntegrator(TwoPortBModel):
    """Ideal current integrator"""

    def __init__(self, Ai=1):

        Ai = LaplaceDomainExpression(Ai)
        super(IdealCurrentIntegrator, self).__init__(
            BMatrix.current_integrator(Ai))
        self.args = (Ai, )


class OpampInverter(TwoPortBModel):
    """Opamp inverter"""

    def __init__(self, R1, R2):

        R1 = ConstantDomainExpression(R1)
        R2 = ConstantDomainExpression(R2)
        # FIXME for initial voltages.
        super(OpampInverter, self).__init__(
            AMatrix(((-R1.Z / R2.Z, 0),
                     (-1 / R2.Z, 0).Bparams)))
        self.args = (R1, R2)


class OpampIntegrator(TwoPortBModel):
    """Inverting opamp integrator"""

    def __init__(self, R1, C1):

        R1 = ConstantDomainExpression(R1)
        C1 = ConstantDomainExpression(C1)
        # FIXME for initial voltages.
        super(OpampIntegrator, self).__init__(
            AMatrix(((-R1.Z / C1.Z, 0), (-1 / C1.Z, 0))).Bparams)
        self.args = (R1, C1)


class OpampDifferentiator(TwoPortBModel):
    """Inverting opamp differentiator"""

    def __init__(self, R1, C1):

        R1 = ConstantDomainExpression(R1)
        C1 = ConstantDomainExpression(C1)
        # FIXME for initial voltages.
        super(OpampDifferentiator, self).__init__(
            AMatrix(((-R1.Z * C1.Z, 0), (-R1.Z, 0))).Bparams)
        self.args = (R1, C1)


class TSection(TwoPortBModel):
    """T (Y) section
    ::

           +---------+       +---------+
         --+   OP1   +---+---+   OP3   +---
           +---------+   |   +---------+
                       +-+-+
                       |   |
                       |OP2|
                       |   |
                       +-+-+
                         |
         ----------------+-----------------

      The Z matrix for a resistive T section is
      [ R1 + R2, R2     ]
      [      R2, R2 + R3]
    """

    def __init__(self, OP1, OP2, OP3):

        _check_oneport_args((OP1, OP2, OP3))
        self.tp = Series(OP1).chain(Shunt(OP2)).chain(Series(OP3))

        super(TSection, self).__init__(self.tp)
        self.args = (OP1, OP2, OP3)

    def Pisection(self):

        ZV = WyeDelta(self.args[0].Z, self.args[1].Z, self.args[2].Z)
        VV = WyeDelta(self.args[0].V, self.args[1].V, self.args[2].V)
        OPV = [(ZV1 + VV1).cpt() for ZV1, VV1 in zip(ZV, VV)]

        return PiSection(*OPV)


class TwinTSection(TwoPortBModel):
    """Twin T section
    ::

              +---------+       +---------+
           +--+   OP1a  +---+---+   OP3a  +--+
           |  +---------+   |   +---------+  |
           |              +-+-+              |
           |              |   |              |
           |              |OP2a              |
           |              |   |              |
           |              +-+-+              |
           |                |                |
           |                v                |
           |  +---------+       +---------+  |
         --+--+   OP1b  +---+---+   OP3b  +--+--
              +---------+   |   +---------+
                          +-+-+
                          |   |
                          |OP2b
                          |   |
                          +-+-+
                            |
         -------------------+--------------------

    """

    def __init__(self, OP1a, OP2a, OP3a, OP1b, OP2b, OP3b):

        _check_oneport_args((OP1a, OP2a, OP3a, OP1b, OP2b, OP3b))
        self.tp = TSection(OP1a, OP2a, OP3a).parallel(
            TSection(OP1b, OP2b, OP3b))
        super(TwinTSection, self).__init__(self.tp)
        self.args = (OP1a, OP2a, OP3a, OP1b, OP2b, OP3b)


class BridgedTSection(TwoPortBModel):
    """Bridged T section
        ::

                       +---------+
           +-----------+   OP4   +-----------+
           |           +---------+           |
           |                                 |
           |  +---------+       +---------+  |
         --+--+   OP1b  +---+---+   OP3b  +--+--
              +---------+   |   +---------+
                          +-+-+
                          |   |
                          |OP2b
                          |   |
                          +-+-+
                            |
         -------------------+--------------------

         """

    def __init__(self, OP1, OP2, OP3, OP4):

        _check_oneport_args((OP1, OP2, OP3, OP4))
        self.tp = TSection(OP1, OP2, OP3).parallel(Series(OP4))
        super(TwinTSection, self).__init__(self.tp)
        self.args = (OP1, OP2, OP3, OP4)


class PiSection(TwoPortBModel):
    """Pi (delta) section
    ::

                  +---------+
        -----+----+   OP2    +---+-----
             |    +---------+   |
           +-+-+              +-+-+
           |   |              |   |
           |OP1|              |OP3|
           |   |              |   |
           +-+-+              +-+-+
             |                  |
        -----+------------------+-----

    """

    def __init__(self, OP1, OP2, OP3):

        self.tp = Shunt(OP1).chain(Series(OP2)).chain(Shunt(OP3))
        super(PiSection, self).__init__(self.tp)
        self.args = (OP1, OP2, OP3)

    def Tsection(self):

        ZV = DeltaWye(self.args[0].Z, self.args[1].Z, self.args[2].Z)
        VV = DeltaWye(self.args[0].V, self.args[1].V, self.args[2].V)
        OPV = [(ZV1 + VV1).cpt() for ZV1, VV1 in zip(ZV, VV)]
        return TSection(*OPV)


class LSection(TwoPortBModel):
    """L Section
    ::

           +---------+
         --+   OP1   +---+----
           +---------+   |
                       +-+-+
                       |   |
                       |OP2|
                       |   |
                       +-+-+
                         |
         ----------------+----
    """

    def __init__(self, OP1, OP2):

        _check_oneport_args((OP1, OP2))
        self.tp = Series(OP1).chain(Shunt(OP2))
        super(LSection, self).__init__(self.tp)
        self.args = (OP1, OP2)


class Ladder(TwoPortBModel):
    """(Unbalanced) ladder network with alternating Series and Shunt
    networks chained
    ::

           +---------+       +---------+
         --+   OP1   +---+---+ args[1] +---
           +---------+   |   +---------+
                       +-+-+
                       |   |
                       |   | args[0]
                       |   |
                       +-+-+
                         |
         ----------------+-----------------
    """

    def __init__(self, OP1, *args):

        _check_oneport_args((OP1, ) + args)

        self.tp = Series(OP1)

        for m, arg in enumerate(args):

            if m & 1:
                self.tp = self.tp.chain(Series(arg))
            else:
                self.tp = self.tp.chain(Shunt(arg))

        super(Ladder, self).__init__(self.tp)
        self.args = (OP1, ) + args

    def simplify(self):

        if len(self.args) == 1:
            return Series(self.args[0])
        elif len(self.args) == 2:
            return LSection(*self.args)
        elif len(self.args) == 3:
            return TSection(*self.args)
        return self

        # A Ladder of voltage sources and current sources
        # collapses to a single Lsection comprised of the total
        # voltage and total current.


class GeneralTransmissionLine(TwoPortBModel):
    """General transmission line

    Z0    characteristic impedance (ohms)
    gamma propagation constant (1/m)
    l     transmission line length (m)
    """

    reactive = True

    def __init__(self, Z0='Z0(s)', gamma='gamma(s)', l='l'):

        Z0 = LaplaceDomainExpression(Z0)
        gamma = LaplaceDomainExpression(gamma)
        l = ConstantDomainExpression(l)

        H = exp(gamma * l).expr

        B11 = sym.S.Half * (H + 1 / H)
        B12 = sym.S.Half * (1 / H - H) * Z0
        B21 = sym.S.Half * (1 / H - H) / Z0
        B22 = sym.S.Half * (H + 1 / H)

        B = BMatrix(((B11, B12), (B21, B22))).simplify()

        super(GeneralTransmissionLine, self).__init__(B)
        self.Z0 = Z0
        self.gamma = gamma
        self.l = l
        self.args = (Z0, gamma, l)


class GeneralTxLine(GeneralTransmissionLine):
    pass


class TL(GeneralTransmissionLine):
    """General transmission line"""

    # This class is used for netlists.
    pass


class LosslessTransmissionLine(GeneralTransmissionLine):
    """Losslees transmission line
        Z0 (real) characteristic impedance (ohms)
        c  propagation speed (m/s)
        l  transmission line length (m)
        """

    def __init__(self, Z0='Z0', c='c', l='l'):

        gamma = s / c

        super(LosslessTransmissionLine, self).__init__(Z0, gamma, l)
        self.args = (Z0, c, l)


class LosslessTxLine(LosslessTransmissionLine):
    pass


class TLlossless(LosslessTransmissionLine):
    """Lossless transmission line"""

    # This class is used for netlists.
    pass


class TransmissionLine(GeneralTransmissionLine):
    """Transmission line

    R series resistance/unit length
    L series inductance/unit length
    G shunt conductance/unit length
    C shunt capacitance/unit length
    l transmission line length
    """

    def __init__(self, R='R', L='L', G='G', C='C', l='l'):

        Z = R + s * L
        Y = G + s * C
        gamma = sqrt(Z * Y).expr
        Z0 = sqrt(Z / Y).expr

        super(TransmissionLine, self).__init__(Z0, gamma, l)
        self.args = (R, L, G, C, l)


class TxLine(TransmissionLine):
    pass
