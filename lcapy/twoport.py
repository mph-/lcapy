"""
This module supports simple linear two-port networks.

Copyright 2014--2021 Michael Hayes, UCECE
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
from .functions import Eq, MatMul


# This needs to be generalised for superpositions.
# For now, it uses only the s-domain.

# TODO:
# 1. Defer the choice of the two-port model.  For example, a T section
# would store the three sub-networks rather than generating a B matrix.
# The appropriate model would be generated when desired.  This would
# avoid the inversion of singular matrices. The downside is that each
# object would require methods to generate each type of two-port model.
#
# Some multiport networks, such as a shunt R, have a singular Z matrix.
# Thus switching to the Y matrix and back to the Z matrix produces a
# bogus result.  The same thing occurs for a series R; this has a
# singular Y matrix.
#
# 2. Fix handling of buffered two ports (amplifier / delay).


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
           'Series', 'Shunt', 'IdealTransformer', 'IdealGyrator',
           'VoltageFollower', 'VoltageAmplifier',
           'IdealVoltageAmplifier', 'IdealDelay',
           'IdealVoltageDifferentiator', 'IdealVoltageIntegrator',
           'CurrentFollower', 'IdealCurrentAmplifier',
           'IdealCurrentDifferentiator', 'IdealCurrentIntegrator',
           'OpampInverter', 'OpampIntegrator', 'OpampDifferentiator',
           'TSection', 'TwinTSection', 'BridgedTSection', 'PiSection',
           'LSection', 'Ladder', 'GeneralTxLine', 'LosslessTxLine',
           'TxLine', 'AMatrix', 'BMatrix', 'GMatrix', 'HMatrix',
           'SMatrix', 'TMatrix', 'YMatrix', 'ZMatrix')

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
        

class TwoPortMatrix(Matrix, TwoPortMixin):

    # The following default properties are fallbacks when other conversions have
    # not been defined.

    @property
    def Aparams(self):
        return AMatrix(self.Bparams.inv())

    @property
    def Bparams(self):
        return BMatrix(self.Aparams.inv())

    @property
    def Gparams(self):
        return GMatrix(self.Hparams.inv())

    @property
    def Hparams(self):
        return HMatrix(self.Gparams.inv())

    @property
    def Sparams(self):
        return self.Aparams.Sparams

    @property
    def Tparams(self):
        return self.Sparams.Tparams

    @property
    def Yparams(self):
        return YMatrix(self.Zparams.inv())

    @property
    def Zparams(self):
        return ZMatrix(self.Yparams.inv())

    def pdb(self):
        import pdb; pdb.set_trace()
        return self
    

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

        # Inverse
        det = self.det().expr
        if det == 0:
            warn('Producing dodgy B matrix')
        return BMatrix(((self._A22 / det, -self._A12 / det),
                        (-self._A21 / det, self._A11 / det)))

    @property
    def Hparams(self):

        if self._A22 == 0:
            warn('Producing dodgy H matrix')
        det = self.det().expr            
        return HMatrix(((self._A12 / self._A22, det / self._A22),
                        (-1 / self._A22, self._A21 / self._A22)))


    @property
    def Sparams(self):
        Z0 = LaplaceDomainImpedance('Z_0').as_expr()
        d = self._A12 + Z0 * (self._A11 + self._A22) + Z0**2 * self._A21
        return SMatrix((((self._A12 + Z0 * (self._A11 - self._A22) - Z0**2 * self._A21) / d,
                         (2 * Z0 * (self._A11 * self._A22 - self._A12 * self._A21)) / d),
                         (2 * Z0 / d, (self._A12 - Z0 * (self._A11 - self._A22) - Z0**2 * self._A21) / d)))


    @property
    def Yparams(self):

        # This produces a bogus Y matrix when A12 is zero (say for a
        # shunt element).   Note, it doesn't use A21.
        if self._A12 == 0:
            warn('Producing dodgy Y matrix')
        det = self.det().expr            
        return YMatrix(((self._A22 / self._A12, -det / self._A12),
                        (-1 / self._A12, self._A11 / self._A12)))

    @property
    def Zparams(self):

        # This produces a bogus Z matrix when A21 is zero (say for a
        # series element).   Note, it doesn't use A12.
        if self._A21 == 0:
            warn('Producing dodgy Z matrix')
        det = self.det().expr            
        return ZMatrix(((self._A11 / self._A21, det / self._A21),
                        (1 / self._A21, self._A22 / self._A21)))

    @property
    def Z1oc(self):
        """open-circuit input impedance"""
        # Z11
        return LaplaceDomainImpedance(self._A11 / self._A21)

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

    def chain(self, OP):

        return self * OP

    def cascade(self, OP):

        return self.chain(OP)


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
        # Inverse
        det = self.det().expr
        return AMatrix(((self._B22 / det, -self._B12 / det),
                        (-self._B21 / det, self._B11 / det)))

    @property
    def Bparams(self):
        # Perhaps we should make a copy?
        return self

    @property
    def Gparams(self):

        det = self.det().expr        
        return GMatrix(((-self._B21 / self._B22, -1 / self._B22),
                        (det / self._B22, -self._B12 / self._B22)))

    @property
    def Hparams(self):

        det = self.det().expr        
        return HMatrix(((-self._B12 / self._B11, 1 / self._B11),
                        (-det / self._B11, -self._B21 / self._B11)))

    @property
    def Yparams(self):

        det = self.det().expr
        return YMatrix(((-self._B11 / self._B12, 1 / self._B12),
                        (det / self._B12, -self._B22 / self._B12)))

    @property
    def Zparams(self):

        det = self.det().expr        
        return ZMatrix(((-self._B22 / self._B21, -1 / self._B21),
                        (-det / self._B21, -self._B11 / self._B21)))

    @property
    def Z1oc(self):
        """open-circuit input impedance"""
        # Z11
        return LaplaceDomainImpedance(-self._B22 / self._B21)

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

        return cls(((1 / Ar, -1 / (Ar * Yin),
                    ( -1 / (Ar * Zout)), -1 / (Ar * Yin * Zout * (Af * Ar - 1)))))

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

        alpha = ConstantDomainExpression(alpha)

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
                        (self._G11 / self._G21, det / self._G21)))

    @property
    def Bparams(self):
        # return self.Hparams.Bparams
        det = self.det().expr        
        return BMatrix(((-det / self._G12),
                        (self._G22 / self._G12, self._G11 / self._G12, -1 / self._G12)))

    @property
    def Gparams(self):
        # Perhaps we should make a copy?
        return self

    @property
    def Hparams(self):
        return HMatrix(self.inv())

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
                        (-self._H22 / self._H21, -1 / self._H21)))

    @property
    def Bparams(self):
        det = self.det().expr
        return BMatrix(((1 / self._H12, -self._H11 / self._H12),
                        (-self._H22 / self._H12, det / self._H12)))

    @property
    def Hparams(self):
        # Perhaps we should make a copy?
        return self

    @property
    def Yparams(self):
        det = self.det().expr        
        return YMatrix(((1 / self._H11, -self._H12 / self._H11),
                        (self._H21 / self._H11, det / self._H11)))

    @property
    def Zparams(self):
        det = self.det().expr        
        return ZMatrix(((det / self._H22, self._H12 / self._H22),
                        (-self._H21 / self._H22, 1 / self._H22)))


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
                        (1 / self._T22, -self._T21 / self._T22)))
    
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
                        (-det / self._Y21, -self._Y11 / self._Y21)))

    @property
    def Bparams(self):
        det = self.det().expr        
        return BMatrix(((-self._Y11 / self._Y12, 1 / self._Y12),
                        (det / self._Y12, -self._Y22 / self._Y12)))

    @property
    def Hparams(self):
        det = self.det().expr
        return HMatrix(((1 / self._Y11, -self._Y12 / self._Y11),
                        (self._Y21 / self._Y11, det / self._Y11)))

    @property
    def Yparams(self):
        # Perhaps we should make a copy?
        return self

    @property
    def Zparams(self):
        # Inverse
        det = self.det().expr
        return ZMatrix(((self._Y22 / det, -self._Y12 / det),
                        (-self._Y21 / det, self._Y11 / det)))


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
    def Zoc(self):
        return LaplaceDomainImpedanceMatrix((self._Z11, self._Z22))

    @property
    def Aparams(self):
        det = self.det().expr
        return AMatrix(((self._Z11 / self._Z21, det / self._Z21),
                        (1 / self._Z21, self._Z22 / self._Z21)))

    @property
    def Bparams(self):
        det = self.det().expr
        return BMatrix(((self._Z22 / self._Z12, -det / self._Z12),
                        (-1 / self._Z12, self._Z11 / self._Z12)))

    @property
    def Hparams(self):
        det = self.det().expr
        return HMatrix(((det / self._Z22, self._Z12 / self._Z22),
                        (-self._Z21 / self._Z22, 1 / self._Z22)))

    @property
    def Yparams(self):
        # Inverse
        det = self.det().expr
        return YMatrix(((self._Z22 / det, -self._Z12 / det),
                        (-self._Z21 / det, self._Z11 / det)))

    @property
    def Zparams(self):
        # Perhaps we should make a copy?
        return self

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


class TwoPort(Network, TwoPortMixin):

    """
    General class for two-port networks.  Two-port networks are
    constrained to have the same current at each port (but flowing in
    opposite directions).  This is called the port condition.
    """

    def _add_elements(self):
        raise ValueError('Cannot generate netlist for two-port objects')

    def netlist(self, form='horizontal', evalf=None):
        raise ValueError('Cannot generate netlist for two-port objects')

    def _check_twoport_args(self):

        # This is an interim measure until Par2, Ser2, etc. generalised.
        if len(self.args) != 2:
            raise ValueError('Only two args supported for %s' %
                             self.__class__.__name__)
        for arg1 in self.args:
            if not isinstance(arg1, TwoPort):
                raise ValueError('%s not a TwoPort' % arg1)

    @property
    def Aparams(self):
        """Return chain parameters"""
        return self._M.Aparams

    @property
    def Bparams(self):
        """Return inverse chain parameters"""
        return self._M.Bparams

    @property
    def Gparams(self):
        """Return inverse hybrid parameters"""
        return self._M.Gparams

    @property
    def Hparams(self):
        """Return hybrid parameters"""
        return self._M.Hparams

    @property
    def Sparams(self):
        """Return scattering parameters"""
        return self._M.Sparams

    @property
    def Tparams(self):
        """Return scattering transfer parameters"""
        return self._M.Tparams    

    @property
    def Yparams(self):
        """Return admittance parameters"""
        return self._M.Yparams

    @property
    def Zparams(self):
        """Return impedance parameters"""
        return self._M.Zparams

    @property
    def I1a(self):
        return LaplaceDomainCurrent(-self.V2b / self._B12)

    @property
    def V1a(self):
        # CHECKME
        return LaplaceDomainVoltage(-self.I2b / self._B21)

    @property
    def I1g(self):
        return LaplaceDomainCurrent(-self.I2b / self._B22)

    @property
    def V2g(self):
        return LaplaceDomainVoltage(self.V2b - self._B12 / self._B22 * self.I2b)

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
        """Return input impedance with the output port open circuit"""
        return LaplaceDomainImpedance(1 / self.Z1oc)

    @property
    def Y2oc(self):
        """Return output impedance with the input port open circuit"""
        return LaplaceDomainAdmittance(1 / self.Z2oc)

    @property
    def Ysc(self):
        """Return admittance vector with ports short circuit"""
        return self.Yparams.Ysc

    @property
    def Y1sc(self):
        """Return input admittance with output port short circuit"""
        return LaplaceDomainAdmittance(self.Ysc[0])

    @property
    def Y2sc(self):
        """Return output admittance with output port short circuit"""
        return LaplaceDomainAdmittance(self.Ysc[1])

    @property
    def Zoc(self):
        """Return impedance vector with ports open circuit"""
        return self.Zparams.Zoc

    @property
    def Z1oc(self):
        """Return input impedance with the output port open circuit"""
        return LaplaceDomainImpedance(self.Zoc[0])

    @property
    def Z2oc(self):
        """Return output impedance with the input port open circuit"""
        return LaplaceDomainImpedance(self.Zoc[1])

    @property
    def Zsc(self):
        """Return impedance vector with ports short circuit"""
        return LaplaceDomainImpedanceMatrix((LaplaceDomainImpedance(1 / self.Y1sc),
                                             LaplaceDomainImpedance(1 / self.Y2sc)))

    @property
    def Z1sc(self):
        """Return input impedance with the output port short circuit"""
        return LaplaceDomainImpedance(1 / self.Y1sc)

    @property
    def Z2sc(self):
        """Return output impedance with the input port short circuit"""
        return LaplaceDomainImpedance(1 / self.Y2sc)

    def Vgain(self, inport=1, outport=2):
        """Return voltage gain for specified ports with internal
        sources zero"""

        # Av  = G21 = 1 / A11 = -det(B) / B22 = Z21 / Z11 =  Y21 / Y22
        # Av' = H12 = 1 / B11 =  |A| / A22 = Z12 / Z22 = -Y12 / Y11

        if inport == outport:
            return LaplaceDomainTransferFunction(1)
        if inport == 1 and outport == 2:
            return LaplaceDomainTransferFunction(1 / self._A11)
        if inport == 2 and outport == 1:
            return LaplaceDomainTransferFunction(1 / self._B11)
        raise ValueError('bad port values')

    def Igain(self, inport=1, outport=2):
        """Return current gain for specified ports with internal
         sources zero"""

        # Ai  = H21 = -1 / A22 = -det(B) / B11 = -Z21 / Z22 = Y21 / Y11
        # Ai' = G12 =  1 / B22 =  |A| / A11 = -Z12 / Z11 = Y12 / Y22

        if inport == outport:
            return LaplaceDomainTransferFunction(1)
        if inport == 1 and outport == 2:
            return LaplaceDomainTransferFunction(-1 / self._A22)
        if inport == 2 and outport == 1:
            return LaplaceDomainTransferFunction(-1 / self._B22)
        raise ValueError('bad port values')

    @property
    def Vgain12(self):
        """Return V2 / V1 for I2 = 0 (forward voltage gain) with
        internal sources zero

        Av = G21 = 1 / A11 = -det(B) / B22 = Z21 / Z11 =  Y21 / Y22
        """

        return self.Vgain(1, 2)

    @property
    def Vtransfer(self):
        """Return V2 / V1 for I2 = 0 (forward voltage gain) with
        internal sources zero  (see Vgain12)"""

        return self.Vgain12

    @property
    def Igain12(self):
        """Return I2 / I1 for V2 = 0 (forward current gain) with
        internal sources zero

        Ai = H21 = -1 / A22 = -det(B) / B11 = -Z21 / Z22 = Y21 / Y11
        """

        return self.Igain(1, 2)

    @property
    def Itransfer(self):
        """Return I2 / I1 for V2 = 0 (forward current gain) with
        internal sources zero  (sett Igain12)"""

        return self.Igain12

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
    def Bmodel(self):

        return TwoPortBModel(self.Bparams, self.V2b, self.I2b)

    @property
    def Hmodel(self):

        return TwoPortHModel(self.Hparams, self.V1h, self.I2h)

    @property
    def Ymodel(self):

        if self.is_shunt:
            warn('Converting a shunt two-port to a Y model is dodgy...')
        return TwoPortYModel(self.Yparams, self.I1y, self.I2y)

    @property
    def Zmodel(self):

        if self.is_series:
            warn('Converting a series two-port to a Z model is dodgy...')
        return TwoPortZModel(self.Zparams, self.V1z, self.V2z)

    def chain(self, TP):
        """Return the model with, TP, appended (cascade or
        chain connection)"""

        if not issubclass(TP.__class__, TwoPort):
            raise TypeError('Argument not', TwoPort)

        return Chain(self, TP)

    def append(self, TP):
        """Return the model with, TP, appended"""

        return self.chain(TP)

    def prepend(self, TP):
        """Return the model with, TP, prepended"""

        return TP.chain(self)

    def cascade(self, TP):
        """Return the model with, TP, appended"""

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

    def bridge(self, TP):
        """Bridge the ports with a one-port element"""

        if not issubclass(TP.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        # FIXME
        return self.parallel(Series(TP))

    def load(self, TP):
        """Apply a one-port load and return a Thevenin (one-port) object"""

        if not issubclass(TP.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        foo = self.chain(Shunt(TP))
        return Z(foo.Z1oc) + V(foo.V1oc)

    def source(self, TP):
        """Apply a one-port source and return a Thevenin (one-port) object"""

        if not issubclass(TP.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        foo = Shunt(TP).chain(self)
        return Z(foo.Z2oc) +  V(foo.V2oc)

    def short_circuit(self, port=2):
        """Apply a short-circuit to specified port and return a
        one-port object"""

        p = port - 1
        Yval = self.Yparams[1 - p, 1 - p]
        Ival = self.Isc[1 - p]

        return Y(Yval) | I(Ival)

    def open_circuit(self, port=2):
        """Apply a open-circuit to specified port and return a
        one-port object"""

        p = port - 1
        Zval = self.Zparams[1 - p, 1 - p]
        Vval = self.Voc[1 - p]

        return Z(Zval) + V(Vval)

    def simplify(self):

        if self.Bparams == sym.eye(2):
            # Have a pair of wires... perhaps could simplify
            # to an LSection comprised of a V and I but
            # may have a weird voltage expression.
            pass
        return self


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

    # A disadvantage of the Y and Z matrices is that they become
    # singular for some simple networks.  For example, the Z matrix is
    # singular for a shunt element and the Y matrix is singular for a
    # series element.  The A and B matrices do not seem to have this
    # problem, however, they cannot be extended to three or more ports.
    #

    def __init__(self, B, V2b=None, I2b=None):

        if V2b is None:
            V2b = LaplaceDomainVoltage(0)
        if I2b is None:
            I2b = LaplaceDomainCurrent(0)            

        if issubclass(B.__class__, TwoPortBModel):
            B, V2b, I2b = B._M, B._V2b, B._I2b

        if not isinstance(B, BMatrix):
            raise ValueError('B not BMatrix')

        if not isinstance(V2b, LaplaceDomainVoltage):
            raise ValueError('V2b not LaplaceDomainVoltage')

        if not isinstance(I2b, LaplaceDomainCurrent):
            raise ValueError('I2b not LaplaceDomainCurrent')

        super(TwoPortBModel, self).__init__()
        self._M = B
        self._V2b = V2b
        self._I2b = I2b

    @property
    def Bparams(self):
        """Return chain matrix"""
        return self._M

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


class TwoPortGModel(TwoPort):

    """
    """

    def __init__(self, G, I1g=None, V2g=None):

        if I1g is None:
            I1g = LaplaceDomainCurrent(0)
        if V2g is None:
            V2g = LaplaceDomainVoltage(0)            
        
        if issubclass(G.__class__, TwoPortGModel):
            G, I1g, V2g = G._M, G._I1g, G._V2g

        if not isinstance(G, GMatrix):
            raise ValueError('G not GMatrix')

        if not isinstance(I1g, LaplaceDomainCurrent):
            raise ValueError('I1g not LaplaceDomainCurrent')

        if not isinstance(V2g, LaplaceDomainVoltage):
            raise ValueError('V2g not LaplaceDomainVoltage')

        super(TwoPortGModel, self).__init__()
        self._M = G
        self._V1g = I1g
        self._I2g = V2g

    @property
    def Gparams(self):
        """Return hybrid matrix"""
        return self._M

    @property
    def V2b(self):
        """Return V2b"""

        # FIXME
        return LaplaceDomainVoltage(self.I1g / self.Gparams._G12)

    @property
    def I2b(self):
        """Return I2b"""

        # FIXME
        return LaplaceDomainCurrent(self.G._G22 / self.G._G12 * self.I1g) - self.V2g

    @property
    def I1g(self):
        return self._V1g

    @property
    def V2g(self):
        return self._I2g


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

    def __init__(self, H, V1h=None, I2h=None):

        if V1h is None:
            V1h = LaplaceDomainVoltage(0)
        if I2h is None:
            I2h = LaplaceDomainCurrent(0)            
        
        if issubclass(H.__class__, TwoPortHModel):
            H, V1h, I2h = H._M, H._V1h, H._I2h

        if not isinstance(H, HMatrix):
            raise ValueError('H not HMatrix')

        if not isinstance(V1h, LaplaceDomainVoltage):
            raise ValueError('V1h not LaplaceDomainVoltage')

        if not isinstance(I2h, LaplaceDomainCurrent):
            raise ValueError('I2h not LaplaceDomainCurrent')

        super(TwoPortHModel, self).__init__()
        self._M = H
        self._V1h = V1h
        self._I2h = I2h

    @property
    def Hparams(self):
        """Return hybrid matrix"""
        return self._M

    @property
    def V2b(self):
        """Return V2b"""

        return LaplaceDomainVoltage(self.V1h / self.Hparams._H12)

    @property
    def I2b(self):
        """Return I2b"""

        return LaplaceDomainCurrent(self.H._H22 / self.H._H12 * self.V1h) - self.I2h

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

    def __init__(self, Y, I1y=None, I2y=None):

        if I1y is None:
            I1y = LaplaceDomainCurrent(0)
        if I2y is None:
            I2y = LaplaceDomainCurrent(0)            
        
        if issubclass(Y.__class__, TwoPortYModel):
            Y, I1y, I2y = Y._M, Y._I1y, Y._I2y

        if not isinstance(Y, YMatrix):
            raise ValueError('Y not YMatrix')

        if not isinstance(I1y, LaplaceDomainCurrent):
            raise ValueError('I1y not LaplaceDomainCurrent')
        if not isinstance(I2y, LaplaceDomainCurrent):
            raise ValueError('I2y not LaplaceDomainCurrent')

        super(TwoPortYModel, self).__init__()
        self._M = Y
        self._I1y = I1y
        self._I2y = I2y

    @property
    def Yparams(self):
        """Return admittance matrix"""
        return self._M

    @property
    def I2b(self):
        return LaplaceDomainCurrent(-self.I1y * self._Y11 * self._Y22 / self._Y12) - self.I2y

    @property
    def V2b(self):
        return LaplaceDomainVoltage(self.I1y * self._Y11 / self._Y22)

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

    def __init__(self, Z, V1z=None, V2z=None):

        if V1z is None:
            V1z = LaplaceDomainVoltage(0)
        if V2z is None:
            V2z = LaplaceDomainVoltage(0)            
        
        if issubclass(Z.__class__, TwoPortZModel):
            Z, V1z, V2z = Z._M, Z._V1z, Z._V2z

        if not isinstance(Z, ZMatrix):
            raise ValueError('Z not ZMatrix')

        if not isinstance(V1z, LaplaceDomainVoltage):
            raise ValueError('V1z not LaplaceDomainVoltage')
        if not isinstance(V2z, LaplaceDomainVoltage):
            raise ValueError('V2z not LaplaceDomainVoltage')

        super(TwoPortZModel, self).__init__()
        self._M = Z
        self._V1z = V1z
        self._V2z = V2z

    @property
    def Zparams(self):
        """Return impedance matrix"""
        return self._M

    @property
    def I2b(self):
        return LaplaceDomainCurrent(self.V1z / self._Z12)

    @property
    def V2b(self):
        return self.V2z - LaplaceDomainVoltage(self.V1z * self._Z22 / self._Z12)

    @property
    def I1y(self):

        Zdet = self.Zparams.det().expr
        return LaplaceDomainCurrent(-self.V1z * self._Z22 / Zdet - self.V2z * self._Z12 / Zdet)

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


class Chain(TwoPortBModel):

    """Connect two-port networks in a chain (aka cascade)"""

    def __init__(self, *args):

        self.args = args
        self._check_twoport_args()

        arg1 = args[-1]
        B = arg1.Bparams

        foo = Vector(arg1.V2b, arg1.I2b)

        for arg in reversed(args[0:-1]):

            foo += B * Vector(arg.V2b, arg.I2b)
            B = B * arg.Bparams

        super(Chain, self).__init__(B, LaplaceDomainVoltage(foo[0, 0]), LaplaceDomainCurrent(foo[1, 0]))

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

        self.args = args
        self._check_twoport_args()

        # This will fail with a Shunt as an argument since it does
        # not have a valid Y model.
        # We can special case this.
        if isinstance(args[0], Shunt) or isinstance(args[1], Shunt):
            print('Warning: need to handle a Shunt in parallel')

        arg = args[0]
        I1y = arg.I1y
        I2y = arg.I2y
        Y = arg.Yparams

        for arg in args[1:]:
            I1y += arg.I1y
            I2y += arg.I2y
            Y += arg.Yparams

        super(Par2, self).__init__(Y, I1y, I2y)

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

        self.args = args
        self._check_twoport_args()

        # Need to be more rigorous.
        if isinstance(self.args[1], (Series, LSection, TSection)):
            print('Warning: This can violate the port condition')

        arg = args[0]
        V1z = arg.V1z
        V2z = arg.V2z
        Z = arg.Zparams

        for arg in args[1:]:
            V1z += arg.V1z
            V2z += arg.V2z
            Z += arg.Zparams

        super(Ser2, self).__init__(Z, V1z, V2z)

    def simplify(self):

        if isinstance(self.args[0], Shunt) and isinstance(self.args[1], Shunt):
            return Shunt(
                (self.args[0].args[0] + self.args[1].args[0]).simplify())

        return self


class Hybrid2(TwoPortHModel):

    """Connect two-port networks in hybrid configuration (inputs in
    series, outputs in parallel)"""

    def __init__(self, *args):

        self.args = args
        self._check_twoport_args()

        arg = args[0]
        V1h = arg.V1h
        I2h = arg.I2h
        H = arg.Hparams

        for arg in args[1:]:
            V1h += arg.V1h
            I2h += arg.I2h
            H += arg.Hparams

        super(Hybrid2, self).__init__(H, V1h, I2h)


class InverseHybrid2(TwoPortGModel):

    """Connect two-port networks in inverse hybrid configuration (outputs in
    series, inputs in parallel)"""

    def __init__(self, *args):

        self.args = args
        self._check_twoport_args()

        arg = args[0]
        I1g = arg.I1g
        V2g = arg.V2g
        G = arg.Gparams

        for arg in args[1:]:
            I1g += arg.I1g
            V2g += arg.V2g
            G += arg.Gparams

        super(Hybrid2, self).__init__(G, I1g, V2g)


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

        self.OP = OP
        self.args = (OP, )
        _check_oneport_args(self.args)
        super(Series, self).__init__(BMatrix.Zseries(OP.Z.laplace()), LaplaceDomainVoltage(OP.Voc.laplace()), LaplaceDomainCurrent(0))


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

        self.OP = OP
        self.args = (OP, )
        _check_oneport_args(self.args)
        super(Shunt, self).__init__(BMatrix.Yshunt(OP.Y.laplace()), LaplaceDomainVoltage(0),
                                    LaplaceDomainCurrent(OP.Isc.laplace()))


class IdealTransformer(TwoPortBModel):

    """Ideal transformer voltage gain alpha, current gain 1 / alpha"""

    def __init__(self, alpha=1):

        self.alpha = ConstantDomainExpression(alpha)
        self.args = (alpha, )
        super(IdealTransformer, self).__init__(BMatrix.transformer(alpha))


class TF(IdealTransformer):
    pass


class IdealGyrator(TwoPortBModel):

    """Ideal gyrator with gyration resistance R.

    A gyrator converts a voltage to current and a current to voltage.
    Cascaded gyrators act like a transformer"""

    def __init__(self, R=1):

        self.R = ConstantDomainExpression(R)
        self.args = (R, )
        super(IdealGyrator, self).__init__(BMatrix.gyrator(R))


class VoltageFollower(TwoPortBModel):

    """Voltage follower"""

    def __init__(self):

        self.args = ()
        super(VoltageFollower, self).__init__(BMatrix.voltage_amplifier(1))


class VoltageAmplifier(TwoPortBModel):

    """Voltage amplifier"""

    def __init__(self, Av=1, Af=0, Yin=0, Zout=0):

        Av = LaplaceDomainExpression(Av)
        Af = LaplaceDomainExpression(Af)
        Yin = LaplaceDomainExpression(Yin)
        Zout = LaplaceDomainExpression(Zout)

        self.args = (Av, Af, Yin, Zout)
        super(VoltageAmplifier, self).__init__(
            BMatrix.voltage_amplifier(Av, Af, Yin, Zout))


class IdealVoltageAmplifier(VoltageAmplifier):

    """Ideal voltage amplifier"""

    def __init__(self, Av=1):

        Av = LaplaceDomainExpression(Av)
        super(IdealVoltageAmplifier, self).__init__(
            BMatrix.voltage_differentiator(Av))
        self.args = (Av, )


class IdealDelay(TwoPortBModel):

    """Ideal buffered delay"""

    def __init__(self, delay=0):

        delay = ConstantDomainExpression(delay)
        super(IdealDelay, self).__init__(
            BMatrix.voltage_amplifier(sym.exp(-s * delay)))
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

        self.args = (OP1, OP2, OP3)
        _check_oneport_args(self.args)
        super(TSection, self).__init__(
            Series(OP1).chain(Shunt(OP2)).chain(Series(OP3)))

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

        self.args = (OP1a, OP2a, OP3a, OP1b, OP2b, OP3b)
        _check_oneport_args(self.args)

        super(TwinTSection, self).__init__(
            TSection(OP1a, OP2a, OP3a).parallel(TSection(OP1b, OP2b, OP3b)))


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

        self.args = (OP1, OP2, OP3, OP4)
        _check_oneport_args(self.args)

        super(TwinTSection, self).__init__(
            TSection(OP1, OP2, OP3).parallel(Series(OP4)))


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

        super(PiSection, self).__init__(
            Shunt(OP1).chain(Series(OP2)).chain(Shunt(OP3)))
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

        self.args = (OP1, OP2)
        _check_oneport_args(self.args)

        super(LSection, self).__init__(Series(OP1).chain(Shunt(OP2)))


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

        self.args = (OP1, ) + args
        _check_oneport_args(self.args)

        TP = Series(OP1)

        for m, arg in enumerate(args):

            if m & 1:
                TP = TP.chain(Series(arg))
            else:
                TP = TP.chain(Shunt(arg))

        super(Ladder, self).__init__(TP)

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


class GeneralTxLine(TwoPortBModel):

    """General transmission line

    Z0 is the (real) characteristic impedance (ohms)
    gamma is the propagation constant (1/m)
    l is the transmission line length (m)
    """

    def __init__(self, Z0, gamma, l):

        Z0 = LaplaceDomainExpression(Z0)
        gamma = LaplaceDomainExpression(gamma)
        l = ConstantDomainExpression(l)

        H = sym.exp(gamma * l)

        B11 = 0.5 * (H + 1 / H)
        B12 = 0.5 * (1 / H - H) * Z0
        B21 = 0.5 * (1 / H - H) / Z0
        B22 = 0.5 * (H + 1 / H)

        super(GeneralTxLine, self).__init__(BMatrix(((B11, B12), (B21, B22))))
        self.args = (Z0, gamma, l)


class LosslessTxLine(GeneralTxLine):

    """Losslees transmission line
        Z0 is the (real) characteristic impedance (ohms)
        c is the propagation speed (m/s)
        l is the transmission line length (m)
        """

    def __init__(self, Z0, c=1.5e8, l=1):

        gamma = s / c

        super(LosslessTxLine, self).__init__(Z0, gamma, l)


class TxLine(GeneralTxLine):

    """Transmission line

    R series resistance/metre
    L series inductance/metre
    G shunt conductance/metre
    C shunt capacitance/metre
    l is the transmission line length
    """

    def __init__(self, R, L, G, C, l=1):

        Z = R + s * L
        Y = G + s * C
        gamma = sym.sqrt(Z * Y)
        Z0 = sym.sqrt(Z / Y)

        super(TxLine, self).__init__(Z0, gamma, l)
