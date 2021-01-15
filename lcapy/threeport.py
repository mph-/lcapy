"""
This module supports simple linear three-port networks.  It is
experimental and needs a rethink.

Copyright 2014--2021 Michael Hayes, UCECE

"""

from __future__ import division
from warnings import warn
import sympy as sym
from .sexpr import LaplaceDomainVoltage, LaplaceDomainTransferFunction
from .sexpr import LaplaceDomainCurrent
from .smatrix import LaplaceDomainVoltageMatrix, LaplaceDomainCurrentMatrix
from .smatrix import LaplaceDomainImpedanceMatrix, LaplaceDomainAdmittanceMatrix
from .cexpr import cexpr
from .oneport import OnePort
from .twoport import YMatrix, ZMatrix, TwoPortZModel, Series, TwoPort

__all__ = ('Opamp', )


class ThreePortMatrix(sym.Matrix):

    def __new__(cls, *args):

        if len(args) == 9:
            args = ((args[0], args[1], args[2]),
                    (args[3], args[4], args[5]),
                    (args[6], args[7], args[8]))
            return super(ThreePortMatrix, cls).__new__(cls, args)

        return super(ThreePortMatrix, cls).__new__(cls, *args)

    @property
    def Z(self):
        return ZMatrix3(self.Y.inv())

    @property
    def Y(self):
        return YMatrix3(self.Z.inv())


class ZMatrix3(ThreePortMatrix):

    """
    +-  -+     +-             -+   +-  -+
    | V1 |     | Z11  Z12  Z13 |   | I1 |
    | V2 |  =  | Z21  Z22  Z23 |   | I2 |
    | V3 |     | Z31  Z32  Z33 |   | I3 |
    +-  -+     +-             -+   +-  -+

    Z = inv(Y)
    """

    @property
    def Z(self):
        # Perhaps we should make a copy?
        return self


class YMatrix3(ThreePortMatrix):

    """
    +-  -+     +-             -+   +-  -+
    | I1 |  =  | Y11  Y12  Y13 |   | V1 |
    | I2 |     | Y21  Y22  Y23 |   | V2 |
    | I3 |     | Y31  Y32  Y33 |   | V3 |
    +-  -+     +-             -+   +-  -+

    Y = inv(Z)
    """

    @property
    def Y(self):
        # Perhaps we should make a copy?
        return self


class ThreePort(object):

    """

    +-  -+     +-             -+   +-  -+     +-   -+
    | V1 |     | Z11  Z12  Z13 |   | I1 |     | V1z |
    | V2 |  =  | Z21  Z22  Z23 |   | I2 |  +  | V2z |
    | V3 |     | Z31  Z32  Z33 |   | I3 |     | V3z |
    +-  -+     +-             -+   +-  -+     +-   -+

    The A, B, G, and H models are invalid for multiports.

    Unfortunately, the Z model can blow up for simple networks.

    """

    def __init__(self, Z, Vz=None):

        if Vz is None:
            Vz = LaplaceDomainVoltageMatrix((0, 0, 0))

        if not isinstance(Z, ZMatrix3):
            raise ValueError('Z not ZMatrix3')

        if not isinstance(Vz, LaplaceDomainVoltageMatrix):
            raise ValueError('Vz not LaplaceDomainVoltageMatrix')

        self._M = Z
        self._Vz = Vz

    @property
    def Voc(self):
        """Return voltage vector with all ports open-circuited
        (i.e., In = 0)"""
        return self._Vz

    @property
    def Isc(self):
        """Return current vector with all ports short-circuited
        (i.e., Vn = 0)"""
        Y = self.Y
        Voc = self.Voc

        Isc = LaplaceDomainCurrentMatrix([Voc[m] * Y[m, m] for m in range(len(Voc))])
        return Isc

    @property
    def Y(self):
        """Return admittance matrix"""
        return YMatrix3(self._M.Y)

    @property
    def Z(self):
        """Return impedance matrix"""
        return self._M

    @property
    def Yoc(self):
        """Return admittance vector with ports open circuit"""
        Z = self.Z
        return LaplaceDomainAdmittanceMatrix([1 / Z[m, m] for m in range(Z.shape[0])])

    @property
    def Ysc(self):
        """Return admittance vector with ports short circuit"""
        Y = self.Y
        return LaplaceDomainAdmittanceMatrix([Y[m, m] for m in range(Y.shape[0])])

    @property
    def Zoc(self):
        """Return impedance vector with ports open circuit"""
        Z = self.Z
        return LaplaceDomainImpedanceMatrix([Z[m, m] for m in range(Z.shape[0])])

    @property
    def Zsc(self):
        """Return impedance vector with ports short circuit"""
        Y = self.Y
        return LaplaceDomainImpedanceMatrix([1 / Y[m, m] for m in range(Y.shape[0])])

    def _portcheck(self, port):

        if port not in (1, 2, 3):
            raise ValueError('Invalid port ' + port)

    def Vgain(self, inport=1, outport=2):
        """Return voltage gain for specified ports with internal
        sources zero"""

        self._portcheck(inport)
        self._portcheck(outport)

        p1 = inport - 1
        p2 = outport - 1

        return LaplaceDomainTransferFunction(self.Z[p2, p1] / self.Z[p1, p1])

    def Igain(self, inport=1, outport=2):
        """Return voltage gain for specified ports with internal
        sources zero"""

        self._portcheck(inport)
        self._portcheck(outport)

        p1 = inport - 1
        p2 = outport - 1

        return LaplaceDomainTransferFunction(self.Y[p2, p1] / self.Y[p1, p1])

    def Vresponse(self, V, inport=1, outport=2):
        """Return voltage response for specified applied voltage and
        specified ports"""

        self._portcheck(inport)
        self._portcheck(outport)

        p1 = inport - 1
        p2 = outport - 1

        H = self.Z[p2, p1] / self.Z[p1, p1]
        return LaplaceDomainVoltage(self.Voc[p2] + (V - self.Voc[p1]) * H)

    def Iresponse(self, I, inport=1, outport=2):
        """Return current response for specified current voltage and
        specified ports"""

        self._portcheck(inport)
        self._portcheck(outport)

        p1 = inport - 1
        p2 = outport - 1

        Y = self.Y
        Isc = self.Isc

        return LaplaceDomainCurrent(Isc[p2] + (I - Isc[p1]) * Y[p2, p1] / Y[p1, p1])

    def attach_parallel(self, OP, port=2):
        """Attach one-port in parallel to specified port"""

        if not issubclass(OP.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        self._portcheck(port)

        p = port - 1

        Y = self.Y
        Y[p, p] += OP.Ys
        Isc = self.Isc
        Isc[p] += OP.Isc
        Z = Y.Z
        Voc = LaplaceDomainVoltageMatrix([LaplaceDomainVoltage(Isc[m] * Z[m, m]) for m in range(len(Isc))])
        return ThreePort(Z, Voc)

    def bridge(self, OP, inport=1, outport=2):
        """Bridge the specified ports with a one-port element"""

        self._portcheck(inport)
        self._portcheck(outport)

        # Create two-port series element.
        s = Series(OP)

        # The impedance matrix for a series element is infinite.

        Y3 = YMatrix3(((0, 0, 0), (0, 0, 0), (0, 0, 0)))
        Y2 = s.Y
        p1 = inport - 1
        p2 = outport - 1

        Y3[p1, p1] = Y2[0, 0]
        Y3[p2, p2] = Y2[1, 1]
        Y3[p1, p2] = Y2[0, 1]
        Y3[p2, p1] = Y2[1, 0]

        Y = self.Y + Y3
        Isc = self.Isc
        Isc[p1] -= OP.Isc
        Isc[p2] += OP.Isc
        Z = Y.Z
        Voc = LaplaceDomainVoltageMatrix([LaplaceDomainVoltage(Isc[m] * Z[m, m]) for m in range(len(Isc))])
        return ThreePort(Y.Z, Voc)

    def parallel(self, MP, port=None):
        """Return the model with, MP, in parallel"""

        if issubclass(MP.__class__, OnePort):
            return self.attach_parallel(MP, port)

        if issubclass(MP.__class__, TwoPort):
            # We could special case a series or shunt network here.
            raise NotImplementedError('TODO')

        if not issubclass(MP.__class__, ThreePort):
            raise TypeError('Argument not ', ThreePort)

        Y = self.Y + MP.Y
        Isc = self.Isc + MP.Isc
        Z = Y.Z
        Voc = LaplaceDomainVoltageMatrix([LaplaceDomainVoltage(Isc[m] * Z[m, m]) for m in range(len(Isc))])
        return ThreePort(Z, Voc)

    def series(self, MP, port=None):
        """Return the model with, MP, in series"""

        if issubclass(MP.__class__, OnePort):
            raise NotImplementedError('TODO')

        if issubclass(MP.__class__, TwoPort):
            raise NotImplementedError('TODO')

        if not issubclass(MP.__class__, ThreePort):
            raise TypeError('Argument not ', ThreePort)

        warn('Will this ever work?')

        Z = self.Z + MP.Z
        Voc = self.Voc + MP.Voc

        return ThreePort(Z, Voc)

    def terminate(self, OP, port=2):
        """Connect one-port in parallel to specified port and return a
        two-port object"""

        return self.attach_parallel(OP, port).open_circuit(port)

    def short_circuit(self, port=2):
        """Apply a short-circuit to specified port and return a
        two-port object"""

        # Remove the unwanted port from the Ymatrix.
        Y = self.Y.copy()
        Y.row_del(port - 1)
        Y.col_del(port - 1)
        Y = YMatrix(Y)

        # CHECKME, perhaps use Isc?
        Voc = self.Voc.copy()
        Voc.row_del(port - 1)

        return TwoPortZModel(Y.Z, LaplaceDomainVoltage(Voc[0]), LaplaceDomainVoltage(Voc[1]))

    def open_circuit(self, port=2):
        """Apply a open-circuit to specified port and return a
        two-port object"""

        # Remove the unwanted port from the Zmatrix.
        Z = self.Z.copy()
        Z.row_del(port - 1)
        Z.col_del(port - 1)
        Z = ZMatrix(Z)

        Voc = self.Voc.copy()
        Voc.row_del(port - 1)

        return TwoPortZModel(Z, LaplaceDomainVoltage(Voc[0]), LaplaceDomainVoltage(Voc[1]))


class Opamp(ThreePort):

    """
    Create an ideal(ish) opamp
    ::

            |\
            |  \
        1 --+ +  \
            |      \
            |        \
            |          +--- 3
            |        /
            |      /
        2 --+ -  /
            |  /
            |/

        Each port voltage, Vn, is referenced to a common ground.
        Port 1: non-inverting input
        Port 2: inverting input
        Port 3: output

    """

    def __init__(self, Rd=1e9, Ro=1e-6, A=100000, Rp=1e9, Rm=1e9):

        # If Ro=0, then Z matrix singular.

        Rd, Ro, A, Rp, Rm = [cexpr(arg) for arg in (Rd, Ro, A, Rp, Rm)]

        Ra = Rp * (Rd + Rm) / (Rp + Rd + Rm)
        Rb = Rm * (Rd + Rp) / (Rp + Rd + Rm)

        Z = ZMatrix3(((Rp + Rd, Rd, 0),
                      (Rd, Rm + Rd, 0),
                      (A * Ra, -A * Rb, Ro)))
        super(Opamp, self).__init__(Z)
        self.args = (Rd, Ro, A, Rp, Rm)
