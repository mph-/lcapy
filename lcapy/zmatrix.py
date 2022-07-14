"""This module implements the ZDomainMatrix class for a matrix of
z-domain expressions.

Copyright 2021 Michael Hayes, UCECE

"""

from .matrix import Matrix


class ZDomainMatrix(Matrix):
    from .zexpr import ZDomainExpression
    _typewrap = ZDomainExpression

    def IZT(self, **assumptions):

        def func(expr):
            return expr.IZT(**assumptions)

        return DiscreteTimeDomainMatrix(self.applyfunc(func))

    def canonical(self):

        return self.applyfunc(self._typewrap.canonical)

    def general(self):

        return self.applyfunc(self._typewrap.general)

    def mixedfrac(self):

        return self.applyfunc(self._typewrap.mixedfrac)

    def partfrac(self):

        return self.applyfunc(self._typewrap.partfrac)

    def timeconst(self):

        return self.applyfunc(self._typewrap.timeconst)

    def ZPK(self):

        return self.applyfunc(self._typewrap.ZPK)


class ZDomainVoltageMatrix(ZDomainMatrix):
    from .zexpr import ZDomainVoltage
    _typewrap = ZDomainVoltage


class ZDomainCurrentMatrix(ZDomainMatrix):
    from .zexpr import ZDomainCurrent
    _typewrap = ZDomainCurrent


class ZDomainAdmittanceMatrix(ZDomainMatrix):
    from .zexpr import ZDomainAdmittance
    _typewrap = ZDomainAdmittance


class ZDomainImpedanceMatrix(ZDomainMatrix):
    from .zexpr import ZDomainImpedance
    _typewrap = ZDomainImpedance


from .nmatrix import DiscreteTimeDomainMatrix  # nopep8
