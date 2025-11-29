"""
This module supports simple linear four-port networks.  It is
experimental and needs a rethink.

Copyright 2025 Michael Hayes, UCECE

"""

from __future__ import division
from warnings import warn
import sympy as sym
from .exprclasses import ConstantDomainExpression
from .network import Network

__all__ = ('Transformer4', )


class Transformer4(Network):
    """Ideal transformer with four windings.  If the turns is negative,
    the dot is at the bottom of the winding.

    """

    def __init__(self, Ns2=1, Ns1=1, Np2=1, Np1=1):

        self.Ns2 = ConstantDomainExpression(Ns2)
        self.Ns1 = ConstantDomainExpression(Ns1)
        self.Np2 = ConstantDomainExpression(Np2)
        self.Np1 = ConstantDomainExpression(Np1)

        alpha3 = self.Np2 / self.Np1
        alpha2 = self.Ns2 / self.Np1
        alpha1 = self.Ns1 / self.Np1

        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.alpha3 = alpha3
        self.args = (Ns2, Ns1, Np2, Np1)


class TFsscss(Transformer4):
    pass
