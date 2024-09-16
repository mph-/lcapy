import sympy.core
from sympy import Mul

import lcapy.state
from sympy.physics.units.prefixes import PREFIXES, Prefix
from typing import Union

from lcapy import ConstantFrequencyResponseDomainExpression as cfrde
from lcapy import LaplaceDomainImpedance as ldi
from lcapy import ConstantTimeDomainVoltage as ctdv
from lcapy import TimeDomainVoltage as tdv
from lcapy import TimeDomainCurrent as tdc



class SIUnitPrefixer:
    def __init__(self):
        self.prefixes: dict[int, Prefix] =\
            {PREFIXES[prefixKey]._exponent: PREFIXES[prefixKey] for prefixKey in PREFIXES.keys()}

    @staticmethod
    def _findExponentFloatInt(value: Union[float, int]) -> int:
        if value == 0:
            return 0

        exponent = 0
        _value = abs(value)

        while _value >= 10:
            _value /= 10
            exponent += 1

        while _value < 1:
            _value *= 10
            exponent -= 1

        return exponent

    @staticmethod
    def _findExponentMul(value: Union[Mul, cfrde]) -> int:
        """
        this function assumes all symbols to be 1, to determine the prefix based on the numerical value in the
        expression if it receives a type it can not handle it returns 0
        """
        sub_dict = {sympy.sin: 1, sympy.cos: 1}

        for freeSymbol in value.free_symbols:
            sub_dict[freeSymbol] = 1

        try:
            if isinstance(value, Mul):
                if value.is_real:
                    _value = float(value.evalf(subs=sub_dict))
                else:
                    return 0
            elif isinstance(value, (cfrde, ctdv, tdv, tdc)):
                if value.expr.is_real:
                    _value = float(value.expr.evalf(subs=sub_dict))
                else:
                    return 0
            else:
                return 0
        except TypeError:
            return 0

        return SIUnitPrefixer._findExponentFloatInt(_value)

    def _findSIPrefix(self, exponent) -> Prefix:
        return self.prefixes[min(self.prefixes.keys(), key=lambda x: abs(x-exponent))]

    def getSIPrefix(self, value: Union[float, int, Mul, cfrde]) -> Prefix:
        if isinstance(value, (Mul, cfrde, lcapy.TimeDomainVoltage, lcapy.TimeDomainCurrent)):
            return self._findSIPrefix(self._findExponentMul(value))
        else:
            return self._findSIPrefix(self._findExponentFloatInt(value))

    def getSIPrefixedValue(self, value: Union[float, int, Mul, cfrde], minExponent=3):
        """
        add the nearest unit prefix to float, int, lcapy.ConstantFrequencyResponseDomainExpression or sympy.Mul
        prefixes are sympy.physics.units.prefixes.PREFIXES
        """

        if isinstance(value, (Mul, float, sympy.core.numbers.Rational)):
            expr = value
        elif isinstance(value, (cfrde, ctdv, tdv, tdc, ldi)):
            value = value
            expr = value.expr_with_units
        elif isinstance(value, (int, sympy.core.numbers.Integer)):
            value = expr = float(value)
        elif value is None:
            return None
        else:
            raise TypeError(f"value has to be type float, int or sympy.Mul or lcapy.ConstantFrequencyResponseDomainExpression not {type(value)}")

        prefix = self.getSIPrefix(value)
        exp = prefix._exponent

        # if this function would return value, the evalf() would remove the unit while converting to float
        if abs(exp) >= minExponent:
            return 1.0 * expr * 10**(-exp) * prefix
        else:
            return expr
