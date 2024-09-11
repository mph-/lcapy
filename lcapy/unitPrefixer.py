import sympy.core
from sympy import Mul

import lcapy.state
from lcapy import ConstantFrequencyResponseDomainExpression as cfrde
from lcapy import ConstantTimeDomainVoltage as ctdv
from sympy.physics.units.prefixes import PREFIXES, Prefix
from typing import Union


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
        sub_dict = {}

        for freeSymbol in value.free_symbols:
            sub_dict[freeSymbol] = 1
        try:
            if isinstance(value, Mul):
                if value.is_real:
                    _value = float(value.evalf(subs=sub_dict))
                else:
                    return 0
            elif isinstance(value, cfrde):
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
        if isinstance(value, (Mul, cfrde)):
            return self._findSIPrefix(self._findExponentMul(value))
        else:
            return self._findSIPrefix(self._findExponentFloatInt(value))

    def getSIPrefixedValue(self, value: Union[float, int, Mul, cfrde], minExponent=3):
        """
        add the nearest unit prefix to float, int, lcapy.ConstantFrequencyResponseDomainExpression or sympy.Mul
        prefixes are sympy.physics.units.prefixes.PREFIXES
        """

        if isinstance(value, (Mul, float)):
            value = expr = value
        elif isinstance(value, (cfrde, ctdv)):
            value = value
            expr = value.expr_with_units
        elif isinstance(value, int):
            value = expr = float(value)
        elif value is None:
            return None
        else:
            raise TypeError("value has to be type float, int or sympy.Mul or "
                            f"lcapy.ConstantFrequencyResponseDomainExpression not {type(value)}")

        prefix = self.getSIPrefix(value)
        exp = prefix._exponent

        if abs(exp) >= minExponent:
            return 1.0 * expr * 10**(-exp) * prefix
        else:
            # if this returns value the evalf() function to convert to float removes the unit
            return expr
