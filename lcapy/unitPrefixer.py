import sympy.core
from sympy import Mul

import lcapy.state
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
    def _findExponentMul(value: lcapy.Expr, forImag=True) -> int:
        """
        this function assumes all symbols to be 1, to determine the prefix based on the numerical value in the
        expression if it receives a type it can not handle it returns 0. If it receives an expression that is an
        addition or multiplication one prefix may not be suitable so no prefix is chosen and the function returns 0
        :param forImag, if imag is True the prefix is also determined for expressions that only have an imaginary part
        an expression with imaginary part and real part 0 is returned
        """

        if value.expr.is_Add:
            return 0

        sub_dict = {sympy.sin: 1, sympy.cos: 1}

        for freeSymbol in value.free_symbols:
            sub_dict[freeSymbol] = 1

        try:
            if value.expr.is_real:
                _value = float(value.expr.evalf(subs=sub_dict))
            elif sympy.re(value.expr) == 0 and not sympy.im(value.expr) == 0 and forImag:
                _value = float(sympy.im(value.expr).evalf(subs=sub_dict))
            else:
                return 0
        except TypeError:
            return 0

        return SIUnitPrefixer._findExponentFloatInt(_value)

    def _findSIPrefix(self, exponent) -> Prefix:
        return self.prefixes[min(self.prefixes.keys(), key=lambda x: abs(x-exponent))]

    def getSIPrefix(self, value: Union[float, int, Mul, lcapy.Expr]) -> Prefix:
        if isinstance(value, (Mul, lcapy.Expr)):
            return self._findSIPrefix(self._findExponentMul(value))
        else:
            return self._findSIPrefix(self._findExponentFloatInt(value))

    def _getSIgetSIPrefixedValue(self, value: Union[float, int, Mul, lcapy.Expr], minExponent=3) \
            -> tuple[lcapy.Expr, Prefix, int]:
        if isinstance(value, lcapy.Expr):
            value = value
        else:
            value = lcapy.expr(value)

        expr = value.expr_with_units
        prefix = self.getSIPrefix(value)
        exp = prefix._exponent

        return lcapy.expr(expr), prefix, exp

    def getSIPrefixedExpr(self, value: Union[float, int, Mul, lcapy.Expr], minExponent=3) -> lcapy.Expr:
        """
        add the nearest unit prefix to float, int, lcapy.ConstantFrequencyResponseDomainExpression or sympy.Mul
        prefixes are sympy.physics.units.prefixes.PREFIXES
        """

        expr, prefix, exp = self._getSIgetSIPrefixedValue(value, minExponent)

        # if this function would return value, the evalf() would remove the unit while converting to float
        if abs(exp) >= minExponent:
            return lcapy.expr(1 * expr * 10**(-exp)) * prefix
        else:
            return lcapy.expr(expr)

    def getSIPrefixedMul(self, value: Union[float, int, Mul, lcapy.Expr], minExponent=3) -> Mul:
        """
        add the nearest unit prefix to float, int, lcapy.ConstantFrequencyResponseDomainExpression or sympy.Mul
        prefixes are sympy.physics.units.prefixes.PREFIXES
        """

        expr, prefix, exp = self._getSIgetSIPrefixedValue(value, minExponent)

        # if this function would return value, the evalf() would remove the unit while converting to float
        if abs(exp) >= minExponent:
            return 1.0 * expr.expr_with_units * 10**(-exp) * prefix
        else:
            return expr.expr_with_units
