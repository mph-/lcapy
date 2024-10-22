from sympy.printing import print_latex
from sympy.printing import latex
from sympy import Float
from typing import Union
from sympy import Mul
from lcapy import Expr
from lcapy.unitPrefixer import SIUnitPrefixer


class JsonExportBase:
    def __init__(self, precision: int):
        self.precision = precision
        self.prefixer = SIUnitPrefixer()

    def _latexRealNumber(self, value: Union[Mul, Expr], prec=None, addPrefix: bool = True) -> str:
        if prec is None:
            prec = self.precision

        if addPrefix:
            toPrint = 1.0 * self.prefixer.getSIPrefixedMul(value)
        else:
            if isinstance(value, Expr):
                toPrint = 1.0 * value.expr_with_units
            else:
                toPrint = 1.0 * value

        for val in list(toPrint.atoms(Float)):
            toPrint = toPrint.evalf(subs={val: str(round(val, prec))})
        latexString = latex(toPrint, imaginary_unit="j")
        return latexString

    @staticmethod
    def _latexComplexNumber(value: Union[Mul, Expr]):

        test = latex(value.evalf(n=3, chop=True))
        return test

    def latexWithPrefix(self, value: Union[Mul, Expr], prec=None, addPrefix: bool = True) -> str:
        if value.is_Add:
            return self._latexComplexNumber(value)
        else:
            return self._latexRealNumber(value, prec, addPrefix)

    def latexWithoutPrefix(self, value: Expr, prec=None) -> str:
        if value.is_Add:
            return self._latexComplexNumber(value)
        else:
            return self._latexRealNumber(value, prec, addPrefix=False)

    def _getValueFieldKeys(self, *args: str) -> list[str]:
        """
        finds fields that include the strings of args in their name to automatically convert them to a latex string
        on export. All fields are converted to lowercase so this functino is not case-sensitive.
        :return: list of keys<str> that have the name of the fields that match the criteria
        """

        keys = list(self.__dict__.keys())
        valueFiledKeys = []
        for key in keys:
            lcKey = key.lower()
            if any(arg.lower() in lcKey for arg in args):
                valueFiledKeys.append(key)

        return valueFiledKeys
