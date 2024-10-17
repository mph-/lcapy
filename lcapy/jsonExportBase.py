from sympy.printing import latex
from sympy import Float
from typing import Union
from sympy import Mul
from lcapy import Expr


class JsonExportBase:
    def latexWithPrefix(self, value: Union[Mul, Expr], prec=None, addPrefix: bool = True) -> str:
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

    def latexWithoutPrefix(self, value: Expr, prec=None) -> str:
        return self.latexWithPrefix(value, prec, addPrefix=False)

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
