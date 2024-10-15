from sympy.printing import latex
from sympy import Float
from lcapy import Expr


class JsonExportBase:
    def latexWithPrefix(self, value: Expr, prec=None, addPrefix: bool = True) -> str:
        if prec is None:
            prec = self.precision

        if addPrefix:
            toPrint = 1.0 * self.prefixer.getSIPrefixedMul(value)
        else:
            toPrint = 1.0 * value.expr_with_units

        for val in list(toPrint.atoms(Float)):
            toPrint = toPrint.evalf(subs={val: str(round(val, prec))})
        latexString = latex(toPrint, imaginary_unit="j")
        return latexString

    def latexWithoutPrefix(self, value: Expr, prec=None) -> str:
        return self.latexWithPrefix(value, prec, addPrefix=False)

    def _getValueFieldKeys(self) -> list[str]:
        keys = list(self.__dict__.keys())
        valueFiledKeys = []
        for key in keys:
            if key.find("val") >= 0 or key.find("Val") >= 0:
                valueFiledKeys.append(key)

        return valueFiledKeys
