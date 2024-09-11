from typing import Union
from lcapy.unitPrefixer import SIUnitPrefixer


class JsonExportStepValues:
    def __init__(self, name1, name2, newName,
                 relation,
                 value1, value2, result,
                 latexEquation,
                 convVal1, convVal2, convResult):

        prefixer = SIUnitPrefixer()

        self.name1: str = name1
        self.name2: str = name2
        self.newName: str = newName
        self.relation: str = relation
        self.value1 = prefixer.getSIPrefixedValue(value1).evalf(n=3)
        self.value2 = prefixer.getSIPrefixedValue(value2).evalf(n=3)
        self.result = prefixer.getSIPrefixedValue(result).evalf(n=3)
        self.latexEquation: str = latexEquation

        if convVal1:
            self.convVal1 = prefixer.getSIPrefixedValue(convVal1).evalf(n=3)
        else:
            self.convVal1 = None
        if convVal2:
            self.convVal2 = prefixer.getSIPrefixedValue(convVal2).evalf(n=3)
        else:
            self.convVal2 = None
        if convResult:
            self.convResult = prefixer.getSIPrefixedValue(convResult).evalf(n=3)
        else:
            self.convResult = None

        self.hasConversion: bool = bool(bool(convVal1) or bool(convVal2) or bool(convResult))

    def toDict(self) -> dict:
        return {
            "name1": self.name1,
            "name2": self.name2,
            "newName": self.newName,
            "relation": self.relation,
            "value1": self.value1,
            "value2": self.value2,
            "result": self.result,
            "latexEquation": self.latexEquation,
            "hasConversion": self.hasConversion,
            "convVal1": self.convVal1,
            "convVal2": self.convVal2,
            "convResult": self.convResult
        }

