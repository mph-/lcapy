import sympy
from sympy import latex
from lcapy import ConstantFrequencyResponseDomainExpression as cfrde

import lcapy
from lcapy import state
from sympy import Mul
from lcapy.impedanceConverter import ImpedanceToComponent
from lcapy.impedanceConverter import ValueToComponent
from lcapy.impedanceConverter import getSourcesFromCircuit, getOmegaFromCircuit
from lcapy.netlistLine import NetlistLine
from sympy.physics.units import Hz
from sympy import parse_expr
from lcapy import omega0
from lcapy.jsonExportStepValues import JsonExportStepValues
from lcapy.unitWorkAround import UnitWorkAround as uwa
from lcapy.unitPrefixer import SIUnitPrefixer


class JsonExport:
    def __init__(self, precision=4):
        self.name1 = None
        self.name2 = None
        self.newName = None
        self.thisStep = None
        self.lastStep = None
        self.step = None
        self.cpt1 = None  # component1
        self.cpt2 = None  # component2
        self.cptRes = None  # componentResult
        self.valCpt1 = None  # valueComponent1
        self.valCpt2 = None  # valueComponent2
        self.valCptRes = None  # valueComponentResult
        self.convValCpt1 = None  # convertedValueComponent1 -> converted from Impedance to R,L or C if possible
        self.convValCpt2 = None  # convertedValueComponent2 -> converted from Impedance to R,L or C if possible
        self.convValCptRes = None  # convertedValueComponentResult -> converted from Impedance to R,L or C if possible
        self.cvc1Type = None  # convertedValueComponent1Type
        self.cvc2Type = None  # convertedValueComponent2Type
        self.cvcrType = None  # convertedValueComponentResultType
        self.omega_0 = None

        self.precision = precision

        self.prefixer = SIUnitPrefixer()
        self.prefixedLatexStr = lambda x: latex(
            self.prefixer.getSIPrefixedValue(x).evalf(n=self.precision),
            imaginary_unit="j"
        )
        self.latexStr = lambda x: latex(x.evalf(n=self.precision), imaginary_unit="j")

    def updateObjectValues(self, step, solution: 'lcapy.Solution'):
        self.name1 = solution[step].cpt1
        self.name2 = solution[step].cpt2
        self.newName = solution[step].newCptName
        self.thisStep = solution[step]
        self.lastStep = solution[step].lastStep
        self.omega_0 = getOmegaFromCircuit(solution[step].circuit, getSourcesFromCircuit(solution[step].circuit))

        if not self._isInitialStep():
            self.cpt1 = self.lastStep.circuit[self.name1]
            self.cpt2 = self.lastStep.circuit[self.name2]
            self.cptRes = self.thisStep.circuit[self.newName]

            self.valCpt1 = str(solution.getElementSpecificValue(self.cpt1))
            self.valCpt2 = str(solution.getElementSpecificValue(self.cpt2))
            self.valCptRes = str(solution.getElementSpecificValue(self.cptRes))

            self.convValCpt1, self.cvc1Type = ValueToComponent(self.valCpt1, omega_0=self.omega_0)
            self.convValCpt2, self.cvc2Type = ValueToComponent(self.valCpt2, omega_0=self.omega_0)
            self.convValCptRes, self.cvcrType = ValueToComponent(self.valCptRes, omega_0=self.omega_0)

    def getDictForStep(self, step, solution: 'lcapy.Solution'):
        self.updateObjectValues(step, solution)

        if self._isInitialStep():
            return self._handleInitialStep()

        elif self._checkEssentialInformation():
            raise ValueError(f"missing information in {step}: "
                             f"{self.name1}, {self.name2}, {self.newName}, {self.thisStep}, {self.lastStep}")

        else:
            state.show_units = True

            if self._allValuesConvertableToComponent():
                if self._isSameType():
                    values = self._handleSameTypeAndConvertibleToComponent()
                else:
                    values = self._handleDifferentTypeAndConvertibleToComponent()

            elif self._resultConvertibleToComponent():
                values = self._handleResultConvertibleToComponent()

            else:
                values = self._handleNoConversionPossible()

            for key in ["value1", "value2", "result", "convVal1", "convVal2", "convResult"]:
                if values[key]:
                    values[key] = self.latexStr(values[key])

            return values

    def _isInitialStep(self) -> bool:
        return not (self.name1 and self.name2 and self.newName and self.lastStep) and self.thisStep

    def _handleInitialStep(self) -> dict:

        # this is the initial step which is used as an overview of the circuit
        as_dict = {}
        state.show_units = True

        for cptName in self.thisStep.circuit.elements.keys():
            cpt = self.thisStep.circuit.elements[cptName]
            if cpt.type == "V":
                as_dict[cptName] = latex(
                    uwa.addUnit(
                        NetlistLine(str(cpt)).value,
                        cpt.type
                    )
                )

                if cpt.has_ac:
                    if cpt.args[2] is not None:
                        as_dict["omega_0"] = latex(self.prefixer.getSIPrefixedValue(
                            parse_expr(str(cpt.args[2]), local_dict={"pi": sympy.pi}) * Hz)
                        )
                        try:
                            self.omega_0 = float(cpt.args[2])
                        except ValueError:
                            self.omega_0 = str(cpt.args[2])
                    else:
                        as_dict["omega_0"] = latex(omega0)
                        self.omega_0 = "omega_0"
                elif cpt.has_dc:
                    as_dict["omega_0"] = latex(Mul(0) * Hz)
                else:
                    raise AssertionError("Voltage Source is not ac or dc")

            elif not cpt.type == "W":
                cCpt = NetlistLine(ImpedanceToComponent(str(cpt), omega_0=self.omega_0))
                as_dict[cCpt.type + cCpt.typeSuffix] = latex(
                    self.prefixer.getSIPrefixedValue(
                        uwa.addUnit(
                            cCpt.value,
                            cCpt.type
                        )
                    )
                )
        return as_dict

    def _checkEssentialInformation(self) -> bool:
        """
        this function makes sure that all information that is needed to compute a solution step is available,
        exception is the initial step that does acquire the information it needs.
        :returns: true if information is available, false otherwise
        """
        return not (self.name1 or self.name2 or self.newName or self.lastStep or self.thisStep or self.step)

    def _allValuesConvertableToComponent(self) -> bool:
        return (not self.valCpt1 == self.convValCpt1
                and not self.valCpt2 == self.convValCpt2
                and not self.valCptRes == self.convValCptRes)

    def _isSameType(self) -> bool:
        return self.cvc1Type == self.cvc2Type

    def _handleSameTypeAndConvertibleToComponent(self) -> dict:
        eqVal1 = uwa.addUnit(self.convValCpt1, self.cvc1Type)
        eqVal2 = uwa.addUnit(self.convValCpt2, self.cvc2Type)
        eqRes = uwa.addUnit(self.convValCptRes, self.cvcrType)
        compType = self.cvc1Type
        assert compType in ["R", "L", "C"]

        equation = self.makeLatexEquation(eqVal1, eqVal2, eqRes, self.thisStep.relation, compType)

        return JsonExportStepValues(self.name1, self.name2, self.newName, self.thisStep.relation,
                                    eqVal1, eqVal2, eqRes, equation,
                                    convVal1=None, convVal2=None, convResult=None).toDict()

    def _handleDifferentTypeAndConvertibleToComponent(self) -> dict:
        eqVal1 = uwa.addUnit(self.valCpt1, "Z")
        eqVal2 = uwa.addUnit(self.valCpt2, "Z")
        eqRes = uwa.addUnit(self.valCptRes, "Z")
        compType = self.cptRes.type
        assert compType == "Z"
        convValCptRes = uwa.addUnit(self.convValCptRes, self.cvcrType)

        equation = self.makeLatexEquation(eqVal1, eqVal2, eqRes, self.thisStep.relation, compType)

        return JsonExportStepValues(self.name1, self.name2, self.newName, self.thisStep.relation,
                                    eqVal1, eqVal2, eqRes, equation,
                                    convVal1=None, convVal2=None, convResult=convValCptRes).toDict()

    def _resultConvertibleToComponent(self) -> bool:
        return (self.valCpt1 == self.convValCpt1
                and self.valCpt2 == self.convValCpt2
                and not self.valCptRes == self.convValCptRes)

    def _handleResultConvertibleToComponent(self) -> dict:
        eqVal1 = uwa.addUnit(self.valCpt1, "Z")
        eqVal2 = uwa.addUnit(self.valCpt2, "Z")
        eqRes = uwa.addUnit(self.valCptRes, "Z")
        compType = self.cpt1.type
        assert compType == "Z"
        convValCptRes = uwa.addUnit(self.convValCptRes, self.cvcrType)

        equation = self.makeLatexEquation(eqVal1, eqVal2, eqRes, self.thisStep.relation, compType)

        return JsonExportStepValues(self.name1, self.name2, self.newName, self.thisStep.relation,
                                    eqVal1, eqVal2, eqRes, equation,
                                    convVal1=None, convVal2=None, convResult=convValCptRes).toDict()

    def _handleNoConversionPossible(self) -> dict:
        eqVal1 = uwa.addUnit(self.valCpt1, "Z")
        eqVal2 = uwa.addUnit(self.valCpt2, "Z")
        eqRes = uwa.addUnit(self.valCptRes, "Z")
        compType = self.cpt1.type
        assert compType == "Z"

        equation = self.makeLatexEquation(eqVal1, eqVal2, eqRes, self.thisStep.relation, compType)

        return JsonExportStepValues(self.name1, self.name2, self.newName, self.thisStep.relation,
                                    eqVal1, eqVal2, eqRes, equation,
                                    convVal1=None, convVal2=None, convResult=None).toDict()

    def makeLatexEquation(self, exp1: cfrde, exp2: cfrde, expRslt: cfrde, cptRelation: str, compType: str) \
            -> str:

        if compType not in ["R", "L", "C", "Z"]:
            raise ValueError(f"{compType} is unknown, component type has to be R, L, C or Z")

        # inverse sum means 1/x1 + 1/x2 = 1/_xresult e.g parallel resistor
        parallelRel = {"R": "inverseSum", "C": "sum", "L": "inverseSum", "Z": "inverseSum"}
        rowRel = {"R": "sum", "C": "inverseSum", "L": "sum", "Z": "sum"}

        state.show_units = True
        if cptRelation == "parallel":
            useFunc = parallelRel[compType]
        elif cptRelation == "series":
            useFunc = rowRel[compType]
        else:
            raise AttributeError(
                f"Unknown relation between elements {cptRelation}. Known relations are: parallel, series"
            )

        if not compType == "Z":
            expStr1 = self.prefixedLatexStr(exp1)
            expStr2 = self.prefixedLatexStr(exp2)
            expStrRslt = self.prefixedLatexStr(expRslt)
        else:
            expStr1 = self.latexStr(exp1)
            expStr2 = self.latexStr(exp2)
            expStrRslt = self.latexStr(expRslt)

        if useFunc == "inverseSum":
            equation = "\\frac{1}{" + expStr1 + "} + \\frac{1}{" + expStr2 + "} = " + expStrRslt
        elif useFunc == "sum":
            equation = expStr1 + " + " + expStr2 + " = " + expStrRslt
        else:
            raise NotImplementedError(f"Unknown function {useFunc}")

        return equation
