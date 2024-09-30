import lcapy
from sympy import latex
from lcapy import state
from enum import Enum
from lcapy.netlistLine import NetlistLine
from lcapy import t
from lcapy.impedanceConverter import ValueToComponent
from lcapy.unitPrefixer import SIUnitPrefixer
from lcapy.unitWorkAround import UnitWorkAround as uwa
from lcapy.componentRelation import ComponentRelation
from lcapy.impedanceConverter import getSourcesFromCircuit, getOmegaFromCircuit


class JsonVCValueExport:
    """
    Jason Volt Current Value Export
    Creates a json-File with information about the Voltage and Current values for one step
    The format is handy to display the simplification on the Web-UI.
    Takes a step <string> that is part of a Solution <lcapy.Solution> Object. The available steps can be accessed by
    <lcapy.Solution>.available_steps
    represents all information in the circuit in combination with lcapy.JsonCompValueExport. Has to be split up because
    in the user based mode not all information can be known when those files are generated
    """
    def __init__(self, precision=4):
        # this class automatically prefixes every field that includes val or Val in the name and transforms it to
        # a latex string before exporting the dictionary
        self.circuit: 'lcapy.Circuit' = None
        self.simpCircuit: 'lcapy.Circuit' = None
        self.prefixer = SIUnitPrefixer()
        self.precision = precision
        self.omega_0 = None

        self.oldNames = {'CompName': None, 'Uname': None, 'Iname': None}
        self.names1 = {'CompName': None, 'Uname': None, 'Iname': None}
        self.names2 = {'CompName': None, 'Uname': None, 'Iname': None}
        self.oldValues = {'Z': None, 'U': None, 'I': None}  # values to Component with oldName
        self.values1 = {'Z': None, 'U': None, 'I': None}  # values to Component with name1
        self.values2 = {'Z': None, 'U': None, 'I': None}  # values to Component with name2
        self.convOldValue = None  # converted Value of Z from oldValues
        self.convValue1 = None  # converted Value of Z from values1
        self.convValue2 = None  # converted Value of Z from values2
        self.relation: ComponentRelation = ComponentRelation.none
        self.equation = {'CompName1': None, 'CompName2': None}

        self.valueFieldKeys = self._getValueFieldKeys()

    def getDictForStep(self, step: str, solution: 'lcapy.Solution') -> dict:
        self._updateObjectValues(step, solution)

        as_dict = {
            'oldNames': list(self.oldNames.values()),
            'names1': list(self.names1.values()),
            'names2': list(self.names2.values()),
            'oldValues': list(self.oldValues.values()),
            'values1': list(self.values1.values()),
            'values2': list(self.values2.values()),
            'convOldValue': [self.convOldValue],
            'convValue1': [self.convValue1],
            'convValue2': [self.convValue2],
            'relation': [self.relation.to_string()],
            'equation': list(self.equation.values())
        }

        try:
            for key in self.valueFieldKeys:
                assert isinstance(as_dict[key], list)
                for idx, val in enumerate(as_dict[key]):
                    if as_dict[key][idx] is not None:
                        as_dict[key][idx] = self.latexWithPrefix(val)
        except KeyError:
            raise AssertionError(f"A filed which name includes val or Val is not in the export dict. Key: {key}")

        return as_dict

    def _updateObjectValues(self, step: str, solution: 'lcapy.Solution'):
        self.solStep: 'lcapy.solutionStep' = solution[step]
        self.simpCircuit: 'lcapy.Circuit' = solution[step].circuit  # circuit with less elements (n elements)
        self.omega_0 = getOmegaFromCircuit(self.simpCircuit, getSourcesFromCircuit(self.simpCircuit))

        # self.olfName, self.name1, self.name2 ist to access the component in the lcapy circuits
        # self.oldNames['CompName'], self.names1['CompName'], self.names2['CompName'] ist the adjusted value for the
        # export and is transformed from e.g. Z1 to R1, L1 or C1
        self.oldName = self.oldNames['CompName'] = solution[step].newCptName
        self.name1 = self.names1['CompName'] = solution[step].cpt1
        self.name2 = self.names2['CompName'] = solution[step].cpt2

        if not self._isInitialStep():
            self.circuit: 'lcapy.Circuit' = solution[step].lastStep.circuit  # circuit with more elements (n+1 elements)

            self._updateSuffix()
            self._updateUnames()
            self._updateInames()
            self._updateZandConvValues()
            self._updateU()
            self._updateI()
            self._updateCompRel()
            self._updateEquations()
            self._addPrefixes()

    def latexWithPrefix(self, value):
        prefixedValue = self.prefixer.getSIPrefixedValue(value)
        evalValue = prefixedValue.evalf(n=self.precision)
        latexString = latex(evalValue, imaginary_unit="j")
        return latexString

    def _updateSuffix(self):
        self.suffixOldName = NetlistLine(str(self.simpCircuit[self.oldName])).typeSuffix
        self.suffixName1 = NetlistLine(str(self.circuit[self.name1])).typeSuffix
        self.suffixName2 = NetlistLine(str(self.circuit[self.name2])).typeSuffix

    def _updateUnames(self):
        self.oldNames['Uname'] = 'U' + self.suffixOldName
        self.names1['Uname'] = 'U' + self.suffixName1
        self.names2['Uname'] = 'U' + self.suffixName2

    def _updateInames(self):
        self.oldNames['Iname'] = 'I' + self.suffixOldName
        self.names1['Iname'] = 'I' + self.suffixName1
        self.names2['Iname'] = 'I' + self.suffixName2

    def _updateZandConvValues(self):
        self.oldValues['Z'], self.convOldValue, compType = self._checkForConversion(self.simpCircuit[self.oldName].Z)
        self.oldNames['CompName'] = compType + self.suffixOldName
        self.values1['Z'], self.convValue1, compType = self._checkForConversion(self.circuit[self.name1].Z)
        self.names1['CompName'] = compType + self.suffixName1
        self.values2['Z'], self.convValue2, compType = self._checkForConversion(self.circuit[self.name2].Z)
        self.names2['CompName'] = compType + self.suffixName2

    def _checkForConversion(self, value) -> tuple:
        convValue, convCompType = ValueToComponent(value, self.omega_0)
        if 'Z' == convCompType:
            return value, None, 'Z'
        else:
            return value, uwa.addUnit(convValue, convCompType), convCompType

    def _updateU(self):
        self.oldValues['U'] = self.simpCircuit[self.oldName].V(t)
        self.values1['U'] = self.circuit[self.name1].V(t)
        self.values2['U'] = self.circuit[self.name2].V(t)

    def _updateI(self):
        self.oldValues['I'] = self.simpCircuit[self.oldName].I(t)
        self.values1['I'] = self.circuit[self.name1].I(t)
        self.values2['I'] = self.circuit[self.name2].I(t)

    def _updateCompRel(self):
        if self.solStep.relation == ComponentRelation.parallel:
            self.relation = ComponentRelation.parallel
        elif self.solStep.relation == ComponentRelation.series:
            self.relation = ComponentRelation.series
        else:
            self.relation = ComponentRelation.none

    def _makeLatexEquationU(self, valueZ, valueI, resultU) -> str:
        return f"{self.latexWithPrefix(valueZ)} • {self.latexWithPrefix(valueI)} = {self.latexWithPrefix(resultU)}"

    def _makeLatexEquationI(self, valueZ, valueU, resultI) -> str:
        return f"{self.latexWithPrefix(valueU)} / {self.latexWithPrefix(valueZ)} = {self.latexWithPrefix(resultI)}"

    def _updateEquations(self):
        if self.relation == ComponentRelation.parallel:
            # parallel components voltage is the same
            self.equation['CompName1'] = self._makeLatexEquationI(self.values1['Z'], self.values1['U'], self.values1['I'])
            self.equation['CompName2'] = self._makeLatexEquationI(self.values2['Z'], self.values2['U'], self.values2['I'])
        elif self.relation == ComponentRelation.series:
            self.equation['CompName1'] = self._makeLatexEquationU(self.values1['Z'], self.values1['I'], self.values1['U'])
            self.equation['CompName2'] = self._makeLatexEquationU(self.values2['Z'], self.values2['I'], self.values2['U'])
        else:
            if not self.relation == ComponentRelation.none:
                raise AttributeError(f"Unknown component relation type: {self.relation} supported types listed in "
                                     f"lcapy.componentRelation")
            else:
                self.equation['CompName1'] = ""
                self.equation['CompName2'] = ""

    def _isInitialStep(self) -> bool:
        assert isinstance(self.solStep, lcapy.SolutionStep)
        return not (self.solStep.cpt1 and self.solStep.cpt2
                    and self.solStep.newCptName and self.solStep.lastStep) and self.solStep

    def _getValueFieldKeys(self):
        keys = list(self.__dict__.keys())
        valueFiledKeys = []
        for key in keys:
            if key.find("val") >= 0 or key.find("Val") >= 0:
                valueFiledKeys.append(key)

        return valueFiledKeys

    def _addPrefixes(self):

        for key in self.valueFieldKeys:
            values = self.__dict__[key]
            if isinstance(values, list):
                for i in range(0, len(values)):
                    values[i] = self.prefixer.getSIPrefixedValue(values[i])
            if isinstance(values, dict):
                innerKeys = list(values.keys())
                for innerKey in innerKeys:
                    values[innerKey] = self.prefixer.getSIPrefixedValue(values[innerKey])
            else:
                self.__dict__[key] = self.prefixer.getSIPrefixedValue(values)
        return
