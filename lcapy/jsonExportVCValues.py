import lcapy
from sympy import latex
from lcapy import state
from enum import Enum
from lcapy.netlistLine import NetlistLine
from lcapy import t
from lcapy.impedanceConverter import ValueToComponent
from lcapy.unitPrefixer import SIUnitPrefixer
from lcapy.unitWorkAround import UnitWorkAround as uwa


# ToDo create Enum ComponentRelation and use everywhere for component relations
class ComponentRelation(Enum):
    parallel = 'parallel'
    series = 'series'
    none = None

    def __eq__(self, other):
        if isinstance(other, str):
            return self.value == other
        elif isinstance(other, ComponentRelation):
            return self.value == other.value
        else:
            raise TypeError

    def to_string(self):
        return self.__str__()

    def __str__(self):
        return str(self.value)


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
    def __init__(self):
        # this class automatically prefixes every field that includes val or Val in the name and transforms it to
        # a latex string before exporting the dictionary
        self.circuit: 'lcapy.Circuit' = None
        self.simpCircuit: 'lcapy.Circuit' = None
        self.prefixer = SIUnitPrefixer()

        self.oldName = {'CompName': None, 'Uname': None, 'Iname': None}
        self.name1 = {'CompName': None, 'Uname': None, 'Iname': None}
        self.name2 = {'CompName': None, 'Uname': None, 'Iname': None}
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
            'oldName': list(self.oldName.values()),
            'name1': list(self.name1.values()),
            'name2': list(self.name2.values()),
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
                        as_dict[key][idx] = latex(val)
        except KeyError:
            raise AssertionError(f"A filed which name includes val or Val is not in the export dict. Key: {key}")

        return as_dict

    def _updateObjectValues(self, step: str, solution: 'lcapy.Solution'):
        self.solStep: 'lcapy.solutionStep' = solution[step]
        self.simpCircuit: 'lcapy.Circuit' = solution[step].circuit  # circuit with less elements (n elements)

        self.oldName['CompName'] = solution[step].newCptName
        self.name1['CompName'] = solution[step].cpt1
        self.name2['CompName'] = solution[step].cpt2

        if not self._isInitialStep():
            self.circuit: 'lcapy.Circuit' = solution[step].lastStep.circuit  # circuit with more elements (n+1 elements)

            self._updateUnames()
            self._updateInames()
            self._updateZandConvValues()
            self._updateU()
            self._updateI()
            self._addPrefixes()
            self._updateCompRel()
            self._updateEquations()

    def _updateUnames(self):
        self.oldName['Uname'] = 'U' + NetlistLine(str(self.simpCircuit[self.oldName['CompName']])).typeSuffix
        self.name1['Uname'] = 'U' + NetlistLine(str(self.circuit[self.name1['CompName']])).typeSuffix
        self.name2['Uname'] = 'U' + NetlistLine(str(self.circuit[self.name2['CompName']])).typeSuffix

    def _updateInames(self):
        self.oldName['Iname'] = 'I' + NetlistLine(str(self.simpCircuit[self.oldName['CompName']])).typeSuffix
        self.name1['Iname'] = 'I' + NetlistLine(str(self.circuit[self.name1['CompName']])).typeSuffix
        self.name2['Iname'] = 'I' + NetlistLine(str(self.circuit[self.name2['CompName']])).typeSuffix

    def _updateZandConvValues(self):
        self.oldValues['Z'], self.convOldValue = self._checkForConversion(self.simpCircuit[self.oldName['CompName']].Z)
        self.values1['Z'], self.convValue1 = self._checkForConversion(self.circuit[self.name1['CompName']].Z)
        self.values2['Z'], self.convValue2 = self._checkForConversion(self.circuit[self.name2['CompName']].Z)

    @staticmethod
    def _checkForConversion(value) -> tuple:
        convValue, convCompType = ValueToComponent(value)
        if 'Z' == convCompType:
            return value, None
        else:
            return value, uwa.addUnit(convValue, convCompType)

    def _updateU(self):
        self.oldValues['U'] = self.simpCircuit[self.oldName['CompName']].V(t)
        self.values1['U'] = self.circuit[self.name1['CompName']].V(t)
        self.values2['U'] = self.circuit[self.name2['CompName']].V(t)

    def _updateI(self):
        self.oldValues['I'] = self.simpCircuit[self.oldName['CompName']].I(t)
        self.values1['I'] = self.circuit[self.name1['CompName']].I(t)
        self.values2['I'] = self.circuit[self.name2['CompName']].I(t)

    def _updateCompRel(self):
        if self.solStep.relation == ComponentRelation.parallel:
            self.relation = ComponentRelation.parallel
        elif self.solStep.relation == ComponentRelation.series:
            self.relation = ComponentRelation.series
        else:
            self.relation = ComponentRelation.none

    @staticmethod
    def _makeLatexEquationU(valueZ, valueI, resultU) -> str:
        return f"{latex(valueZ)} • {latex(valueI)} = {latex(resultU)}"

    @staticmethod
    def _makeLatexEquationI(valueZ, valueU, resultI) -> str:
        return f"{latex(valueU)} / {latex(valueZ)} = {latex(resultI)}"

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
