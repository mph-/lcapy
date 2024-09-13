import lcapy
from sympy import latex
from lcapy import state
from enum import Enum
from lcapy.netlistLine import NetlistLine
from lcapy import t
from lcapy.impedanceConverter import ValueToComponent

# ToDo create Enum ComponentRelation and use everywhere for component relations
class ComponentRelation(Enum):
    parallel = 'parallel'
    series = 'series'
    none = None


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
        self.circuit: 'lcapy.Circuit' = None
        self.simpCircuit: 'lcapy.Circuit' = None

        self.oldName = {'CompName': None, 'Uname': None, 'Iname': None}
        self.name1 = {'CompName': None, 'Uname': None, 'Iname': None}
        self.name2 = {'CompName': None, 'Uname': None, 'Iname': None}
        self.oldValues = {'RLCZ': None, 'U': None, 'I': None}  # values to Component with oldName
        self.values1 = {'RLCZ': None, 'U': None, 'I': None}  # values to Component with name1
        self.values2 = {'RLCZ': None, 'U': None, 'I': None}  # values to Component with name2
        self.convOldValue = None  # converted Value of RLCZ from oldValues
        self.convValue1 = None  # converted Value of RLCZ from values 1
        self.convValue2 = None  # converted Value of RLCZ from values 2
        self.relation: ComponentRelation = ComponentRelation.none
        self.equation = {'CompName1': None, 'CompName2': None}

    def getDictForStep(self, step: str, solution: 'lcapy.Solution') -> dict:
        self._updateObjectValues(step, solution)
        return {}

    def _updateObjectValues(self, step: str, solution: 'lcapy.Solution'):
        self.circuit = solution[step].lastStep.circuit
        self.simpCircuit = solution[step].circuit

        self.oldName['CompName'] = solution[step].newCptName
        self.name1['CompName'] = solution[step].cpt1
        self.name2['CompName'] = solution[step].cpt2

        self._updateUnames()
        self._updateInames()
        self._updateRLCZandConvValues()
        self._updateU()
        self._updateI()

    def _updateUnames(self):
        self.oldName['Uname'] = 'U' + NetlistLine(str(self.circuit[self.oldName])).typeSuffix
        self.name1['Uname'] = 'U' + NetlistLine(str(self.circuit[self.name1])).typeSuffix
        self.name2['Uname'] = 'U' + NetlistLine(str(self.circuit[self.name2])).typeSuffix

    def _updateInames(self):
        self.oldName['Iname'] = 'I' + NetlistLine(str(self.circuit[self.oldName])).typeSuffix
        self.name1['Iname'] = 'I' + NetlistLine(str(self.circuit[self.name1])).typeSuffix
        self.name2['Iname'] = 'I' + NetlistLine(str(self.circuit[self.name2])).typeSuffix

    def _updateRLCZandConvValues(self):
        self.oldValues['RLCZ'], self.convOldValue = self._checkForConversion(self.circuit[self.oldName].Z)
        self.name1['RLCZ'] = self._checkForConversion(self.circuit[self.name1].Z)
        self.name2['RLCZ'] = self._checkForConversion(self.circuit[self.name2].Z)

    @staticmethod
    def _checkForConversion(value) -> tuple:
        convValue, convCompType = ValueToComponent(value)
        if 'Z' == convCompType:
            return value, None
        else:
            return convValue, value

    def _updateU(self):
        self.oldValues['U'] = self.circuit[self.oldName].V(t)
        self.name1['U'] = self.circuit[self.name1].V(t)
        self.name2['U'] = self.circuit[self.name2].V(t)

    def _updateI(self):
        self.oldValues['I'] = self.circuit[self.oldName].I(t)
        self.name1['I'] = self.circuit[self.name1].I(t)
        self.name2['I'] = self.circuit[self.name2].I(t)
