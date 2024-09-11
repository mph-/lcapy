import os.path
from enum import Enum

import sympy

import json
from lcapy import ConstantDomainExpression
from .solutionStep import SolutionStep
from ordered_set import OrderedSet
from typing import Iterable
from warnings import warn
from lcapy import state
from lcapy.mnacpts import R, L, C, Z
from lcapy import DrawWithSchemdraw
from lcapy.jsonExport import JsonExport
from lcapy.unitWorkAround import UnitWorkAround as uwa
from typing import Union


class Solution:
    def __init__(self, steps: list[tuple]):
        """
        This class simplifies the Solution and handles the access to all data that is necessary to create a step-by-step
        solution to a circuit. The input is of this class is the output of simplify_Stepwise
        :param steps:
        """
        self._attributes = {}
        self.available_steps = []
        self.mapKey = dict([("initialCircuit", "step0")])
        # convert the steps returned from simplify_stepwise to SolutionSteps
        # the simplify function cant return SolutionSteps because it imports lcapy and therefor results in a circular
        # import
        solSteps = []
        for step in steps:
            solSteps.append(SolutionStep(step))

        if not solSteps:
            raise AttributeError('can`t create Solution from empty list\n'
                                 'make sure parameter steps isn`t an empty list')

        # name and add first Circuit
        self.available_steps.append("step0")
        self.__setitem__("step0", solSteps[0])
        if len(solSteps) >= 2:
            self["step0"].nextStep = solSteps[1]
        else:
            warn("Solution only contains initial circuit and no simplification Steps")
            return

        for i in range(1, len(solSteps)):
            curStep = "step" + str(i)
            self.available_steps.append(curStep)

            self.__setitem__(curStep, solSteps[i])
            self[curStep].lastStep = solSteps[i - 1]

            # the list index is only to len(list) -1 accessible
            if i + 1 <= len(solSteps) - 1:
                self[curStep].nextStep = solSteps[i + 1]

        for solText, solStep in self._nextSolutionText(returnSolutionStep=True):
            solStep.solutionText = solText

    def __getitem__(self, key):
        try:
            return self._attributes[key]
        except KeyError:
            if key in self.mapKey.keys():
                return self._attributes[self.mapKey[key]]
            else:
                raise AttributeError

    def __setitem__(self, key, value):
        self._attributes[key] = value

    def __getattr__(self, key):
        try:
            return self._attributes[key]
        except KeyError:
            if key in self.mapKey.keys():
                return self._attributes[self.mapKey[key]]
            else:
                raise AttributeError

    def __setattr__(self, key, value):
        if key.startswith('_'):
            super().__setattr__(key, value)
        else:
            self.__setitem__(key, value)

    def addKeyMapping(self, accessKey, realKey: str):
        """
        If a KeyError is thrown a dictionary is searched. If it is in the dictionary the __getItem__ or __getAttr__
        it is tried again with the specified key in the Dictionary
        :param accessKey: key that you want to use
        :param realKey: key name from this class step0, step1, step2, ...
        :return: nothing / void
        """
        if realKey not in self.available_steps:
            raise KeyError("realKey doesn't exist")

        self.mapKey[accessKey] = realKey

    def removeKeyMapping(self, accessKey):
        if accessKey in self.mapKey.keys():
            self.mapKey.pop(accessKey)
        else:
            warn("accessKey not in mapKey")

    def getAvailableSteps(self, skip: set) -> OrderedSet[str]:
        """
        Returns all available steps
        :param skip: steps to skip
        :return: OrderedSet[str]
        """
        if skip:
            return OrderedSet(self.available_steps) - skip
        else:
            return OrderedSet(self.available_steps)

    @staticmethod
    def getElementSpecificValue(element: Union[R, C, L, Z], unit=False) -> ConstantDomainExpression:
        """
        accesses the value resistance, capacitance, inductance, or impedance of an element based on its type
        :param element: mnacpts.R | mnacpts.C | mnacpts.L | mnacpts.Z
        :param unit: if True the Unit (ohm, F, H) are added to the str
        :return: lcapy.ConstantDomainExpression
        """
        if unit:
            return uwa.addUnit(Solution.getElementSpecificValue(element), element.type)

        state.show_units = False
        if isinstance(element, R):
            returnVal = element.R
        elif isinstance(element, C):
            returnVal = element.C
        elif isinstance(element, L):
            returnVal = element.L
        elif isinstance(element, Z):
            returnVal = element.Z
        else:
            raise NotImplementedError(f"{type(element)} "
                                      f"not supported edit Solution.getElementSpecificValue() to support")

        return returnVal

    def solutionText(self, step: str) -> str:
        """
        returns the solution text of the given step
        :param step: step0, step1, step2, step<n> ..., getAvailableSteps returns all valid steps
        :return: solutionText for the given step
        """
        assert isinstance(self, Solution)
        assert isinstance(self[step], SolutionStep)

        if self[step].isInitialStep:
            return "\ninitialCircuit / step0"

        name1 = self[step].cpt1
        name2 = self[step].cpt2
        newName = self[step].newCptName
        thisStep = self[step]
        lastStep = self[step].lastStep

        solText = "\n-------------------------------"
        solText += f"\nSimplified {name1} and {name2} to {newName}"
        solText += f"\nthe components are in {thisStep.relation}"
        solText += f"\n{name1}: {self.getElementSpecificValue(lastStep.circuit[name1], unit=True)}"
        solText += f"\n{name2}: {self.getElementSpecificValue(lastStep.circuit[name2], unit=True)}"
        solText += (f"\n{newName} (Result):" +
                    f"{self.getElementSpecificValue(thisStep.circuit[newName], unit=True)}")
        return solText

    def _nextSolutionText(self, skip: set = None, returnSolutionStep: bool = False) -> str or (str, SolutionStep):
        """
        is used internally, Externally use steps() or completeSolutionText()
        :param skip: steps to skip
        :param returnSolutionStep: bool, default False when true returns str and SolutionStep else only str
        :return: str || str and SolutionStep, based on returnSolutionStep
        """
        for step in self.getAvailableSteps(skip):
            if not returnSolutionStep:
                yield self.solutionText(step)
            else:
                yield self.solutionText(step), self[step]

    def steps(self, skip: set = None) -> Iterable[SolutionStep]:
        """
        yields the steps from the simplification. They can be iterated in a for loop e.g.:
        for step in sol.Steps():
            print(step.solutionText)
            step.circuit.draw()
        :param skip: define a set of steps that should be skipped e.g. the step0 every other step is
        named as follows: step1, step2, ..., step<n>
        :return:
        """
        for step in self.getAvailableSteps(skip=skip):
            yield self[step]

    def completeSolutionText(self, skip: set = None) -> str:
        """
        returns the complete solution text for the Steps saved in the Solution Object it's called on
        :param skip: steps to skip
        :return: string that contains the complete solution text
        """
        solText = ""
        for step in self.getAvailableSteps(skip):
            solText += self.solutionText(step)
        return solText

    @staticmethod
    def check_path(path: str):
        if not os.path.isdir(path) and os.path.isfile(path):
            raise ValueError(f"{path} is a file not a directory")
        elif not os.path.isdir(path) and not path == "":
            os.mkdir(path)

    def draw(self, filename: str = "circuit", path: str = None):
        """
        saves a svg-File for each step in the Solution.
        can raise a value error if path is not a directory
        :param filename: optional filename, files will be named filename_step<n>.svg n = 0,1 ..., len(availableSteps)
        :param path: directory in which to save the json-File in, if None save in current directory
        :return: nothing
        """

        if path is None:
            path = ""

        Solution.check_path(path)
        filename = os.path.splitext(filename)[0]

        for step in self.available_steps:
            self.drawStep(step, filename=filename, path=path)

    def drawStep(self, step, filename=None, path: str = None) -> str:
        """
        draws the circuit for a specific step
        :param step: step0, step1, step2, ..., step<n> ..., self.getAvailableSteps returns all valid steps
        :param filename: optional filename, files will be named filename_step<n>.svg n = 0,1 ..., len(availableSteps)
        :param path: directory in which to save the json-File in, if None path is current directory
        :return: nothing
        """
        if path is None:
            path = ""

        if filename is None:
            filename = self.filename

        DrawWithSchemdraw(self[step].circuit, fileName=filename + f"_{step}.svg").draw(path=path)

        return os.path.join(path, filename + f"_{step}.svg")

    def exportStepAsJson(self, step, path: str = None, filename: str ="circuit", debug: bool = False) -> str:
        """
        saves a step as a .json File with the following information:
        name1 and name2 -> names of the simplified components
        newName -> name of the simplified component/ new component
        relation -> if the simplification was parallel or in series
        value1 and value2 -> the value of the component e.g. 10 ohm or 10 F ...
        result -> the value of the new Component
        hasConversion -> if a value is transformed to or from impedance
        convVal1, convVal2, convResult -> the value that was transformed from


        raises a value Error if information is missing in a step use try/except or when Path does not point to a file
        :param debug: if ture print the dictionary that is used for creating the json file
        :param step: a step name e.g. step0, step1, step2, ..., step<n>
        :param path:  path to save the json-File in if None save in current directory
        :param filename: svg-File will be named <filename>_step<n>.svg n = 0 | 1 | ...| len(availableSteps)
        :return: nothing
        """

        if path is None:
            path = ""

        Solution.check_path(path)
        filename = os.path.splitext(filename)[0]

        jsonExport = JsonExport()
        as_dict = jsonExport.getDictForStep(step, self)

        if debug:
            print(as_dict)

        fullPathName = os.path.join(path, filename) + "_" + step + ".json"
        with open(fullPathName, "w", encoding="utf-8") as f:
            json.dump(as_dict, f, ensure_ascii=False, indent=4)

        return fullPathName

    def export(self, path: str = None, filename: str = "circuit", debug: bool = False):
        """
        save a json-File for each step in available_steps.
        Files are named step<n> n = 0, 1 ..., len(availableSteps)
        can raise a value error, see exportStepAsJson for more information

        :param debug: print dictionary that is used to create the json-File
        :param path: directory in which to save the json-File in, if None save in current directory
        :param filename: svg-File will be named <filename>_step<n>.svg n = 0,1 ..., len(availableSteps)
        :return:
        """

        for step in self.available_steps:
            self.exportStepAsJson(step, path=path, filename=filename, debug=debug)
