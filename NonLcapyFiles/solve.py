# for lcapy version: 1.24+inskale.0.17
from lcapy import Circuit, FileToImpedance, DrawWithSchemdraw
from lcapy.solution import Solution
from lcapy.componentRelation import ComponentRelation
from lcapy.solutionStep import SolutionStep
import os
from lcapy.jsonExportCircuitInfo import JsonExportCircuitInfo


def solve_circuit(filename: str, filePath="Circuits/", savePath="Solutions/"):
    cct = Circuit(FileToImpedance(os.path.join(filePath, filename)))
    cct.namer.reset()
    steps = cct.simplify_stepwise()
    sol = Solution(steps)
    sol.draw(path=savePath, filename=filename)
    sol.export(path=savePath, filename=filename)


class SolveInUserOrder:
    def __init__(self, filename: str, filePath=None, savePath=None):
        """
        :param filename: str with filename of circuit to simplify
        :param filePath: str with path to circuit file if not in current directory
        :param savePath: str with path to save the result svg and jason files to
        """
        if filePath is None:
            filePath = ""
        if savePath is None:
            savePath = ""

        self.filename = filename
        self.filePath = filePath
        self.savePath = savePath
        self.circuit = Circuit(FileToImpedance(os.path.join(filePath, filename)))
        self.steps: list[SolutionStep] = [
            SolutionStep(self.circuit, None, None, None, None, None,
                         None, None)
        ]
        self.circuit.namer.reset()

        return

    def simplifyTwoCpts(self, cpts: list) -> tuple[bool, tuple[str, str], str]:
        """
        :param cpts: list with two component name strings to simplify ["R1", "R2"]
        :return tuple with bool if simplification is possible, str with json filename, str with svg filename
        """
        # ToDo this only works as long as only simplifiable components are selected which are represented as a
        # impedance internally in the cirucuit
        cpts[0] = "Z" + cpts[0][1::]
        cpts[1] = "Z" + cpts[1][1::]

        if cpts[1] in self.circuit.in_series(cpts[0]):
            newNet, newCptName = self.circuit.simplify_two_cpts(self.circuit, cpts)
            self.steps.append(SolutionStep(newNet, cpt1=cpts[0], cpt2=cpts[1], newCptName=newCptName,
                                           relation=ComponentRelation.series.value,
                                           solutionText=None, lastStep=None, nextStep=None))
        elif cpts[1] in self.circuit.in_parallel(cpts[0]):
            newNet, newCptName = self.circuit.simplify_two_cpts(self.circuit, cpts)
            self.steps.append(SolutionStep(newNet, cpt1=cpts[0], cpt2=cpts[1], newCptName=newCptName,
                                           relation=ComponentRelation.parallel.value,
                                           solutionText=None, lastStep=None, nextStep=None))
        else:
            return False, ("", ""), ""

        sol = Solution(self.steps)
        newestStep = sol.available_steps[-1]

        jsonName = sol.exportStepAsJson(newestStep, path=self.savePath, filename=os.path.splitext(self.filename)[0])
        svgName = sol.drawStep(newestStep, path=self.savePath, filename=os.path.splitext(self.filename)[0])

        self.circuit = newNet
        return True, jsonName, svgName

    def createInitialStep(self) -> tuple[bool, tuple[str, str], str]:
        """
        create the initial step or step0 of the circuit
        :return tuple with bool if simplification is possible, str with json filename, str with svg filename
        """

        sol = Solution(self.steps)
        nameStep0Json = sol.exportStepAsJson("step0", path=self.savePath, filename=self.filename)
        nameStep0Svg = sol.drawStep('step0', filename=self.filename, path=self.savePath)

        return True, nameStep0Json, nameStep0Svg

    def createCircuitInfo(self) -> str:
        sol = Solution(self.steps)
        return sol.exportCircuitInfo("step0", path=self.savePath, filename=self.filename)

    def createStep0(self) -> tuple[bool, tuple[str, str], str]:
        return self.createInitialStep()

    def getSolution(self):
        return Solution(self.steps)

