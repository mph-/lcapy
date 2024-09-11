from lcapy import Circuit
from lcapy import DrawWithSchemdraw


class SolutionStep:
    def __init__(self, step: tuple):
        self.circuit = step[0]
        self.cpt1 = step[1]
        self.cpt2 = step[2]
        self.newCptName = step[3]
        self.relation = step[4]
        self.isInitialStep = not (step[1] or step[2] or step[3] or step[4])
        self.solutionText = None
        self.lastStep = None
        self.nextStep = None

    def draw(self):
        DrawWithSchemdraw(self.circuit).draw()
