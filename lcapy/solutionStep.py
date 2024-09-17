from lcapy import Circuit
from lcapy import DrawWithSchemdraw
from lcapy.componentRelation import ComponentRelation
from typing import Union


class SolutionStep:
    def __init__(self, step: tuple):
        self.circuit: Circuit = step[0]
        self.cpt1: Union[str, None] = step[1]
        self.cpt2: Union[str, None] = step[2]
        self.newCptName: Union[str, None] = step[3]
        self.relation: ComponentRelation = step[4]
        self.isInitialStep: bool = not (step[1] or step[2] or step[3] or step[4])
        self.solutionText: Union[str, None] = None
        self.lastStep: Union[Circuit, None] = None
        self.nextStep: Union[Circuit, None] = None

    def draw(self):
        DrawWithSchemdraw(self.circuit).draw()
