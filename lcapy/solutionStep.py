from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from lcapy.circuit import Circuit
    from lcapy import DrawWithSchemdraw

from lcapy.componentRelation import ComponentRelation
from typing import Union


class SolutionStep:
    def __init__(self, circuit: Circuit, cpt1: str, cpt2: str, newCptName: str,
                 relation: ComponentRelation, solutionText: str, lastStep: Circuit, nextStep: Circuit,
                 ):
        self.circuit: Circuit = circuit
        self.cpt1: Union[str, None] = cpt1
        self.cpt2: Union[str, None] = cpt2
        self.newCptName: Union[str, None] = newCptName
        self.relation: ComponentRelation = relation
        self.isInitialStep: bool = not (self.cpt1 or self.cpt2 or self.newCptName or self.relation)
        self.solutionText: Union[str, None] = solutionText
        self.lastStep: Union[Circuit, None] = lastStep
        self.nextStep: Union[Circuit, None] = nextStep

    def draw(self):
        DrawWithSchemdraw(self.circuit).draw()
