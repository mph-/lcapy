from lcapy.jsonExportBase import JsonExportBase
from sympy.physics.units import Hz
from lcapy import Circuit
from sympy import parse_expr, pi, Mul
from lcapy import omega0
from lcapy.impedanceConverter import ComponentToImpedance, ValueToComponent
from lcapy.unitWorkAround import UnitWorkAround as uwa
from lcapy.netlistLine import NetlistLine


class JsonExportCircuitInfo(JsonExportBase):
    def __init__(self, precision=3):
        super().__init__(precision)
        self.omega_0 = None

    def getDictForStep(self, step, solution: 'lcapy.Solution') -> dict:
        as_dict = {}
        compTypes = set()

        for cptName in solution[step].circuit.elements.keys():
            cpt = solution[step].circuit.elements[cptName]
            if cpt.type == "V" or cpt.type == "I":
                if cpt.type == "V":
                    as_dict[cptName] = self.latexWithPrefix(cpt.v.expr_with_units)
                else:
                    as_dict[cptName] = self.latexWithPrefix(cpt.i.expr_with_units)

                if cpt.has_ac:
                    if cpt.args[2] is not None:
                        as_dict["omega_0"] = self.latexWithPrefix(
                            parse_expr(str(cpt.args[2]), local_dict={"pi": pi}) * Hz
                        )
                        try:
                            self.omega_0 = float(cpt.args[2])
                        except ValueError:
                            self.omega_0 = str(cpt.args[2])
                    else:
                        as_dict["omega_0"] = self.latexWithPrefix(omega0)
                        self.omega_0 = "omega_0"
                elif cpt.has_dc:
                    as_dict["omega_0"] = self.latexWithPrefix(Mul(0) * Hz)
                else:
                    raise AssertionError("Voltage Source is not ac or dc")

            elif not cpt.type == "W":
                value, compType = ValueToComponent(cpt.Z)
                compTypes.add(compType)
                as_dict[compType + cpt.id] = self.latexWithPrefix(uwa.addUnit(value, compType))

        if len(compTypes) == 1:
            if "R" in compTypes:
                as_dict["componentTypes"] = "R"
            elif "L" in compTypes:
                as_dict["componentTypes"] = "L"
            elif "C" in compTypes:
                as_dict["componentTypes"] = "C"
            else:
                raise ValueError("Unexpected type in set types")
        else:
            as_dict["componentTypes"] = "RLC"

        return as_dict
