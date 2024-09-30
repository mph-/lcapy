import pytest

from lcapy import Circuit
from lcapy import Solution
from lcapy.impedanceConverter import ValueToComponent, FileToImpedance, getSourcesFromCircuit, getOmegaFromCircuit
from os.path import join
from lcapy.componentnamer import ComponentNamer


class TestSimplifyStepwise:
    """tests ValueToComponent from impedanceConverter.py"""

    @staticmethod
    def assertResult(filename: str, compType: str, oneTroughImp=False, oneTroughReal=False,
                     path: str = "./Schematics"):
        ComponentNamer().reset()

        orgCct = Circuit(join(path, filename))
        orgVal1 = Solution.getElementSpecificValue(orgCct[compType+'1'])
        orgVal2 = Solution.getElementSpecificValue(orgCct[compType+'2'])

        cct = Circuit(FileToImpedance(join(path, filename)))
        sol = Solution(cct.simplify_stepwise())

        val1 = sol['step1'].lastStep.circuit.Z1.Z
        val2 = sol['step1'].lastStep.circuit.Z2.Z
        result = sol['step1'].circuit.Zsim1.Z
        if not oneTroughImp:
            assert val1 + val2 == result, "Calculation in Impedance is wrong"
        else:
            assert 1/(1/val1 + 1/val2) == result, "Calculation in Impedance is wrong"

        omega_0 = getOmegaFromCircuit(orgCct, getSourcesFromCircuit(orgCct))
        realRes, _ = ValueToComponent(result, omega_0)

        if not oneTroughReal:
            assert realRes == orgVal1 + orgVal2, "Transformed calculation is wrong"
        else:
            assert realRes == 1/(1/orgVal1.expr + 1/orgVal2.expr), "Transformed calculation is wrong"

    def test_R(self):
        # test for dc
        self.assertResult("R_series_dc", "R", False, False)
        self.assertResult("R_parallel_dc", "R", True, True)

        # test for ac
        self.assertResult("R_series_ac", "R", False, False)
        self.assertResult("R_parallel_ac", "R", True, True)

    def test_L(self):
        # test for dc
        self.assertResult("L_series_dc", "L", False, False)
        self.assertResult("L_parallel_dc", "L", True, True)

        # test for ac
        self.assertResult("L_series_ac", "L", False, False)
        self.assertResult("L_parallel_ac", "L", True, True)

    def test_C(self):
        # test for dc
        self.assertResult("C_series_dc", "C", False, True)
        self.assertResult("C_parallel_dc", "C", True, False)

        # test for ac
        self.assertResult("C_series_ac", "C", False, True)
        self.assertResult("C_parallel_ac", "C", True, False)

    def test_RL(self):
        # test for dc
        self.assertResult("RL_series_dc")
        self.assertResult("RL_parallel_dc", True)

        # test for ac
        self.assertResult("RL_series_ac")
        self.assertResult("RL_parallel_ac", True)

    def test_RC(self):
        # test for dc
        self.assertResult("RC_series_dc")
        self.assertResult("RC_parallel_dc", True)

        # test for ac
        self.assertResult("RC_series_ac")
        self.assertResult("RC_parallel_ac", True)

    def test_CL(self):
        # test for dc
        self.assertResult("CL_series_dc")
        self.assertResult("CL_parallel_dc", True)

        # test for ac
        self.assertResult("CL_series_ac")
        self.assertResult("CL_parallel_ac", True)
