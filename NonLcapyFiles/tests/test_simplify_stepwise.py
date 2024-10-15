import pytest
import sympy

import lcapy
from lcapy import Circuit
from lcapy.solution import Solution
from lcapy.impedanceConverter import StrToComponent, FileToImpedance, getSourcesFromCircuit, getOmegaFromCircuit
from os.path import join
from lcapy.componentnamer import ComponentNamer
from typing import Union


class TestSimplifyStepwise:
    """tests ValueToComponent from impedanceConverter.py"""

    @staticmethod
    def assertResult(filename: str, compType: str, oneTroughImp=False, oneTroughReal=False,
                     path: str = "./Schematics"):
        ComponentNamer().reset()

        orgCct = Circuit(join(path, filename))
        orgVal1 = Solution.getElementSpecificValue(orgCct[compType + '1'])
        orgVal2 = Solution.getElementSpecificValue(orgCct[compType + '2'])

        cct = Circuit(FileToImpedance(join(path, filename)))
        sol = Solution(cct.simplify_stepwise())

        val1 = sol['step1'].lastStep.circuit.Z1.Z
        val2 = sol['step1'].lastStep.circuit.Z2.Z
        result = sol['step1'].circuit.Zsim1.Z
        if not oneTroughImp:
            assert result == val1 + val2, "Calculation in Impedance is wrong"
        else:
            assert result == 1 / (1 / val1 + 1 / val2), "Calculation in Impedance is wrong"

        omega_0 = getOmegaFromCircuit(orgCct, getSourcesFromCircuit(orgCct))
        realRes, realResCompType = StrToComponent(result, omega_0)
        assert realResCompType == compType, "Converted type is not as expected"

        if not oneTroughReal:
            calcResult = orgVal1 + orgVal2
            if isinstance(calcResult, lcapy.Expr):
                assert realRes == pytest.approx(calcResult.expr), "Transformed calculation is wrong"
            else:
                assert realRes == pytest.approx(calcResult)
        else:
            calcResult = 1 / (1 / orgVal1.expr + 1 / orgVal2.expr)
            if isinstance(calcResult, lcapy.Expr):
                assert realRes == pytest.approx(calcResult.expr), "Transformed calculation is wrong"
            else:
                assert realRes == pytest.approx(calcResult)

    @staticmethod
    def assertResultImp(filename: str, compType: Union[str, None], oneTroughImp: bool = False,
                        path: str = "./Schematics"):
        ComponentNamer().reset()
        cct = Circuit(FileToImpedance(join(path, filename)))
        sol = Solution(cct.simplify_stepwise())

        val1 = sol['step1'].lastStep.circuit.Z1.Z
        val2 = sol['step1'].lastStep.circuit.Z2.Z
        result = sol['step1'].circuit.Zsim1.Z

        if compType:
            _, resultCompType = StrToComponent(result)
            assert resultCompType == compType, "Converted type is not as expected"

        if not oneTroughImp:
            assert result == val1 + val2, "Calculation in Impedance is wrong"
        else:
            assert result == 1 / (1 / val1 + 1 / val2), "Calculation in Impedance is wrong"

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
        self.assertResultImp("RL_series_dc", "Z", False)
        self.assertResultImp("RL_parallel_dc", "Z", True)

        # test for ac
        self.assertResultImp("RL_series_ac", "Z", False)
        self.assertResultImp("RL_parallel_ac", "Z", True)

    def test_RC(self):
        # test for dc
        self.assertResultImp("RC_series_dc", "Z", False)
        self.assertResultImp("RC_parallel_dc", "Z", True)

        # test for ac
        self.assertResultImp("RC_series_ac", "Z", False)
        self.assertResultImp("RC_parallel_ac", "Z", True)

    def test_CL(self):
        # test for dc
        self.assertResultImp("CL_series_dc", None, False)
        self.assertResultImp("CL_parallel_dc", None, True)

        # test for ac
        self.assertResultImp("CL_series_ac", None, False)
        self.assertResultImp("CL_parallel_ac", None, True)
