from lcapy import Circuit
from lcapy import Solution
from lcapy.impedanceConverter import ValueToComponent, FileToImpedance
from os.path import join
from lcapy.componentnamer import ComponentNamer


class TestImpedanceConverter:
    """tests ValueToComponent from impedanceConverter.py"""
    @staticmethod
    def getCompTypeSolStep1(sol) -> str:
        """
        :param sol: Solution that has 'step1'
        :return: Z value of Zsim1 from sol['step1'].circuit
        """
        cir = sol['step1'].circuit
        value = cir.Zsim1.Z
        _, compType = ValueToComponent(value)
        return compType

    @staticmethod
    def assertType(filename: str, compType: str, path: str = "./Schematics"):
        ComponentNamer().reset()
        cct = Circuit(FileToImpedance(join(path, filename)))
        sol = Solution(cct.simplify_stepwise())
        assert TestImpedanceConverter.getCompTypeSolStep1(sol) == compType

    @staticmethod
    def assert2Types(filename: str, compType1: str, compType2: str, path: str = "./Schematics"):
        ComponentNamer().reset()
        cct = Circuit(FileToImpedance(join(path, filename)))
        sol = Solution(cct.simplify_stepwise())
        gotType = TestImpedanceConverter.getCompTypeSolStep1(sol)
        assert gotType == compType1 or gotType == compType2

    def test_R(self):
        # test for dc
        self.assertType("R_serie_dc", "R")
        self.assertType("R_parallel_dc", "R")

        # test for ac
        self.assertType("R_serie_ac", "R")
        self.assertType("R_parallel_ac", "R")

    def test_L(self):
        # test for dc
        self.assertType("L_serie_dc", "L")
        self.assertType("L_parallel_dc", "L")

        # test for ac
        self.assertType("L_serie_ac", "L")
        self.assertType("L_parallel_ac", "L")

    def test_C(self):
        # test for dc
        self.assertType("C_serie_dc", "C")
        self.assertType("C_parallel_dc", "C")

        # test for ac
        self.assertType("C_serie_ac", "C")
        self.assertType("C_parallel_ac", "C")

    def test_RL(self):
        # test for dc
        self.assertType("RL_serie_dc", "Z")
        self.assertType("RL_parallel_dc", "Z")

        # test for ac
        self.assertType("RL_serie_ac", "Z")
        self.assertType("RL_parallel_ac", "Z")

    def test_RC(self):
        # test for dc
        self.assertType("RC_serie_dc", "Z")
        self.assertType("RC_parallel_dc", "Z")

        # test for ac
        self.assertType("RC_serie_ac", "Z")
        self.assertType("RC_parallel_ac", "Z")

    def test_CL(self):
        # test for dc
        self.assert2Types("CL_serie_dc", "C", "L")
        self.assert2Types("CL_parallel_dc", "C", "L")

        # test for ac
        self.assert2Types("CL_serie_ac", "C", "L")
        self.assert2Types("CL_parallel_ac", "C", "L")
