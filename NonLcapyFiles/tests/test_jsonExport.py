import os
import shutil

from lcapy.jsonExportCircuitInfo import JsonExportCircuitInfo
from lcapy.jsonExportCompValue import JsonCompValueExport
from lcapy.jsonExportVCValues import JsonVCValueExport
import NonLcapyFiles.solve as solve
import json
from lcapy import Circuit, FileToImpedance
from lcapy.solution import Solution


class TestJsonExport:

    def __del__(self):
        if os.path.isdir(".\\tempTest"):
            shutil.rmtree(".\\tempTest")

    @staticmethod
    def makeTestDir(folderName: str = "tempTest"):
        if not os.path.isdir(".\\tempTest"):
            os.mkdir(folderName)

    @staticmethod
    def readJson(filePath):
        with open(filePath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        return data

    def helperJsonExportCircuitInfo(self, fileName: str, values: list[str], compType: str):
        sol = solve.SolveInUserOrder(fileName, filePath="Schematics", savePath="/tempTest")
        filePath = sol.createCircuitInfo()
        data = self.readJson(filePath)
        for key in values:
            assert key in data.keys()

        assert data["componentTypes"] == compType

    def test_JasonExportCircuitInfo(self):
        self.makeTestDir()

        self.helperJsonExportCircuitInfo("R_parallel_dc",
                                         ["R1", "R2", "V1", "componentTypes", "omega_0"],
                                         "R")

        self.helperJsonExportCircuitInfo("R_parallel_ac",
                                         ["R1", "R2", "V1", "componentTypes", "omega_0"],
                                         "R")

        self.helperJsonExportCircuitInfo("L_parallel_ac",
                                         ["L1", "L2", "V1", "componentTypes", "omega_0"],
                                         "L")

        self.helperJsonExportCircuitInfo("C_parallel_ac",
                                         ["C1", "C2", "V1", "componentTypes", "omega_0"],
                                         "C")

        self.helperJsonExportCircuitInfo("Circuit_resistors_I",
                                         ["R1", "R2", "R3", "R4", "R5", "I1"],
                                         "R")

        self.helperJsonExportCircuitInfo("RC_series_ac",
                                         ["R1", "C2", "V1", "componentTypes", "omega_0"],
                                         "RLC")
