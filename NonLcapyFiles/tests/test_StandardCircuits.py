from NonLcapyFiles.solve import solve_circuit, SolveInUserOrder
import os
import random

# string is filename, integer is number of steps that shall be created
# the initial step has to be included
filenames = [("Circuit_inductors.txt", 5),
             ("Circuit_resistors.txt", 5),
             ("Circuit_capacitors.txt", 5),
             ("Circuit_mixed_2pi30.txt", 5),
             ("Circuit_mixed_omega0.txt", 5),
             ("Circuit_mixed_30.txt", 5),
             ("Circuit_mixed.txt", 5),
             ("Circuit_resistors_I", 5),
             ("Circuit_resistors_I_ac", 5)]


def clearDir(path):
    if os.path.exists(path) and os.path.isdir(path):
        toRemove = os.listdir(path)
        for remove in toRemove:
            os.remove(os.path.join(path, remove))


def check_for_solutions(solSteps: int, filename: str, path: str = "../Solutions", filesPerStep: int = 3):
    # each step produces 2 json files and 1 svg file
    assert len(os.listdir(path)) == solSteps*filesPerStep, f"{filename} didn't produce as many files as expected"


def solveInUserOrder(filename):
    clearDir("../Solutions")

    test = SolveInUserOrder(filename=filename, savePath="../Solutions/", filePath="../StandardCircuits/")
    test.createInitialStep()

    test.simplifyTwoCpts(["Z4", "Z5"])
    test.simplifyTwoCpts(["Z1", "Zsim1"])
    test.simplifyTwoCpts(["Z2", "Z3"])
    test.simplifyTwoCpts(["Zsim2", "Zsim3"])


def test_solve_circuits():
    for filename in filenames:
        clearDir("../Solutions")
        solve_circuit(filename[0], filePath="../StandardCircuits/", savePath="../Solutions")
        check_for_solutions(filename=filename[0], solSteps=filename[1])
        print(f"{filename[0]+'...':40s}success")


def test_SolveInUserOrder():
    for filename in filenames:
        clearDir("../Solutions")
        solveInUserOrder(filename[0])
        check_for_solutions(filename=filename[0], solSteps=filename[1])
        print(f"{filename[0]+'...':40s}success")