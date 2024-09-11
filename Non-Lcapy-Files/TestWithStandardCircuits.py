import warnings
import solve
import os
import argparse

filenames = ["Circuit_inductors.txt",
             "Circuit_resistors.txt",
             "Circuit_capacitors.txt",
             "Circuit_mixed_2pi30.txt",
             "Circuit_mixed_omega0.txt",
             "Circuit_mixed_30.txt",
             "Circuit_mixed.txt"]


def clearDir(path):
    if os.path.exists(path) and os.path.isdir(path):
        toRemove = os.listdir(path)
        for remove in toRemove:
            os.remove(os.path.join(path, remove))


def check_for_solutions(path):
    if len(os.listdir("Solutions")) == 0:
        raise AssertionError(f"{filename} did not produce an output")
    elif len(os.listdir("Solutions")) == 1:
        warnings.warn(f"{filename} only produced one solution")


def solveInUserOrder(filename):
    clearDir("Solutions")

    test = solve.SolveInUserOrder(filename=filename, savePath="Solutions/", filePath="StandardCircuits/")
    test.createInitialStep()

    test.simplifyTwoCpts(["Z4", "Z5"])
    test.simplifyTwoCpts(["Z1", "Zsim1"])
    test.simplifyTwoCpts(["Z2", "Z3"])
    test.simplifyTwoCpts(["Zsim2", "Zsim3"])


print("Try <solve.solve_circuit>")
for filename in filenames:
    clearDir("Solutions")
    solve.solve_circuit(filename, filePath="StandardCircuits/")
    check_for_solutions("Solutions")
    print(f"{filename+'...':40s}success")

print("Try <solve.SolveInUserOrder")
for filename in filenames:
    clearDir("Solutions")
    solveInUserOrder(filename)
    check_for_solutions("Solutions")
    print(f"{filename+'...':40s}success")

exit(0)
