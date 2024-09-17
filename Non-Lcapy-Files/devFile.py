import warnings
import solve
import os
from lcapy import Circuit
from lcapy import t
from lcapy import state
from sympy import latex
import sympy

from lcapy import units

state.show_units = False

cct = Circuit('''V1 0 1 ac {10} {0} {100}; up
W 1 2; right
R1 2 3 {1000}; down
W 3 0; left
''')
voltage = cct.R1.V(t)
current = cct.R1.I(t)
voltageUnit = cct.R1.V(t).units
currentUnit = cct.R1.I(t).units
state.show_units = True
print("DevVoltage: " + str(latex(voltage.expr_with_units)))
print("DevCurrent: " + str(latex(current.expr_with_units)))

# exit("Process finished with exit code 0 - early Exit with exit statement")

filenames = ["Circuit_inductors.txt",  # 0
             "Circuit_resistors.txt",  # 1
             "Circuit_capacitors.txt",  # 2
             "Circuit_mixed_2pi30.txt",  # 3
             "Circuit_mixed_omega0.txt",  # 4
             "Circuit_mixed_30.txt",  # 5
             "Circuit_mixed.txt"]  # 6

solve.solve_circuit(filenames[1], filePath="StandardCircuits")

# cct = Circuit("StandardCircuits/Circuit_resistors.txt")
# print(cct.I(cct.R1))
# print(cct['R1'].I)
# print(cct.R1.I)
# print(cct.R1.V)
