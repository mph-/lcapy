import lcapy
from sympy import latex
from lcapy.unitPrefixer import SIUnitPrefixer
from lcapy import omega0
from lcapy import resistance

lcapy.state.show_units = True

# I know this is a terrible test, but it does for now

prefixer = SIUnitPrefixer()
val = latex(prefixer.getSIPrefixedValue(10))
if not val == "10.0":
    raise AssertionError(f"value should be 10.0, is {val}")
print(f"{'10':20} -> {val+'...':30} passed")

val = latex(prefixer.getSIPrefixedValue(100))
if not val == "100.0":
    raise AssertionError(f"value should be 100.0, is {val}")
print(f"{'100':20} -> {val+'...':30} passed")

val = latex(prefixer.getSIPrefixedValue(1000))
if not val == "1.0 \\text{k}":
    raise AssertionError(f"value should be 1.0 \\text{{k}}, is {val}")
print(f"{'1000':20} -> {val+'...':30} passed")

val = latex(prefixer.getSIPrefixedValue(100.0))
if not val == "100.0":
    raise AssertionError(f"value should be 100.0, is {val}")
print(f"{'100.0':20} -> {val+'...':30} passed")

val = latex(prefixer.getSIPrefixedValue(1000.0))
if not val == "1.0 \\text{k}":
    raise AssertionError(f"value should be 1.0 \\text{{k}}, is {val}")
print(f"{'1000.0':20} -> {val+'...':30} passed")

val = latex(prefixer.getSIPrefixedValue(100*omega0))
if not val == "100 \\omega_{0}":
    raise AssertionError(f"value should be 100 \\omega_{{0}}, is {val}")
print(f"{'100*omgea0':20} -> {val+'...':30} passed")

val = latex(prefixer.getSIPrefixedValue(1000*omega0))
if not val == "1.0 \\omega_{0} \\text{k}":
    raise AssertionError(f"value should be 1.0 \\omega_{{0}} \\text{{k}}, is {val}")
print(f"{'1000*omgea0':20} -> {val+'...':30} passed")

val = latex(prefixer.getSIPrefixedValue(1000000*omega0))
if not val == "1.0 \\omega_{0} \\text{M}":
    raise AssertionError(f"value should be 1.0 \\omega_{{0}} \\text{{M}}, is {val}")
print(f"{'1000000*omgea0':20} -> {val+'...':30} passed")

val = latex(prefixer.getSIPrefixedValue(100000*omega0))
if not val == "0.1 \\omega_{0} \\text{M}":
    raise AssertionError(f"value should be 0.1 \\omega_{{0}} \\text{{M}}, is {val}")
print(f"{'100000*omega0':20} -> {val+'...':30} passed")

val = latex(prefixer.getSIPrefixedValue(100.0*omega0))
if not val == "100.0 \\omega_{0}":
    raise AssertionError(f"value should be 100.0 \\omega_{{0}}, is {val}")
print(f"{'100.0*omega0':20} -> {val+'...':30} passed")

val = latex(prefixer.getSIPrefixedValue(1000.0*omega0))
if not val == "1.0 \\omega_{0} \\text{k}":
    raise AssertionError(f"value should be 1.0 \\omega_{0} \\text{{k}}, is {val}")
print(f"{'1000.0*omega0':20} -> {val+'...':30} passed")

value = resistance(1000)
val = latex(prefixer.getSIPrefixedValue(value))
if not val == "1.0 \\text{k} \\Omega":
    raise AssertionError(f"value should be 1.0 \\text{{k}} \\Omega, is {val}")
print(f"{'1000 ohm':20} -> {val+'...':30} passed")