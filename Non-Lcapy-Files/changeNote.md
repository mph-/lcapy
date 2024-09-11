## 1.22+InskaLE.0.9
- Internal class NetlistLine uses lcapy parser instead of own implementation to be more robust

## 1.22+InskaLE.0.8
- write the value of a component (its resistance, capacitance ...) with its unit and unit prefix in latex
format into the svg under the tag value

## 1.22+InskaLE.0.7
- The step0-json now contains all components of the circuit in its initial state (without) simplifications
and the frequency omega_0 in Hz if it is an AC circuit. If it is a DC circuit omega_0 is equal to 0.
- simplification is done in impedance in all cases if the result can be transformed to R, L or C
the calculation and the saved values are transformed accordingly. If components are R, L or C and the
result is a component which is a mixture of those tree the calculation is displayed in impedance to the
user and the transformation of the components as well. Therefore, there are tree new values in the
export-Jason convVal1 convVal2 and convRes these are populated if on of the values value1 value2 or 
result values can be transformed. Also, there is an omega_0 value for AC analysis.
- Prefixes for values are added when the calculation is not in impedance. The Prefixes available are
available by importing: from sympy.physics.units.prefixes import PREFIXES
Values are round to max. 3 decimal places.

## 1.22+InskaLE.0.6
- show dots on branch knots, when there are more than two occurrences of the same node in the circuit.

## 1.22+InskaLE.0.5
- Fix error that the inverse sum is not used when needed.

## 1.22+InskaLE.0.4
- return filenames from draw and export functions. Add function to draw a standalone step.

## 1.22+InskaLE.0.3
- removed unnecessary print of latex expression

## 1.22+InskaLE.0.2
- The json-File is extended by latexEquation. The calculation from the simplification of two
elements. The expression is converted from sympy to a latex string

## V1.22+inskale.0.1
- Creates a step-by-step Solution where only to Elements are simplified
in one step and exports the Steps into a json-File. The Circuits
(more specific the Netlists) can be drawn with the additional package
Schemdraw and don't need LaTeX installation. Only Supported Elements are
Resistors, Inductors, Wires, Voltage Sources, Impedances (special form of Resistors)