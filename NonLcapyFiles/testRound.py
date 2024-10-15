
import sympy
from sympy.physics.units.prefixes import PREFIXES, Prefix
from sympy.printing import latex
from sympy import Float
from sympy import Rational

pref = PREFIXES['k']
a = sympy.Mul(11.3333333333) * pref

b = Rational(1133, 100)

for floatVal in list(a.atoms(Float)):
    expr = a.evalf(subs={floatVal: str(round(floatVal, 3))})
    print(f"Atoms: {a.atoms()} | floatVal: {floatVal:.10f} | newExpr: {expr} | roundedVal: {round(floatVal, 3):.10f}")
    print(latex(expr))
    print(latex(b))

exit()
