from lcapy import Circuit, t, s
a = Circuit("circuit-RLC-ivp1.sch")
a.R.V(s)
VR, defs = a.R.V(s).parameterize(zeta=False)
VR
defs

