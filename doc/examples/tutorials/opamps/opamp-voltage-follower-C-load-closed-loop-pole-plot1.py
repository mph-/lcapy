from matplotlib.pyplot import savefig
from lcapy import Circuit, j, f, s, degrees, pi

a = Circuit('opamp-voltage-follower-C-load.sch')
H = a.transfer(2, 0, 1, 0)

f1 = 10
f2 = 1e6
A0 = 1e6

A = A0 * (1 / (1 + j * f / f1)) * (1 / (1 + j * f / f2))
A = A.subs(f, s / (j * 2 * pi))

Ro = 40
C = 100e-9

H = H.subs({'A': A, 'Ro': Ro, 'C': C})

H.plot()

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')
