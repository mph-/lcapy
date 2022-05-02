from matplotlib.pyplot import savefig
from lcapy import Circuit, j, f, s, degrees, pi

b = Circuit('opamp-voltage-follower-RC-load-open-loop.sch')
G = b.transfer(2, 0, 4, 0)

f1 = 10
f2 = 1e6
A0 = 1e6

A = A0 * (1 / (1 + j * f / f1)) * (1 / (1 + j * f / f2))
A = A.subs(f, s / (j * 2 * pi))

Ro = 40
R = 20
C = 100e-9

G = G.subs({'A': A, 'Ro': Ro, 'R': R, 'C': C})

ax = G.bode_plot((1, 10e6), plot_type='dB-degrees')
ax[0].set_ylim(-30, 130)
ax[1].set_ylim(-240, 0)

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')
