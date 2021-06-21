from lcapy import f, sqrt
In = 1e-12 * sqrt(250 / f + 1)
ax = (In * 1e12).plot((0.1, 10e3), loglog=True)
ax.grid(True, 'both')
ax.set_ylabel('Noise current density pA$/\sqrt{\mathrm{Hz}}$')
ax.set_ylim(0.1, 100)

from matplotlib.pyplot import savefig
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')
