from matplotlib.pyplot import savefig

from lcapy import f, sqrt
V = 1e-9 * sqrt(3.5 / f + 1)
ax = (V * 1e9).plot((-1, 4), loglog=True)
ax.grid(True, 'both')
ax.set_ylabel('Noise voltage density nV$\sqrt{\mathrm{Hz}}$')
ax.set_ylim(0.1, 10)

from matplotlib.pyplot import savefig
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')
