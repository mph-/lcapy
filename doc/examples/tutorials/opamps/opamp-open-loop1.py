from matplotlib.pyplot import savefig
from lcapy import j, f, degrees, pi

f1 = 10
f2 = 1e6
A0 = 1e6

A = A0 * (1 / (1 + j * f / f1)) * (1 / (1 + j * f / f2))

ax = A.bode_plot((1, 10e6), plot_type='dB-degrees')
ax[0].set_ylim(-30, 130)
ax[1].set_ylim(-240, 0)

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')
