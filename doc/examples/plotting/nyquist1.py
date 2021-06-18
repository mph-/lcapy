from matplotlib.pyplot import savefig
from lcapy import *

L = 10 * (s + 1) * (s + 2) / ((s - 3) * (s - 4))

L(f, causal=True).nyquist_plot((1e-6, 100), npoints=1000)

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')

