from matplotlib.pyplot import savefig
from lcapy import *

H = 10 * (s + 1) * (s + 2) / ((s - 3) * (s - 4))

H.nyquist_plot((-100, 100))

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')

