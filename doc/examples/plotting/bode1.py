from matplotlib.pyplot import savefig
from lcapy import *

H = 1 / (s**2 + 20*s + 10000)

H.bode_plot((1, 1000))

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')

