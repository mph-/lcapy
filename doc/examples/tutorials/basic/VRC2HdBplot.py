from lcapy import *

a = Circuit('VRC2.sch')

H = a.P1.transfer('P2')

H(jw).dB.plot((0.1, 10), log_frequency=True)

from matplotlib.pyplot import savefig
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')





