from lcapy import *

a = Circuit('VRC2.sch')

H = a.P1.transfer('P2')

H(jw).magnitude.plot((0, 10))

from matplotlib.pyplot import savefig
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')





