from lcapy import *

a = Circuit('VRC1step.sch')

v = a.R.v

ax = v.plot((-1, 10))

from matplotlib.pyplot import savefig
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')





