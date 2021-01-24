from lcapy import *

a = Circuit('VRC1.sch')

v = a.R.v.subs(omega0, 3)

ax = v.plot((-1, 10))

from matplotlib.pyplot import savefig
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')





