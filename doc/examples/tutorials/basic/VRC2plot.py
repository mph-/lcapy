from lcapy import *

a = Circuit('VRC2.sch')

H = a.P1.transfer('P2')

v_i = voltage(sin(3 * t) * u(t))
V_i = v_i(s)
V_o = H * V_i
v_o = V_o(t)

ax = v_i.plot((-1, 10), label='input')
ax = v_o.plot((-1, 10), axes=ax, label='output')
ax.legend()

from matplotlib.pyplot import savefig
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')





