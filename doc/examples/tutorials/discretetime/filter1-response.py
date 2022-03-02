from matplotlib.pyplot import savefig
from lcapy import Circuit, f

a = Circuit('filter1.sch')
H = a.transfer('P1', 'P2')
Hv = H.subs({'R1':22, 'C1':100e-9, 'R2':1e6, 'C2':1e-9})
Hv.bode_plot((0.1, 10e3))

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')
