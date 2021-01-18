from lcapy import Circuit, f

a = Circuit(""" 
R 1 0; down 
W 1 2; right 
C 2 0_2; down 
W 0 0_2; right""")   

b = a.noisy()

Vn = b.C.V.n

Vns = Vn.subs({'R':10e3, 'C':1e-9, 'T':293, 'k_B':1.38e-23})   

from numpy import linspace
vf = linspace(0, 100e3, 200)

from matplotlib.pyplot import savefig
(Vns(f) * 1e9).plot(vf, plot_type='mag', ylabel='ASD (nV/rootHz)')

savefig(__file__.replace('.py', '.png'), bbox_inches='tight')

