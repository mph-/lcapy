from lcapy import Circuit, sqrt, f, oo

Rs = 30
G = 1000

R1 = 100
R2 = (G - 1) * R1

Vn = 1e-9 * sqrt(3.5 / f + 1)
In = 1e-12 * sqrt(250 / f + 1)

T = 273 + 20
k_B = 1.38e-23

a = Circuit('opamp-noninverting-amplifier.sch')

an = a.noisy()

Vno = an[8].V.n(f)

Vno = Vno.limit('A', oo)

Vnov = Vno.subs({'R1':R1, 'R2':R2, 'In1':0, 'In2':0, 'Vn':Vn, 'Rs':Rs,
                 'k_B':k_B, 'T':0})

Vnoi = Vno.subs({'R1':R1, 'R2':R2, 'In1':In, 'In2':In, 'Vn':0, 'Rs':Rs,
                 'k_B':k_B, 'T':0})

Vnor = Vno.subs({'R1':R1, 'R2':R2, 'In1':0, 'In2':0, 'Vn':0, 'Rs':Rs,
                 'k_B':k_B, 'T':T})

Vnot = Vno.subs({'R1':R1, 'R2':R2, 'In1':In, 'In2':In, 'Vn':Vn, 'Rs':Rs,
                 'k_B':k_B, 'T':T})

flim = (0.1, 10e3)
ax = (Vnot * 1e9).plot(flim, loglog=True, label='total')
ax = (Vnov * 1e9).plot(flim, loglog=True, label='Vn', axes=ax)
ax = (Vnoi * 1e9).plot(flim, loglog=True, label='In', axes=ax)
ax = (Vnor * 1e9).plot(flim, loglog=True, label='R', axes=ax)
ax.set_ylabel('Noise voltage density nV/$\sqrt{\mathrm{Hz}}$')
ax.grid(True, 'both')
ax.legend()

from matplotlib.pyplot import savefig
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')
