from matplotlib.pyplot import savefig
from lcapy import s

r1 = 500
i1 = 500

r2 = 5000
i2 = 5000

p1 = r1 + 1j * i1
p1c = r1 - 1j * i1

p2 = r2 + 1j * i2
p2c = r2 - 1j * i2

H = s**2 / (s + p1) / (s + p1c) / (s + p2) / (s + p2c)

H.plot()
savefig(__file__.replace('.py', '.png'), bbox_inches='tight')
