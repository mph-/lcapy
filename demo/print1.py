from lcapy import Vdc, R, L, C, pprint
import numpy as np
from matplotlib.pyplot import figure, savefig, show

a = (Vdc(5) + L(10)) | C(1, 5)
b = a.load(R(5))


print('general')
pprint(b.V.general())

print('canonical')
pprint(b.V.canonical())

print('ZPK')
pprint(b.V.ZPK())

print('partfrac')
pprint(b.V.partfrac())

