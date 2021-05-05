from lcapy import *

alpha = symbol('alpha')
m = symbol('m')
a = symbol('a')
f0 = symbol('f0')
w0 = 2 * pi * f0

sigs = [cos(w0 * n * dt), sin(w0 * n * dt), exp(j * w0 * n * dt),
        nexpr(1), delta(n), delta(n - m), a**n, a**-n, n * a**n, n * a**-n,
        H(n), exp(-n * dt) * H(n)]

for sig in sigs:
    print(':math:`%s \\longleftrightarrow %s`\n' % (sig.latex(), sig(z).latex()))
