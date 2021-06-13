from lcapy import *

alpha = symbol('alpha')
m = symbol('m', integer=True)
f0 = symbol('f0')
w0 = 2 * pi * f0

sigs = [nexpr(1), delta(n), delta(n - m),
        H(n), n * H(n),
        alpha**-n * H(n)]

sigs2 = [n, n**2, 1 / n, rect(n), sinc(n), alpha**-abs(n), sign(n),
         cos(w0 * n * dt), sin(w0 * n * dt), exp(j * w0 * n * dt)]

for sig in sigs:
    print(':math:`%s \\longleftrightarrow %s`\n' % (sig.latex(), sig(k).latex()))
