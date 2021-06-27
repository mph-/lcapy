from lcapy import *

alpha = symbol('alpha')
m = symbol('m', integer=True)
f0 = symbol('f0')
w0 = 2 * pi * f0
N = symbol('N', odd=True)
W = symbol('W')

sigs = [cos(w0 * n * dt), sin(w0 * n * dt), exp(j * w0 * n * dt),
        nexpr(1), delta(n), delta(n - m),
        H(n), n * H(n), sign(n),
        alpha**-n * H(n), rect(n), rect(n / N), sinc(n), sinc(W * n)]

for sig in sigs:
    print(':math:`%s \\longleftrightarrow %s`\n' % (sig.latex(), sig(f).latex()))
