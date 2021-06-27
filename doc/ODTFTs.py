from lcapy import *

alpha = symbol('alpha')
m = symbol('m', integer=True)
W0 = symbol('Omega_0')
N = symbol('N')

sigs = [cos(W0 * n), sin(W0 * n), exp(j * W0 * n),
        nexpr(1), delta(n), delta(n - m),
        H(n), n * H(n), sign(n),
        alpha**-n * H(n), rect(n), rect(n / N), sinc(n)]

for sig in sigs:
    print(':math:`%s \\longleftrightarrow %s`\n' % (sig.latex(), sig(Omega).latex()))
