from lcapy import (symbol, nexpr, delta, H, sign,
                   rect, sinc, sincu, n, sin, cos, exp, j, Omega)

alpha = symbol('alpha')
m = symbol('m', integer=True)
W0 = symbol('Omega_0')
No = symbol('N_o', odd=True)
Ne = symbol('N_e', even=True)
K = symbol('K')

sigs = [nexpr('x(n)'), nexpr('x(a * n)'), nexpr('x(n - m)'),
        cos(W0 * n), sin(W0 * n), exp(j * W0 * n),
        nexpr(1), delta(n), delta(n - m),
        H(n), n * H(n), sign(n),
        alpha**-n * H(n), rect(n), rect(n / No), rect(n / Ne),
        sinc(n), sinc(K * n), sinc(K * n)**2, sincu(K * n)]

for sig in sigs:
    print(':math:`%s \\longleftrightarrow %s`\n' %
          (sig.latex(), sig(Omega).latex()))
