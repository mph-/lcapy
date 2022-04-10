from lcapy import (symbol, nexpr, cos, sin, exp, j, sinc,
                   rect, sincu, n, dt, delta, sign, H, pi, f)

alpha = symbol('alpha')
m = symbol('m', integer=True)
f0 = symbol('f0')
w0 = 2 * pi * f0
No = symbol('N_o', odd=True)
Ne = symbol('N_e', even=True)
K = symbol('K')

sigs = [nexpr('x(n)'), nexpr('x(a * n)'), nexpr('x(n - m)'),
        cos(w0 * n * dt), sin(w0 * n * dt), exp(j * w0 * n * dt),
        nexpr(1), delta(n), delta(n - m),
        H(n), n * H(n), sign(n),
        alpha**-n * H(n), rect(n), rect(n / No), rect(n / Ne),
        sinc(n), sinc(K * n), sinc(K * n)**2, sincu(K * n)]

for sig in sigs:
    print(':math:`%s \\longleftrightarrow %s`\n' %
          (sig.latex(), sig(f).latex()))
