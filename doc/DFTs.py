from lcapy import (symbol, nexpr, n, k, rect, sinc, sign, sin, cos, exp,
                   dt, j, pi, delta, H)

alpha = symbol('alpha')
m = symbol('m', integer=True)
N = symbol('N', integer=True, positive=True)
f0 = symbol('f0')
w0 = 2 * pi * f0

sigs = [nexpr('x(n)'), nexpr('x(a * n)'), nexpr('x(n - m)'),
        nexpr(1), delta(n), delta(n - m),
        H(n), n * H(n),
        alpha**-n * H(n),
        exp(j * 2 * pi * n / N),
        cos(2 * pi * n / N),
        sin(2 * pi * n / N)]

sigs2 = [n, n**2, 1 / n, rect(n), sinc(n), alpha**-abs(n), sign(n),
         cos(w0 * n * dt), sin(w0 * n * dt), exp(j * w0 * n * dt)]

for sig in sigs:
    print(':math:`%s \\longleftrightarrow %s`\n' %
          (sig.latex(), sig(k, N=N).simplify().latex()))
