from lcapy import symbol, nexpr, delta, cos, sin, exp, H, n, dt, pi, j, z

alpha = symbol('alpha')
p = symbol('p')
a = symbol('a')
f0 = symbol('f0')
w0 = 2 * pi * f0

sigs = [nexpr('x(n)'), nexpr('x(n - m)'),
        cos(w0 * n * dt), sin(w0 * n * dt), exp(j * w0 * n * dt),
        nexpr(1), delta(n), delta(n - p), a**n, a**-n, n * a**n, n * a**-n,
        H(n), exp(-n * dt) * H(n)]

for sig in sigs:
    print(':math:`%s \\longleftrightarrow %s`\n' %
          (sig.latex(), sig(z).latex()))
