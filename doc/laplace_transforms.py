from lcapy import (texpr, cos, sin, exp, delta, sign, t, H, ramp,
                   rampstep, rect, tri, symbol, s, pi, j, Derivative,
                   Integral)

alpha = symbol('alpha')
t0 = symbol('t0')
f0 = symbol('f0')
w0 = 2 * pi * f0

sigs = [texpr('x(t)'), texpr('x(a * t)'), texpr('x(t - tau)'),
        cos(w0 * t), sin(w0 * t), exp(j * w0 * t),
        texpr(1), t, t**2, delta(t), delta(t - t0),
        H(t), t * H(t), sign(t),
        exp(-abs(t)), exp(-t) * H(t),
        rect(t - 0.5), tri(t - 1), ramp(t), rampstep(t),
        Derivative('x(t)', t), Derivative('x(t)', t, 2),
        Integral('x(t)', (t, 0, t))]

for sig in sigs:
    print(':math:`%s \\longleftrightarrow %s`\n' %
          (sig.latex(), sig(s).latex()))
