from lcapy import *

alpha = symbol('alpha')
t0 = symbol('t0')
f0 = symbol('f0')
w0 = 2 * pi * f0

sigs = [cos(w0 * t), sin(w0 * t), exp(j * w0 * t),
        texpr(1), t, t**2, 1 / t, delta(t), delta(t - t0),
        H(t), t * H(t), sign(t),
        rect(t), sinc(t), tri(t), trap(t, alpha),
        exp(-abs(t)), exp(-t) * H(t)]

for sig in sigs:
    print(':math:`%s \\longleftrightarrow %s`\n' % (sig.latex(), sig(f).latex()))
