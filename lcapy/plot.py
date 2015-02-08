import numpy as np
from matplotlib.pyplot import figure

def plot_pole_zero(obj, **kwargs):

    poles = obj.poles()
    zeros = obj.zeros()
    try:
        p = np.array([complex(p.evalf()) for p in poles.keys()])
        z = np.array([complex(z.evalf()) for z in zeros.keys()])
    except TypeError:
        raise TypeError('Cannot plot poles and zeros of symbolic expression')

    fig = figure()
    ax = fig.add_subplot(111)
    
    ax.axvline(0, color='0.7')
    ax.axhline(0, color='0.7')
    ax.axis('equal')
    ax.grid()

    a = np.hstack((p, z))
    xmin = a.real.min()
    xmax = a.real.max()
    ymin = a.imag.min()
    ymax = a.imag.max()

    xextra, yextra = 0.0, 0.0
    
    # This needs tweaking for better bounds.
    if len(a) >= 2:
        xextra, yextra = 0.1 * (xmax - xmin), 0.1 * (ymax - ymin)
    if xextra == 0:
        xextra += 1.0
    if yextra == 0:
        yextra += 1.0

    ax.set_xlim(xmin - 0.5 * xextra, xmax + 0.5 * xextra)
    ax.set_ylim(ymin - 0.5 * yextra, ymax + 0.5 * yextra)

    # TODO, annotate with number of times a pole or zero is repeated.
    ax.plot(z.real, z.imag, 'bo', fillstyle='none', ms=10, **kwargs)
    ax.plot(p.real, p.imag, 'bx', fillstyle='none', ms=10, **kwargs)
    

def plot_frequency(obj, t, **kwargs):

    from matplotlib.pyplot import figure
    
    # FIXME, determine useful frequency range...
    if f is None:
        f = (-0.2, 2)
    if isinstance(f, (int, float)):
        f = (0, f)
    if isinstance(f, tuple):
        f = np.linspace(f[0], f[1], 400)

    V = obj(f)

    # TODO, handle different formats; real/imag, mag/phase

    fig = figure()
    ax = fig.add_subplot(111)
    ax.plot(f, abs(V), **kwargs)
    ax.set_xlabel(obj.domain_label)
    ax.set_ylabel(obj.label)
    ax.grid(True)
    

def plot_time(obj, t, **kwargs):

    # FIXME, determine useful time range...
    if t is None:
        t = (-0.2, 2)
    if isinstance(t, (int, float)):
        t = (0, t)
    if isinstance(t, tuple):
        t = np.linspace(t[0], t[1], 400)

    v = obj(t)

    fig = figure()
    ax = fig.add_subplot(111)
    ax.plot(t, v, **kwargs)
    ax.set_xlabel(obj.domain_label)
    ax.set_ylabel(obj.label)
    ax.grid(True)




