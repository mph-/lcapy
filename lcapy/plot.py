"""
This module performs plotting using matplotlib.

Copyright 2014--2018 Michael Hayes, UCECE
"""

import numpy as np
from matplotlib.pyplot import figure

# Perhaps add Formatter classes that will produce the plot data?


def plot_pole_zero(obj, **kwargs):

    poles = obj.poles()
    zeros = obj.zeros()
    try:
        p = np.array([complex(p.evalf()) for p in poles.keys()])
        z = np.array([complex(z.evalf()) for z in zeros.keys()])
    except TypeError:
        raise TypeError('Cannot plot poles and zeros of symbolic expression')


    ax = kwargs.pop('axes', None)
    if ax is None:
        figsize = kwargs.pop('figsize', (8, 6))        
        fig = figure(figsize=figsize)
        ax = fig.add_subplot(111)
    ax.axvline(0, color='0.7')
    ax.axhline(0, color='0.7')
    ax.axis('equal')
    ax.grid()

    a = np.hstack((p, z))
    x_min = a.real.min()
    x_max = a.real.max()
    y_min = a.imag.min()
    y_max = a.imag.max()

    x_extra, y_extra = 0.0, 0.0

    # This needs tweaking for better bounds.
    if len(a) >= 2:
        x_extra, y_extra = 0.1 * (x_max - x_min), 0.1 * (y_max - y_min)
    if x_extra == 0:
        x_extra += 1.0
    if y_extra == 0:
        y_extra += 1.0

    ax.set_xlim(x_min - 0.5 * x_extra, x_max + 0.5 * x_extra)
    ax.set_ylim(y_min - 0.5 * y_extra, y_max + 0.5 * y_extra)

    # TODO, annotate with number of times a pole or zero is repeated.

    # Marker size
    ms = kwargs.pop('ms', 10)
    fillstyle = kwargs.pop('fillstyle', 'none')
    ax.plot(z.real, z.imag, 'bo', fillstyle=fillstyle, ms=ms, **kwargs)
    ax.plot(p.real, p.imag, 'bx', fillstyle=fillstyle, ms=ms, **kwargs)
    return ax


def plot_frequency(obj, f, **kwargs):

    from matplotlib.pyplot import figure

    # FIXME, determine useful frequency range...
    if f is None:
        f = (0, 2)
    if isinstance(f, (int, float)):
        f = (0, f)
    if isinstance(f, tuple):
        f = np.linspace(f[0], f[1], 400)

    if not hasattr(obj, 'part'):
        obj2 = None        
        plot_type = kwargs.pop('plot_type', 'dB_phase')
        if plot_type == 'dB_phase':
            obj1 = obj.magnitude.dB
            if obj.is_complex:
                obj2 = obj.phase
        elif plot_type == 'mag_phase':
            obj1 = obj.magnitude
            if not obj.is_positive:
                obj2 = obj.phase
        elif plot_type == 'real_imag':
            obj1 = obj.real
            obj2 = obj.imag
        else:
            raise ValueError('Unknown plot type: %s' % plot_type)

        if obj2 is None:
            return plot_frequency(obj1, f, **kwargs)
        
        ax = plot_frequency(obj1, f, **kwargs)
        ax2 = ax.twinx()
        kwargs['axes'] = ax2
        kwargs['linestyle'] = '--'
        ax2 = plot_frequency(obj2, f, **kwargs)
        return ax, ax2

    ax = kwargs.pop('axes', None)
    if ax is None:
        figsize = kwargs.pop('figsize', (8, 6))
        fig = figure(figsize=figsize)        
        ax = fig.add_subplot(111)

    V = obj.evaluate(f)

    log_magnitude = kwargs.pop('log_magnitude', False)
    log_frequency = kwargs.pop('log_frequency', False) or kwargs.pop('log_scale', False)
    if kwargs.pop('loglog', False):
        log_magnitude = True 
        log_frequency = True    

    plots = {(True, True) : ax.loglog,
             (True, False) : ax.semilogy,
             (False, True) : ax.semilogx,
             (False, False) : ax.plot}
    
    if obj.part == 'magnitude':    
        plot = plots[(log_magnitude, log_frequency)]
    else:
        plot = plots[(False, log_frequency)]                    

    xlabel = kwargs.pop('xlabel', obj.domain_label)
    ylabel = kwargs.pop('ylabel', obj.label)        
    plot(f, V, **kwargs)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)    
    ax.grid(True)
    return ax


def plot_angular_frequency(obj, omega, **kwargs):

    from matplotlib.pyplot import figure

    # FIXME, determine useful frequency range...
    if omega is None:
        omega = (0, np.pi)
    if isinstance(omega, (int, float)):
        omega = (0, omega)
    if isinstance(omega, tuple):
        omega = np.linspace(omega[0], omega[1], 400)

    return plot_frequency(obj, omega, **kwargs)


def plot_time(obj, t, **kwargs):

    # FIXME, determine useful time range...
    if t is None:
        t = (-0.2, 2)
    if isinstance(t, (int, float)):
        t = (0, t)
    if isinstance(t, tuple):
        t = np.linspace(t[0], t[1], 400)

    v = obj.evaluate(t)

    ax = kwargs.pop('axes', None)
    if ax is None:
        figsize = kwargs.pop('figsize', (8, 6))
        fig = figure(figsize=figsize)
        ax = fig.add_subplot(111)
    xlabel = kwargs.pop('xlabel', obj.domain_label)
    ylabel = kwargs.pop('ylabel', obj.label)
    ax.plot(t, v, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True)
    return ax
