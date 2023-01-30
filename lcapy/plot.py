"""
This module performs plotting using matplotlib.

Copyright 2014--2022 Michael Hayes, UCECE
"""

import numpy as np
from .symbols import pi
from .utils import separate_dirac_delta, factor_const
from sympy import DiracDelta, solve
from warnings import warn


# Perhaps add Formatter classes that will produce the plot data?


def make_axes(figsize=None, axes=None, **kwargs):

    from matplotlib.pyplot import subplots

    if axes is not None:
        if isinstance(axes, tuple):
            # FIXME
            axes = axes[0]
        fig = axes.figure
    elif figsize is not None:
        fig, axes = subplots(1, figsize=figsize, **kwargs)
    else:
        fig, axes = subplots(1, **kwargs)

    return axes


def parse_range(frange, zero=0.1, positive=False):

    # Ensure that have floats.  This deals with SymPy symbols.
    fmin = float(frange[0])
    fmax = float(frange[1])

    if fmax < fmin:
        raise ValueError('Max value smaller than min value')

    if positive and fmin < 0:
        raise ValueError('Min value not positive.')

    if positive and fmin == 0:
        fmin = zero

    return fmin, fmax


def plot_deltas(ax, t, deltas, var, plot_type='real', color='k'):

    vals = []
    for delta in deltas:
        delta = delta.expand(diracdelta=True, wrt=var)
        const, expr = factor_const(delta, var)
        if not (expr.is_Function and expr.func is DiracDelta):
            raise ValueError('Cannot handle %s' % expr)
        t0 = solve(expr.args[0], var)[0]

        const = complex(const)
        if plot_type in ('real', 'default'):
            const = const.real
        elif plot_type == 'imag':
            const = const.imag
        elif plot_type == 'abs':
            const = abs(const)
        elif plot_type == 'dB':
            const = 20 * np.log10(abs(const))
        elif plot_type == 'phase':
            const = np.angle(const)
        elif plot_type == 'radians':
            const = np.angle(const)
        elif plot_type == 'degrees':
            const = np.degrees(np.angle(const))
        else:
            raise ValueError('Unhandled plot type %s' % plot_type)

        if t0 >= min(t) and t0 <= max(t):
            vals.append((t0, const))

    if vals == []:
        return

    vals = np.array(vals)

    tmin, tmax = ax.get_xlim()
    T = tmax - tmin
    ymin, ymax = ax.get_ylim()
    yvalsmax = vals[:, 1].max()

    for t0, y in vals:

        # TODO: determine best scaling
        y1 = y / yvalsmax * ymax / 2

        st = T / 10
        sy = ymax

        ax.arrow(t0, 0, 0, y1, fc=color, ec=color,
                 width=0.025 * st, head_width=0.1 * st,
                 head_length=0.1 * sy, overhang=0,
                 length_includes_head=True, clip_on=False)

        # TODO: nicely format the scale factor
        ax.text(t0 + T / 100, y1 + ymax / 100, '%.*g' % (3, y))


def plot_pole_zero(obj, **kwargs):
    """This is a helper function for a pole zero plot.  It is better
    to use the `plot()` method of a LaplaceDomainExpression."""

    from matplotlib.pyplot import Circle, rcParams

    obj = obj.doit()

    poles = obj.poles()
    zeros = obj.zeros()
    try:
        p = np.array([p.cval for p in poles.keys()])
        z = np.array([z.cval for z in zeros.keys()])
    except TypeError:
        raise TypeError('Cannot plot poles and zeros of symbolic expression')

    ax = make_axes(figsize=kwargs.pop('figsize', None),
                   axes=kwargs.pop('axes', None))

    unitcircle = kwargs.pop('unitcircle', False)

    ax.axvline(0, color='0.7')
    ax.axhline(0, color='0.7')

    a = np.hstack((p, z))

    if unitcircle:
        ax.add_artist(Circle((0, 0), 1, color='blue',
                      linestyle='--', fill=False))

    x_min = a.real.min()
    x_max = a.real.max()
    y_min = a.imag.min()
    y_max = a.imag.max()

    if unitcircle:
        if x_min > -1:
            x_min = -1
        if x_max < 1:
            x_max = 1
        if y_min > -1:
            y_min = -1
        if y_max < 1:
            y_max = 1

    x_extra, y_extra = 0.0, 0.0

    # This needs tweaking for better bounds.
    if len(a) >= 2:
        x_extra, y_extra = 0.1 * (x_max - x_min), 0.1 * (y_max - y_min)
    if x_extra == 0:
        x_extra += 1.0
    if y_extra == 0:
        y_extra += 1.0

    x_min -= 0.5 * x_extra
    x_max += 0.5 * x_extra
    if unitcircle:
        bbox = ax.get_window_extent()
        aspect = bbox.width / bbox.height

        x_min *= aspect
        x_max *= aspect

    ax.axis('equal')
    ax.set_xlim(x_min, x_max)
    # overconstrained so ignored
    #ax.set_ylim(y_min - 0.5 * y_extra, y_max + 0.5 * y_extra)

    def annotate(axes, poles, offset=None):
        if offset is None:
            xmin, xmax = axes.get_xlim()
            offset = (xmax - xmin) / 40

        for pole, num in poles.items():
            if num > 1:
                p = pole.cval
                axes.text(p.real + offset, p.imag + offset, '%d' % num)

    # Marker size.  Unfortunately, the default is too small
    # but if the user wants a size of 6 set in rcParams, they are out of luck...
    if 'ms' not in kwargs and 'markersize' not in kwargs and rcParams['lines.markersize'] == 6:
        kwargs['ms'] = 19

    fillstyle = kwargs.pop('fillstyle', 'none')
    xlabel = kwargs.pop('xlabel', 'Re(%s)' % obj.var)
    ylabel = kwargs.pop('ylabel', 'Im(%s)' % obj.var)
    title = kwargs.pop('title', None)

    lines = ax.plot(z.real, z.imag, 'o', fillstyle=fillstyle, **kwargs)
    annotate(ax, zeros)
    color = kwargs.pop('color', lines[0].get_color())
    ax.plot(p.real, p.imag, 'x', fillstyle=fillstyle, color=color, **kwargs)
    annotate(ax, poles)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

    ax.grid(True)

    return ax


def plotit(ax, obj, f, V, plot_type=None, deltas=None, log_magnitude=False,
           log_frequency=False, norm=False, dbmin=-120, unwrap=False, **kwargs):

    plots = {(True, True): ax.loglog,
             (True, False): ax.semilogy,
             (False, True): ax.semilogx,
             (False, False): ax.plot}

    label = ''
    if hasattr(obj, 'quantity_label'):
        label = obj.quantity_label
    part = obj.part
    units = str(obj.units)

    if plot_type == 'real':
        V = V.real
        part = 'real part'
    elif plot_type == 'imag':
        V = V.imag
        part = 'imag part'
    elif plot_type == 'magnitude':
        V = abs(V)
        part = 'magnitude'
    elif plot_type == 'dB':
        Vabs = abs(V)
        Vmin = 10**(dbmin / 20)
        if np.all(Vabs < Vmin):
            warn('All values below dbmin')
        Vabs[Vabs < Vmin] = Vmin
        V = 20 * np.log10(Vabs)
        part = 'magnitude'
        units = 'dB'
    elif plot_type == 'radians':
        V = np.angle(V)
        if unwrap:
            V = np.unwrap(V)
        part = 'phase'
        units = 'radians'
    elif plot_type == 'degrees':
        V = np.angle(V)
        if unwrap:
            V = np.unwrap(V)
        V = np.degrees(V)
        part = 'phase'
        units = 'degrees'

    if obj.is_magnitude or np.all(V > 0):
        plot = plots[(log_magnitude, log_frequency)]
    else:
        plot = plots[(False, log_frequency)]

    if norm:
        default_xlabel = 'Normalized frequency'
    else:
        default_xlabel = obj.domain_label_with_units

    if label == '':
        default_ylabel = part.capitalize()
    else:
        default_ylabel = label + ' ' + part
    if units != '' and units != '1':
        default_ylabel += ' (' + units + ')'

    xlabel = kwargs.pop('xlabel', default_xlabel)
    ylabel = kwargs.pop('ylabel', default_ylabel)
    ylabel2 = kwargs.pop('ylabel2', default_ylabel)
    second = kwargs.pop('second', False)
    xscale = kwargs.pop('xscale', 1)
    yscale = kwargs.pop('yscale', 1)
    title = kwargs.pop('title', None)

    lines = plot(f * xscale, V * yscale, **kwargs)
    color = kwargs.pop('color', lines[0].get_color())

    if deltas is not None:
        plot_deltas(ax, f * xscale, deltas, obj.var, plot_type, color)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    ylabel = ylabel2 if second else ylabel
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

    ax.grid(True)
    return lines


def plot_frequency(obj, f, plot_type=None, **kwargs):
    """This is a helper function for a frequency response plot.  It is
    better to use the `plot()` method of a FourierDomainExpression."""

    obj = obj.doit()

    # Much of the hoop jumping is to speed up plotting since
    # obj.real can be slow.  Instead we evaluate complex
    # objects and then convert to phase, magnitude, etc.

    norm = kwargs.pop('norm', False)
    npoints = kwargs.pop('npoints', 400)
    log_magnitude = kwargs.pop('log_magnitude', False)
    log_frequency = kwargs.pop(
        'log_frequency', False) or kwargs.pop('log_scale', False)
    if kwargs.pop('loglog', False):
        log_magnitude = True
        log_frequency = True

    dbmin = kwargs.pop('dbmin', -120)
    label = kwargs.pop('label', None)
    nyticks = kwargs.pop('nyticks', None)
    color2 = kwargs.pop('color2', None)
    linestyle2 = kwargs.pop('linestyle2', '--')
    unwrap = kwargs.pop('unwrap', False)

    # FIXME, determine useful frequency range...
    if f is None:
        if norm:
            f = (-0.5, 0.5)
        else:
            f = (0, 2)
    if isinstance(f, (int, float)):
        f = (0, f)
    if isinstance(f, tuple):
        if log_frequency:
            fmin, fmax = parse_range(f, 1e-1, positive=True)
            f = np.geomspace(fmin, fmax, npoints)
        else:
            fmin, fmax = parse_range(f, 1e-1, positive=False)
            f = np.linspace(fmin, fmax, npoints)

    # Objects can have a `part` attribute that is set by methods such
    # as real, imag, phase, magnitude.  If this is defined,
    # `plot_type` is ignored.

    plot1_type = 'default'
    plot2_type = None

    V = obj.evaluate(f)

    types = ['dB-phase', 'dB-radians', 'dB-phase-degrees',
             'dB-degrees', 'mag-phase', 'magnitude-phase',
             'mag-phase-degrees', 'magnitude-phase-degrees',
             'real-imag', 'mag', 'magnitude', 'phase', 'radians',
             'phase-degrees', 'degrees', 'real', 'imag', 'dB', 'abs']
    if plot_type is not None and plot_type not in types:
        raise ValueError('Unknown plot type %s, expecting: %s ' %
                         (plot_type, ', '.join(types)))

    if not V.dtype == complex:
        if plot_type is None:
            plot_type = 'real'

    elif obj.part == '':
        if plot_type is None:
            plot_type = 'dB-phase'

    if plot_type in ('dB_phase', 'dB-phase', 'dB-radians'):
        plot1_type = 'dB'
        plot2_type = 'radians'
    elif plot_type in ('dB_phase_degrees', 'dB-phase-degrees',
                       'dB-degrees'):
        plot1_type = 'dB'
        plot2_type = 'degrees'
    elif plot_type in ('mag_phase', 'magnitude_phase',
                       'mag-phase', 'magnitude-phase'):
        plot1_type = 'magnitude'
        plot2_type = 'radians'
    elif plot_type in ('mag_phase_degrees', 'magnitude_phase_degrees',
                       'mag-phase-degrees', 'magnitude-phase-degrees'):
        plot1_type = 'magnitude'
        plot2_type = 'degrees'
    elif plot_type in ('real_imag', 'real-imag'):
        plot1_type = 'real'
        plot2_type = 'imag'
    elif plot_type in ('mag', 'magnitude'):
        plot1_type = 'magnitude'
    elif plot_type in ('phase', 'radians'):
        plot1_type = 'radians'
    elif plot_type in ('phase-degrees', 'degrees'):
        plot1_type = 'degrees'
    elif plot_type == 'abs':
        plot1_type = 'abs'
    elif plot_type == 'real':
        plot1_type = 'real'
    elif plot_type == 'imag':
        plot1_type = 'imag'
    elif plot_type == 'dB':
        plot1_type = 'dB'
    else:
        raise RuntimeError('Internal error')

    deltas = None
    if kwargs.pop('plot_deltas', True) and obj.has(DiracDelta):
        cls = obj.__class__
        rest, deltas = separate_dirac_delta(obj.sympy.expand())
        obj = cls(rest, **obj.assumptions)
        V = obj.evaluate(f)

    ax = make_axes(figsize=kwargs.pop('figsize', None),
                   axes=kwargs.pop('axes', None))

    if label is None:
        label = plot1_type

    lines = plotit(ax, obj, f, V, plot1_type, deltas,
                   log_frequency=log_frequency,
                   log_magnitude=log_magnitude, norm=norm,
                   label=label, dbmin=dbmin, unwrap=unwrap, **kwargs)

    if plot2_type is None:
        return ax

    if color2 is None:
        color2 = lines[0].get_color()

    kwargs2 = kwargs.copy()
    for key in ('color', 'linestyle', 'xlabel', 'xlabel2', 'ylabel', 'title', 'xscale', 'yscale'):
        kwargs2.pop(key, None)

    # Dummy plot to add label to legend.
    ax.plot([], [], label=plot2_type, color=color2,
            linestyle=linestyle2, **kwargs2)

    ax2 = ax.twinx()
    kwargs['axes'] = ax2
    kwargs.pop('color', None)

    plotit(ax2, obj, f, V, plot2_type, deltas, log_frequency=log_frequency,
           log_magnitude=log_magnitude, norm=norm, second=True, color=color2,
           linestyle=linestyle2, unwrap=unwrap, **kwargs)

    if plot2_type in ('phase', 'radians') and not unwrap:
        ax2.set_ylim(-np.pi, np.pi)
    elif plot2_type in ('phase-degrees', 'degrees') and not unwrap:
        ax2.set_ylim(-180, 180)

    if plot1_type == 'dB':
        Vabs = abs(V)
        Vmin = 10**(dbmin / 20)
        if np.all(Vabs < Vmin):
            warn('All values below dbmin')
        Vabs[Vabs < Vmin] = Vmin

        m = np.isinf(Vabs)
        if any(m):
            warn('Approximating infinite value!')
            Vabs[m] = 1e30

        dB = 20 * np.log10(Vabs)
        dBmin = min(dB)
        dBmax = max(dB)

        ymin = np.floor(dBmin / 10) * 10
        ymax = np.ceil(dBmax / 10) * 10
        yrange = ymax - ymin

        if False:
            # Provide some wiggle room
            ymin -= yrange / 10
            ymax += yrange / 10
            yrange = ymax - ymin

        # This will clobber previous plot values
        ax.set_ylim(ymin, ymax)

        if nyticks is None:
            # An odd number is better to show zero phase.
            for nyticks in (7, 6, 5, 8, 9):
                if yrange / (nyticks - 1) in (20, 10, 5, 2):
                    break

    from matplotlib import ticker

    if nyticks is None:
        nyticks = 5

    ax.yaxis.set_major_locator(ticker.LinearLocator(nyticks))
    ax2.yaxis.set_major_locator(ticker.LinearLocator(nyticks))

    ax.legend()
    return ax, ax2


def plot_bode(obj, f, **kwargs):
    """This is a helper function for a Bode plot.  It is
    better to use the `bode_plot()` method of a FourierDomainExpression."""

    if 'log_frequency' not in kwargs:
        kwargs['log_frequency'] = True
    if 'unwrap' not in kwargs:
        kwargs['unwrap'] = True

    return plot_frequency(obj, f, log_magnitude=True, **kwargs)


def plot_angular_bode(obj, f, **kwargs):
    """This is a helper function for a Bode plot.  It is better to use the
    `bode_plot()` method of an AngularFourierDomainExpression."""

    if 'log_frequency' not in kwargs:
        kwargs['log_frequency'] = True

    return plot_angular_frequency(obj, f, **kwargs)


def plot_angular_frequency(obj, omega, plot_type=None, **kwargs):
    """This is a helper function for a frequency response plot.  It is
    better to use the `plot()` method of an AngularFourierDomainExpression."""

    npoints = kwargs.pop('npoints', 400)
    norm = kwargs.pop('norm', False)
    log_frequency = kwargs.pop(
        'log_frequency', False) or kwargs.pop('log_scale', False)

    # FIXME, determine useful frequency range...
    if omega is None:
        if norm:
            omega = (-np.pi, np.pi)
        else:
            omega = (0, np.pi)
    if isinstance(omega, (int, float)):
        omega = (0, omega)
    if isinstance(omega, tuple):
        if log_frequency:
            wmin, wmax = parse_range(omega, 1e-1, positive=True)
            omega = np.geomspace(wmin, wmax, npoints)
        else:
            wmin, wmax = parse_range(omega, 1e-1, positive=False)
            omega = np.linspace(wmin, wmax, npoints)

    if norm and 'xlabel' not in kwargs:
        kwargs['xlabel'] = 'Normalized angular frequency'

    return plot_frequency(obj, omega, plot_type=plot_type,
                          norm=norm, log_frequency=log_frequency, **kwargs)


def plot_time(obj, t, plot_type=None, **kwargs):
    """This is a helper function for a time plot.  It is
    better to use the `plot()` method of a TimeDomainExpression."""

    obj = obj.doit()
    npoints = kwargs.pop('npoints', 400)

    # FIXME, determine useful time range...
    if t is None:
        t = (-0.2, 2)
    if isinstance(t, (int, float)):
        t = (0, t)
    if isinstance(t, tuple):
        tmin, tmax = parse_range(t, 1e-1, positive=False)
        t = np.linspace(tmin, tmax, npoints)

    deltas = None
    if kwargs.pop('plot_deltas', True) and obj.has(DiracDelta):
        cls = obj.__class__
        rest, deltas = separate_dirac_delta(obj.expr)
        obj = cls(rest, **obj.assumptions)

    v = obj.evaluate(t)

    if plot_type == None:
        plot_type = 'real'
    if plot_type == 'real':
        v = v.real
    elif plot_type == 'imag':
        v = v.imag
    elif plot_type == 'abs':
        v = abs(v)
    else:
        raise ValueError(
            'Invalid plot_type: must be real, imag, real-imag, abs')

    ax = make_axes(figsize=kwargs.pop('figsize', None),
                   axes=kwargs.pop('axes', None))

    xlabel = kwargs.pop('xlabel', obj.domain_label_with_units)
    ylabel = kwargs.pop('ylabel', obj.label_with_units)
    xscale = kwargs.pop('xscale', 1)
    yscale = kwargs.pop('yscale', 1)
    title = kwargs.pop('title', None)

    lines = ax.plot(t * xscale, v * yscale, **kwargs)
    color = kwargs.pop('color', lines[0].get_color())

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

    if deltas is not None:
        # TODO, fix yscale
        plot_deltas(ax, t * xscale, deltas, obj.var, plot_type, color=color)

    ax.grid(True)
    return ax


# make a stem plot for complex values sequences in the complex domain
def plot_sequence_polar(obj, ni=(-10, 10), **kwargs):

    obj = obj.doit()
    npoints = kwargs.pop('npoints', 400)

    # FIXME, determine useful range...
    if ni is None:
        ni = (-20, 20)
    if isinstance(ni, tuple):
        # Use float data type since NumPy barfs for n**(-a) where a is
        # an integer.
        ni = np.arange(ni[0], ni[1] + 1, dtype=float)

    v = obj.evaluate(ni)

    ax = make_axes(figsize=kwargs.pop('figsize', None),
                   axes=kwargs.pop('axes', None),
                   subplot_kw=dict(polar=True))

    phi = np.angle(v)
    mag = abs(v)

    # Plot symbols
    lines = ax.plot(phi, mag, 'o', **kwargs)

    color = kwargs.pop('color', lines[0].get_color())

    # Plot lines from origin
    Nv = len(v)
    ax.plot((phi, phi), (np.zeros(Nv), mag), color=color, **kwargs)
    ax.grid(True)
    return ax


def plot_sequence(obj, ni, plot_type=None, polar=False, **kwargs):
    """This is a helper function for a stem (lollipop) plot.  It is
    better to use the `plot()` method of a DiscreteTimeDomainExpression."""

    from matplotlib.ticker import MaxNLocator

    obj = obj.doit()

    if polar:
        return plot_sequence_polar(obj, ni, **kwargs)

    npoints = kwargs.pop('npoints', 400)

    # FIXME, determine useful range...
    if ni is None:
        ni = (-20, 20)
    if isinstance(ni, tuple):
        # Use float data type since NumPy barfs for n**(-a) where a is
        # an integer.
        ni = np.arange(ni[0], ni[1] + 1, dtype=float)

    v = obj.evaluate(ni)

    if plot_type == None:
        plot_type = 'real'
    if plot_type == 'real':
        v = v.real
    elif plot_type == 'imag':
        v = v.imag
    elif plot_type == 'abs':
        v = abs(v)
    else:
        raise ValueError(
            'Invalid plot_type: must be real, imag, real-imag, abs')

    deltas = None
    if obj.has(DiracDelta):
        cls = obj.__class__
        rest, deltas = separate_dirac_delta(obj.expr)
        obj = cls(rest, **obj.assumptions)

    ax = make_axes(figsize=kwargs.pop('figsize', None),
                   axes=kwargs.pop('axes', None))

    xlabel = kwargs.pop('xlabel', obj.domain_label_with_units)
    ylabel = kwargs.pop('ylabel', obj.label_with_units)
    xscale = kwargs.pop('xscale', 1)
    yscale = kwargs.pop('yscale', 1)
    title = kwargs.pop('title', None)
    color = kwargs.pop('color', None)

    markerline, stemlines, baseline = ax.stem(
        ni * xscale, v * yscale, use_line_collection=True, **kwargs)
    if color is not None:
        stemlines.set_color(color)
        markerline.set_color(color)
    else:
        color = stemlines.get_color()

    # Ensure integer ticks.
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    if deltas is not None:
        # TODO, fix yscale
        plot_deltas(ax, ni * xscale, deltas, obj.var, plot_type, color=color)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

    ax.grid(True)
    return ax


def plot_phasor(obj, **kwargs):
    """This is a helper function for a phasor plot.  It is
    better to use the `plot()` method of a PhasorDomainExpression."""

    obj = obj.doit()

    ax = make_axes(figsize=kwargs.pop('figsize', None),
                   axes=kwargs.pop('axes', None),
                   subplot_kw=dict(polar=True))

    phi = obj.phase.fval
    mag = obj.magnitude.fval

    ax.plot((phi, phi), (0, mag), **kwargs)
    return ax


def plot_nyquist(obj, f, norm=False, **kwargs):
    """This is a helper function for a Nyquist plot.  It is better to use
    the `nyquist_plot()` method of a LaplaceDomainExpression."""

    from matplotlib.pyplot import Circle, rcParams

    npoints = kwargs.pop('npoints', 400)
    # Default to True to get better plots.
    log_frequency = kwargs.pop('log_frequency', False)
    unitcircle = kwargs.pop('unitcircle', True)

    fn = None

    # FIXME, determine useful frequency range...
    if f is None:
        f = (-100, 100)
    if isinstance(f, (int, float)):
        f = (0, f)
    if isinstance(f, tuple):
        fmin, fmax = parse_range(f, 1e-1, positive=False)
        if log_frequency:
            if fmin < 0:
                frange = fmax - fmin
                f = np.geomspace(frange / 1e4, fmax, npoints // 2)
                f = np.hstack((0, f))
                fn = -np.geomspace(frange / 1e4, abs(fmin), npoints // 2)
                fn = np.hstack((0, fn))
            elif fmin == 0:
                f = np.geomspace(fmax / 1e4, fmax, npoints // 2)
                f = np.hstack((0, f))
            else:
                f = np.geomspace(fmin, fmax, npoints)
        else:
            fmin, fmax = parse_range(f, 1e-1, positive=False)
            f = np.linspace(fmin, fmax, npoints)

    V = obj.evaluate(f)
    if not V.dtype == complex:
        raise ValueError('Data not complex')

    obj = obj.doit()

    ax = make_axes(figsize=kwargs.pop('figsize', None),
                   axes=kwargs.pop('axes', None))

    if unitcircle:
        ax.add_artist(Circle((0, 0), 1, color='blue',
                      linestyle='--', fill=False))

    lines = ax.plot(V.real, V.imag)

    if fn is not None:
        V = obj.evaluate(fn)
        color = lines[0].get_color()
        ax.plot(V.real, V.imag, color=color)

    x_min = V.real.min()
    x_max = V.real.max()
    y_min = V.imag.min()
    y_max = V.imag.max()

    if unitcircle:
        if x_min > -1:
            x_min = -1
        if x_max < 1:
            x_max = 1
        if y_min > -1:
            y_min = -1
        if y_max < 1:
            y_max = 1

    x_extra, y_extra = 1.0, 1.0
    x_min -= 0.5 * x_extra
    x_max += 0.5 * x_extra
    if unitcircle:
        bbox = ax.get_window_extent()
        aspect = bbox.width / bbox.height

        x_min *= aspect
        x_max *= aspect

    ax.axis('equal')
    #ax.set_xlim(x_min, x_max)

    xlabel = kwargs.pop('xlabel', 'Re')
    ylabel = kwargs.pop('ylabel', 'Im')
    title = kwargs.pop('title', None)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

    ax.grid(True)

    return ax


def plot_nichols(obj, f, norm=False, **kwargs):
    """This is a helper function for a Nichols plot.  It is better to use
    the `nichols_plot()` method of a LaplaceDomainExpression."""

    from matplotlib.pyplot import Circle, rcParams

    npoints = kwargs.pop('npoints', 400)
    # Default to True to get better plots.
    log_frequency = kwargs.pop('log_frequency', False)

    fn = None

    # FIXME, determine useful frequency range...
    if f is None:
        f = (-100, 100)
    if isinstance(f, (int, float)):
        f = (0, f)
    if isinstance(f, tuple):
        fmin, fmax = parse_range(f, 1e-1, positive=False)
        if log_frequency:
            if fmin < 0:
                frange = fmax - fmin
                f = np.geomspace(frange / 1e4, fmax, npoints // 2)
                f = np.hstack((0, f))
                fn = -np.geomspace(frange / 1e4, abs(fmin), npoints // 2)
                fn = np.hstack((0, fn))
            elif fmin == 0:
                f = np.geomspace(fmax / 1e4, fmax, npoints // 2)
                f = np.hstack((0, f))
            else:
                f = np.geomspace(fmin, fmax, npoints)
        else:
            fmin, fmax = parse_range(f, 1e-1, positive=False)
            f = np.linspace(fmin, fmax, npoints)

    V = obj.evaluate(f)
    if not V.dtype == complex:
        raise ValueError('Data not complex')

    obj = obj.doit()

    ax = make_axes(figsize=kwargs.pop('figsize', None),
                   axes=kwargs.pop('axes', None))

    lines = ax.plot(np.angle(V), 20 * np.log10(abs(V)))

    if fn is not None:
        V = obj.evaluate(fn)
        color = lines[0].get_color()
        ax.plot(np.angle(V), 20 * np.log10(abs(V)), color=color)

    xlabel = kwargs.pop('xlabel', 'Phase (radians)')
    ylabel = kwargs.pop('ylabel', 'dB')
    title = kwargs.pop('title', None)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

    ax.grid(True)

    return ax
