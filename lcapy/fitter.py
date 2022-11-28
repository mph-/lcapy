"""This module provides the Fitter class.  This uses optimization
techniques to find the parameters of an expression that best fits measured data.

Copyright 2022 Michael Hayes, UCECE

"""

from scipy.optimize import brute, fmin, curve_fit
from numpy import iscomplexobj, hstack, zeros


class FitterResult(object):

    def __init__(self, params, rmse):

        self.params = params
        self.rmse = rmse


class Fitter(object):
    """Data fitting."""

    def __init__(self, expr):

        self.expr = expr
        self.symbols = expr.symbols
        self.symbols.pop(str(expr.var))

    def _make_defs(self, params, ranges):

        defs = {}
        for m, r in enumerate(ranges.items()):
            defs[r[0]] = params[m]
        return defs

    def model(self, params, x, ranges):

        defs = self._make_defs(params, ranges)

        return self.expr.subs(defs).evaluate(x)

    def _optimize_brute(self, x, y, ranges=None, Ns=10, finish='fmin', **kwargs):
        """Ranges is a list of tuples, of the form: (min, max) or (min, max,
        numsteps).  If `numsteps` is not specified then `Ns` is used."""

        kwargs.pop('method', None)
        self.verbose = kwargs.pop('verbose', 0)

        if finish in ('none', 'None', ''):
            finish = None
        elif finish == 'fmin':
            finish = fmin

        oranges = []
        for k, r in ranges.items():
            if len(r) == 2:
                # Note, a complex value specifies the number of steps
                # (see numpy.mgrid).
                oranges.append(slice(r[0], r[1], complex(Ns)))
            elif len(r) == 3:
                oranges.append(slice(r[0], r[1], complex(r[2])))
            else:
                raise ValueError('Range %s can only have 2 or 3 values' % r)

        iscomplex = iscomplexobj(y)

        def func(params):

            yp = self.model(params, x, ranges)
            if iscomplex:
                return (abs(y - yp)**2).mean()

            return ((y - yp)**2).mean()

        params, rmse, foo, bar = brute(func, oranges, Ns=Ns, args=(),
                                       finish=finish, full_output=1)

        defs = self._make_defs(params, ranges)
        return FitterResult(defs, rmse)

    def _optimize_curvefit(self, x, y, ranges=None, method='trf', ftol=1e-14, xtol=1e-14,
                           maxfev=1e5, **kwargs):
        """Ranges is a list of tuples, of the form: (min, max)."""

        kwargs.pop('Ns', None)
        kwargs.pop('finish', None)

        bounds_min = zeros(len(ranges))
        bounds_max = zeros(len(ranges))

        for m, r in enumerate(ranges.values()):
            if len(r) in (2, 3):
                bounds_min[m] = r[0]
                bounds_max[m] = r[1]
            else:
                raise ValueError('Range %s can only have 2 or 3 values' % r)

        bounds = (bounds_min, bounds_max)

        # Initial guess.
        p0 = 0.5 * (bounds_min + bounds_max)

        iscomplex = iscomplexobj(y)

        if iscomplex:
            y = hstack((y.real, y.imag))

        def func(x, *params):

            yp = self.model(params, x, ranges)
            if iscomplex:
                yp = hstack((yp.real, yp.imag))
            return yp

        params, cov = curve_fit(func, x, y, p0=p0,
                                bounds=bounds, ftol=ftol, xtol=xtol,
                                maxfev=maxfev, **kwargs)

        yp = func(x, *params)
        rmse = ((y - yp)**2).mean()

        defs = self._make_defs(params, ranges)
        return FitterResult(defs, rmse)

    def optimize(self, x, y, ranges=None, method='trf', **kwargs):

        if method == 'brute':
            return self._optimize_brute(x, y, ranges, **kwargs)
        elif method in ('trf', 'dogbox'):
            return self._optimize_curvefit(x, y, ranges, method, **kwargs)
        else:
            raise ValueError(
                'Unknown method %s: needs to be brute, trf, dogbox' % method)


def fit(expr, x, y, method='trf', ranges=None, Ns=10, **kwargs):

    return Fitter(expr).optimize(x, y, method=method, ranges=ranges, Ns=Ns, **kwargs)


def test1():

    from lcapy import expr, f
    from numpy import arange

    e = expr('a * exp(-t  / tau) * u(t)')

    tv = arange(100)

    vv = e.subs({'a': 1, 'tau': 10}).evaluate(tv)

    results = fit(e, tv, vv, method='brute',
                  ranges={'a': (0, 10), 'tau': (1, 20)})
    print(results.params)

    results = fit(e, tv, vv, method='trf',
                  ranges={'a': (0, 10), 'tau': (1, 20)})
    print(results.params)

    E = e(f)

    fv = arange(100)

    Vv = E.subs({'a': 1, 'tau': 10}).evaluate(fv)

    results = fit(E, fv, Vv, method='brute',
                  ranges={'a': (0, 10), 'tau': (1, 20)})
    print(results.params)

    results = fit(E, fv, Vv, method='trf',
                  ranges={'a': (0, 10), 'tau': (1, 20)})
    print(results.params)
