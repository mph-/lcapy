"""This module provides the Fitter class.  This uses optimization
techniques to find the parameters of an expression that best fits measured data.

Copyright 2022 Michael Hayes, UCECE

"""

from scipy.optimize import brute, fmin, curve_fit, minimize
from numpy import iscomplexobj, hstack, zeros, sqrt


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

    def _make_ranges(self):

        ranges = {}
        for symbol in self.symbols:
            # Perhaps make (-inf, inf) for unbounded
            # but will need to need to fix initial guess.
            ranges[symbol] = (-1e9, 1e9)
        return ranges

    def model(self, params, x, ranges):

        defs = self._make_defs(params, ranges)

        return self.expr.subs(defs).evaluate(x)

    def _optimize_brute(self, x, y, ranges=None, Ns=10, finish='fmin', **kwargs):

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

    def _optimize_minimize(self, x, y, ranges=None, method='Nelder-Mead', **kwargs):

        kwargs.pop('Ns', None)
        kwargs.pop('finish', None)

        bounds = []
        p0 = []

        for m, r in enumerate(ranges.values()):
            if len(r) in (2, 3):
                bounds.append((r[0], r[1]))
                p0.append(0.5 * (r[0] + r[1]))
            else:
                raise ValueError('Range %s can only have 2 or 3 values' % r)

        iscomplex = iscomplexobj(y)

        if iscomplex:
            y = hstack((y.real, y.imag))

        def func(params):

            yp = self.model(params, x, ranges)

            if iscomplex:
                return (abs(y - yp)**2).mean()
            return ((y - yp)**2).mean()

        results = minimize(func, x0=p0, bounds=bounds, method=method, **kwargs)

        params = results.x

        yp = func(params)
        rmse = ((y - yp)**2).mean()

        defs = self._make_defs(params, ranges)
        return FitterResult(defs, rmse)

    def _optimize1(self, x, y, ranges=None, method='trf', **kwargs):

        if method == 'brute':
            return self._optimize_brute(x, y, ranges, **kwargs)
        elif method in ('lm', 'trf', 'dogbox'):
            return self._optimize_curvefit(x, y, ranges, method, **kwargs)
        elif method in ('Nelder-Mead', 'Powell'):
            return self._optimize_minimize(x, y, ranges, method, **kwargs)
        else:
            return self._optimize_minimize(x, y, ranges, method, **kwargs)

    def optimize(self, x, y, ranges=None, method='trf', iterations=1, **kwargs):

        if ranges is None:
            ranges = self._make_ranges()

        iscomplex = iscomplexobj(y)

        if iterations <= 1:
            return self._optimize1(x, y, ranges, method, **kwargs)

        tol = kwargs.get('tol', 1e-6)

        rmse_prev = None
        x1 = x
        y1 = y
        for i in range(iterations):

            result = self._optimize1(x1, y1, ranges, method, **kwargs)

            params = result.params
            yp = self.expr.subs(params).evaluate(x)

            e = abs(y - yp)
            sd = e.std()
            m = e <= 3 * sd
            num_outliers = len(x) - sum(m)
            if num_outliers == 0:
                break

            x1 = x[m]
            y1 = y[m]

            if len(x1) == 0:
                raise ValueError('Outlier removal failed')

            e = abs(y1 - yp[m])
            rmse = sqrt((e**2).mean())

            if False:
                print('Removed %d outliers at iteration %d: sd=%f, rmse=%f' %
                      (num_outliers, i + 1, sd, rmse))

            result.rmse = rmse
            if rmse < tol:
                break
            if rmse_prev is not None and rmse == rmse_prev:
                break
            rmse_prev = rmse

        return result


def fit(expr, x, y, ranges, method='trf', Ns=10, **kwargs):

    return Fitter(expr).optimize(x, y, ranges=ranges, method=method, Ns=Ns, **kwargs)


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
