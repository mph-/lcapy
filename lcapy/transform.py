"""This module performs transformations between domains.

Copyright 2018--2022 Michael Hayes, UCECE

"""

from .impedance import impedance  # nope8p
from .sym import sympify, pi
from .sym import dt
from .symbols import f, s, t, omega, j, jf, jw, jw0, Omega, F, k, n, z
from .symbols import domain_vars, domain_var_ids
from .nexpr import n
from .expr import expr as expr1
from .expr import Expr


def transform(expr, arg, **assumptions):
    """If arg is a domain variable perform domain transformation,
    otherwise perform substitution.

    Note (1 / s)(omega) will fail since 1 / s is assumed not to be
    causal and so the Fourier transform is unknown.  However,
    impedance(1 / s)(omega) will work since an impedance is assumed to
    be causal.  Alternatively, use (1 / s)(omega, causal=True).

    Transforming from s->jomega is fast since it just requires
    a substitution of s with jomega.

    Transforming from s->omega, s->Omega, s->f, or s->F can be slow since
    this requires a inverse Laplace transform followed by a Fourier
    transform.  However, if the expression is causal and the
    expression is lossy when s is replaced by jw, the result can be
    found by substituting jw or 2 * 2 * pi * f for s.  This does not
    apply for an expression such as Z = 1 / (s * C).

    """

    arg = expr1(arg)

    # Handle things like (3 * s)(5 * s)
    if isinstance(expr, arg.__class__) and not isinstance(expr, Superposition):
        return expr.subs(arg)

    # Handle expr(t), expr(s), expr(f), expr(omega), expr(jomega), expr(3j)
    try:
        if arg is t:
            return expr.time(**assumptions)
        elif arg is s:
            return expr.laplace(**assumptions)
        elif arg is f:
            return expr.fourier(**assumptions)
        elif arg is omega:
            return expr.angular_fourier(**assumptions)
        elif arg is Omega:
            return expr.norm_angular_fourier(**assumptions)
        elif arg is F:
            return expr.norm_fourier(**assumptions)
        elif arg is jf:
            return expr.frequency_response(**assumptions)
        elif arg is jw:
            return expr.angular_frequency_response(**assumptions)
        elif arg.has(j):
            # Note, substituting s = jw is valid for the phasor domain
            # but does not always produce the freqguency response since
            # the evaluation point can be outside the region of
            # convergence.  For example, the LT for u(t), t * u(t),
            # and exp(a * t) for positive a, i.e., the impulse
            # responses for marginally stable and unstable systems.
            # For the frequency response use omega.
            return expr.phasor_ratio(omega=arg / j, **assumptions)
        elif arg is n:
            return expr.discrete_time(**assumptions)
        elif arg is k:
            return expr.discrete_frequency(**assumptions)
        elif arg is z:
            return expr.zdomain(**assumptions)
    except AttributeError:
        raise AttributeError('Transform of %s domain to %s domain not implemented' % (
            expr.domain, arg.domain))

    # Handle expr(texpr), expr(sexpr), expr(fexpr), expr(omegaexpr).
    # For example, expr(2 * f).
    result = None
    if isinstance(arg, TimeDomainExpression):
        result = expr.time(**assumptions)
    elif isinstance(arg, LaplaceDomainExpression):
        result = expr.laplace(**assumptions)
    elif isinstance(arg, FourierDomainExpression):
        result = expr.fourier(**assumptions)
    elif isinstance(arg, NormFourierDomainExpression):
        result = expr.norm_fourier(**assumptions)
    elif isinstance(arg, AngularFourierDomainExpression):
        result = expr.angular_fourier(**assumptions)
    elif isinstance(arg, NormAngularFourierDomainExpression):
        result = expr.norm_angular_fourier(**assumptions)
    elif isinstance(arg, AngularFrequencyResponseDomainExpression):
        result = expr.angular_frequency(**assumptions)
    elif arg.has(j):
        result = expr.phasor(omega=arg / j, **assumptions)
    elif isinstance(arg, DiscreteTimeDomainExpression):
        result = expr.discrete_time(**assumptions)
    elif isinstance(arg, DiscreteFourierDomainExpression):
        result = expr.discrete_frequency(**assumptions)
    elif isinstance(arg, ZDomainExpression):
        result = expr.zdomain(**assumptions)
    elif arg.is_constant:
        if not isinstance(expr, Superposition):
            result = expr.time(**assumptions)
        else:
            result = expr
    else:
        raise ValueError(
            'Can only return t, f, s, F, omega, Omega, or jw domains')

    return result.subs(arg, **assumptions)


def call(expr, arg, **assumptions):

    from numpy import ndarray, array

    if isinstance(arg, (tuple, list)):
        return [expr._subs1(expr.var, arg1) for arg1 in arg]

    elif isinstance(arg, ndarray):
        return array([expr._subs1(expr.var, arg1) for arg1 in arg])

    arg = expr1(arg)
    if arg.is_constant:
        return expr.subs(arg)

    return expr.transform(arg, **assumptions)


def select(expr, kind):

    if not isinstance(kind, str):
        return expr.subs(j * kind)

    # If kind is an expr, then will add 'dc', 'time', etc. as symbols.

    if kind == 't':
        return expr.time()
    elif kind in ('dc', 'time'):
        return expr.subs(0)
    elif kind in ('s', 'ivp', 'super', 'laplace'):
        return expr.laplace()
    elif kind == 'f':
        return expr.fourier()
    elif kind == 'omega':
        return expr.angular_fourier()
    elif kind == 'Omega':
        return expr.norm_angular_fourier()
    elif kind == 'F':
        return expr.norm_fourier()
    elif isinstance(kind, str) and kind.startswith('n'):
        return expr.angular_fourier()
    else:
        raise RuntimeError('unknown kind')


from .fexpr import FourierDomainExpression  # nopep8
from .sexpr import LaplaceDomainExpression  # nopep8
from .texpr import TimeDomainExpression  # nopep8
from .omegaexpr import AngularFourierDomainExpression  # nopep8
from .normfexpr import NormFourierDomainExpression  # nopep8
from .normomegaexpr import NormAngularFourierDomainExpression  # nopep8
from .jomegaexpr import AngularFrequencyResponseDomainExpression  # nopep8
from .nexpr import DiscreteTimeDomainExpression  # nopep8
from .kexpr import DiscreteFourierDomainExpression  # nopep8
from .zexpr import ZDomainExpression  # nopep8
from .super import Superposition  # nopep8
from .current import current  # nopep8
from .voltage import voltage  # nopep8
from .admittance import admittance  # nopep8
from .transfer import transfer  # nopep8
