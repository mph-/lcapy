import sympy as sym


def laplace_limits(expr, t, s, tmin, tmax):
    
    F = sym.integrate(expr * sym.exp(-s * t), (t, tmin, tmax))

    if not F.has(sym.Integral):
        return F

    if not F.is_Piecewise:
        raise ValueError('Could not compute Laplace transform for ' + str(expr))

    F, cond = F.args[0]
    if F.has(sym.Integral):
        raise ValueError('Could not compute Laplace transform for ' + str(expr))

    return F

def laplace_0minus(expr, t, s):
    
    t0 = sym.symbols('t0', negative=True, real=True)

    F = laplace_limits(expr, t, s, t0, sym.oo)
    return sym.limit(F, t0, 0)

def laplace_0(expr, t, s):

    return laplace_limits(expr, t, s, 0, sym.oo)

def laplace_term(expr, t, s):

    var = sym.Symbol(str(t))
    expr = expr.replace(var, t)

    if expr.has(sym.function.AppliedUndef):
        # TODO, handle things like 3 * v(t), a * v(t), 3 * t * v(t)
        if not isinstance(expr, sym.function.AppliedUndef):
            raise ValueError('Could not compute Laplace transform for ' + str(expr))
        # Convert v(t) to V(s), etc.
        if expr.func.__name__ == 'u':
            print('Warning, use Heaviside(t) for unit step u(t)')

        name = expr.func.__name__
        name = name[0].upper() + name[1:] + '(s)'
        return sym.sympify(name)

    if expr.has(sym.Heaviside(t)):
        return laplace_0(expr.replace(sym.Heaviside(t), 1), t, s)

    if expr.has(sym.DiracDelta) or expr.has(sym.Heaviside):
        try:
            return laplace_0minus(expr, t, s)
        except ValueError:
            pass

    return laplace_0(expr, t, s)

def laplace_transform(expr, t, s):
    """Compute single sided Laplace transform of expr
    with lower limit 0-.

    Undefined functions such as v(t) are converted to V(s)

    """

    # The variable may have been created with different attributes,
    # say when using sym.sympify('Heaviside(t)') since this will
    # default to assuming that t is complex.  So if the symbol has the
    # same representation, convert to the desired one.

    var = sym.Symbol(str(t))
    expr = expr.replace(var, t)

    terms = expr.as_ordered_terms()
    result = 0

    try:
        for term in terms:
            result += laplace_term(term, t, s)
    except ValueError:
        raise ValueError('Could not compute Laplace transform for ' + str(expr))

    return result

