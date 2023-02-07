"""This module provides common utility functions.

Copyright 2021--2022 Michael Hayes, UCECE
"""

from warnings import warn
import sympy as sym


def factor_const(expr, var):
    """Extract constant factor from expression and return tuple
    of constant and the rest of the expression.

    For example `a * r(var)` returns `a, r(var)`.

    If have a polynomial expression, the leading coefficient is
    returned as the constant, for example: `2 * s + 4` returns `2, s + 2`.

    """

    # Perhaps use expr.as_coeff_Mul() ?

    if expr.is_polynomial():
        poly = sym.Poly(expr, var)
        const = poly.LC()
        if const == 0 or const == 1:
            return 1, expr
        return const, expr / const

    rest = sym.S.One
    const = sym.S.One
    for factor in expr.as_ordered_factors():
        # Cannot use factor.is_constant() since SymPy 1.2, 1.3
        # barfs for Heaviside(t) and DiracDelta(t)
        if not factor.has(var):
            const *= factor
        else:
            rest *= factor
    return const, rest


def term_const(expr, var):
    """Extract constant term from expression and return tuple
    of constant and the rest of the expression."""

    rest = sym.S.One
    const = sym.S.Zero
    for term in expr.as_ordered_terms():
        # Cannot use factor.is_constant() since SymPy 1.2, 1.3
        # barfs for Heaviside(t) and DiracDelta(t)
        if not term.has(var):
            const += term
        else:
            rest += term
    return const, rest


def scale_shift(expr, var):

    if expr == var:
        return sym.S.One, sym.S.Zero

    expr = expr.expand()
    if not expr.as_poly(var).is_linear:
        raise ValueError(
            'Expression not a linear function of %s: %s' % (var, expr))

    scale = expr.coeff(var, 1)
    shift = expr.coeff(var, 0)

    return scale, shift


def as_N_D(expr, var, monic_denominator=False, use_sympy=False):

    if use_sympy:
        return expr.as_numer_denom()

    factors = expr.as_ordered_factors()

    numers = []
    denoms = []
    for factor in factors:
        if factor.is_Function and factor.func == sym.exp:
            # SymPy treats exp(-s * a) as 1 / exp(s * a)
            numer = factor
            denom = sym.S.One
        else:
            numer, denom = factor.as_numer_denom()
        numers.append(numer)
        denoms.append(denom)

    poly_denom = False
    for denom in denoms:
        if denom != 1 and denom.is_polynomial(var):
            poly_denom = True
            break

    if not poly_denom:
        return expr.as_numer_denom()

    N = sym.S.One
    D = sym.S.One

    for numer, denom in zip(numers, denoms):
        N *= numer
        if denom.is_polynomial(var):
            D *= denom
        else:
            N /= denom

    N = N.simplify()

    if monic_denominator:
        Dpoly = sym.Poly(D, var)
        LC = Dpoly.LC()
        D = Dpoly.monic().as_expr()
        N = (N / LC).simplify()

    return N, D


def as_sum_terms(expr, var):

    N, D = as_N_D(expr, var)
    N = N.simplify()

    return [term / D for term in N.expand().as_ordered_terms()]


def as_sum(expr, var):

    result = 0
    for term in as_sum_terms(expr, var):
        result += term
    return result


def merge_common(lists):
    # From www.geeksforgeeks.org

    from collections import defaultdict

    neighbours = defaultdict(set)
    visited = set()
    for each in lists:
        for item in each:
            neighbours[item].update(each)

    def comp(node, neighbours=neighbours, visited=visited, visit=visited.add):

        nodes = set([node])
        next_node = nodes.pop
        while nodes:
            node = next_node()
            visit(node)
            nodes |= neighbours[node] - visited
            yield node

    for node in neighbours:
        if node not in visited:
            yield sorted(comp(node))


def isiterable(arg):

    return hasattr(arg, '__iter__')


def factor_expr(expr, factor):
    """Extract factor from expression or None if expression does
    not have factor."""

    factors = expr.as_ordered_factors()
    if factor not in factors:
        return None

    return expr / factor


def separate_dirac_delta(expr):
    """Separate Dirac delta terms from expression."""

    terms = expr.as_ordered_terms()
    deltas = []
    rest = 0

    for term in terms:
        if term.has(sym.DiracDelta):
            deltas.append(term)
        else:
            rest += term

    return rest, deltas


def split_dirac_delta(expr):
    """Return expression as a list of terms.
    The first term has no DiracDeltas, the second term collates
    the DiracDeltas, the third term collates derivatives of DiracDeltas, etc.

    For example, u(t) + DiractDelta(t, 1) returns [u(t), 0, DiracDelta(t, 1)]
    """

    terms = expr.as_ordered_terms()
    parts = {}
    rest = 0

    # FIXME, DiracDelta needs to be a factor

    for term in terms:
        if term.has(sym.DiracDelta):
            if len(term.args) == 1:
                if 1 not in parts:
                    parts[1] = 0
                parts[1] += term
            else:
                idx = term.args[1] + 1
                if idx not in parts:
                    parts[idx] = 0
                parts[idx] += term
        else:
            parts[0] = term

    maxkey = max(parts.keys())
    result = []
    for key in range(maxkey + 1):
        if key in parts:
            result.append(parts[key])
        else:
            result.append(0)

    return result


def remove_images(expr, var, dt, m1=0, m2=0):

    if m2 == 0 and isinstance(m1, tuple) and len(m1) == 2:
        # Perhaps should warn that this might be deprecated?
        m1, m2 = m1

    remove_all = m1 == 0 and m2 == 0

    const, expr1 = factor_const(expr, var)

    result = sym.S.Zero
    terms = expr1.as_ordered_terms()

    if len(terms) > 1:
        for term in expr1.as_ordered_terms():
            result += remove_images(term, var, dt, m1, m2)
        return const * result

    if not isinstance(expr1, sym.Sum):
        return expr

    sumsym = expr1.args[1].args[0]

    def query(expr):

        return expr.is_Add and expr.has(var) and expr.has(sumsym)

    def value(expr):
        if not expr.is_Add:
            return expr

        if not expr.is_polynomial(var) and not expr.as_poly(var).is_linear:
            return expr
        expr = expr.expand()
        a = expr.coeff(var, 1)
        b = expr.coeff(var, 0)

        if a == 0:
            return expr
        if b / a != -sumsym / dt:
            return expr
        return a * var

    expr1 = expr1.replace(query, value)

    if remove_all:
        return const * expr1.args[0]

    return const * sym.Sum(expr1.args[0], (sumsym, m1, m2))


def pair_conjugates(poles_dict):
    """Return dictionary of conjugate pole pairs and a dictionary of the
    remaining single poles."""

    pole_single_dict = poles_dict.copy()
    pole_pair_dict = {}

    pole_list = list(poles_dict)

    for i, pole in enumerate(pole_list):
        pole_c = sym.conjugate(pole)
        # Check for conjugate pole
        if pole_c in pole_list[i + 1:]:
            pole_single_dict.pop(pole, None)
            pole_single_dict.pop(pole_c, None)

            o1 = poles_dict[pole]
            o2 = poles_dict[pole_c]
            if o1 == o2:
                pole_pair_dict[pole, pole_c] = o1
            elif o1 > o2:
                pole_pair_dict[pole, pole_c] = o2
                pole_single_dict[pole] = o1 - o2
            else:
                pole_pair_dict[pole, pole_c] = o1
                pole_single_dict[pole_c] = o2 - o1

    return pole_pair_dict, pole_single_dict


def similarity_shift(expr, var):
    """Rewrite foo(a * t + b) as foo(t) and return a, b."""

    scale = None
    shift = None
    fail = False

    for expr1 in sym.preorder_traversal(expr):

        if not expr1.is_Function:
            continue

        arg = expr1.args[0]
        if not arg.has(var):
            continue

        poly = arg.as_poly(var)
        if poly is None or not poly.is_linear:
            fail = True
            break

        scale1 = arg.coeff(var, 1)
        shift1 = arg.coeff(var, 0)

        if scale is None:
            scale = scale1
        if shift is None:
            shift = shift1

        if scale != scale1 or shift != shift1:
            fail = True
            break

    if fail or scale is None:
        return expr, 1, 0

    expr2 = expr.replace(var * scale + shift, var)

    return expr2, scale, shift


def expand_functions(expr, var):

    # Try rewriting functions such as rampstep(t - 1)
    # to ramp(t - 1) - ramp(t - 2) to sums of weighted Heavisides...

    const, expr = factor_const(expr, var)

    if expr.is_Function:
        new_expr = expr.rewrite()
        if expr != new_expr:
            return const * expand_functions(new_expr, var)
        return const * new_expr

    terms = expr.as_ordered_terms()
    if len(terms) == 1:
        return const * expr

    expr = 0
    for term in terms:
        expr += expand_functions(term, var)

    return const * expr


def split_parens(s, delimiter=','):
    """Split a string by delimiter except if in ()"""

    parts = []
    bracket_level = 0
    current = []
    for c in (s + delimiter):
        if c == delimiter and bracket_level == 0:
            parts.append(''.join(current))
            current = []
        else:
            if c == '(':
                bracket_level += 1
            elif c == ')':
                bracket_level -= 1
            current.append(c)
    if bracket_level != 0:
        raise ValueError('Mismatched parentheses for ' + s)
    return parts
