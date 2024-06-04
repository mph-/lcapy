from sympy import sqrt, S, I

Zero = S.Zero
One = S.One


def _is_a_minus_b(expr):

    if not expr.is_Add:
        return False

    sub = 0
    add = 0
    for term in expr.as_ordered_terms():
        if ((term.is_Mul and term.args[0] < 0) or
            (term.is_constant() and term.is_negative)):
            sub += 1
        else:
            add += 1

    return add >= 1 and sub >= 1


def _decompose_bar(expr):

    # Look for scale * sqrt(dexpr) where dexpr = aexpr - bexpr

    dexpr = None
    scale = One

    for factor in expr.as_ordered_factors():
        if (factor.is_Pow and factor.args[1] == One / 2 and
            _is_a_minus_b(factor.args[0])):
            if dexpr is not None:
                return None
            dexpr = factor.args[0]
        else:
            scale *= factor
    if dexpr is None:
        return None

    return scale, dexpr


def _decompose_foo(expr):

    # Look for offset + scale * sqrt(dexpr) where dexpr = aexpr - bexpr

    offset = Zero
    scale = One
    dexpr = None

    for term in expr.as_ordered_terms():
        p = _decompose_bar(term)
        if p is None:
            offset += term
        elif dexpr is not None:
            return None
        else:
            scale, dexpr = p

    return offset, scale, dexpr


def _decompose(expr):

    scale = One
    scale2 = One
    offset = One
    dexpr = None

    for factor in expr.as_ordered_factors():
        if factor.is_Add:
            p = _decompose_foo(factor)
            if p is None:
                scale *= factor
            elif dexpr is not None:
                return None
            else:
                offset, scale2, dexpr = p
        else:
            scale *= factor

    if dexpr is None:
        dexpr = Zero
    return scale * offset, scale * scale2, dexpr


class QuadraticRoot:

    def __init__(self, offset, scale, dexpr, damping=None):

        # If the expression cannot be decomposed, then
        # scale = 1 and dexpr = 0.

        self.offset = offset
        self.scale = scale
        self.dexpr = dexpr
        self.damping = damping

    @property
    def expr(self):

        return self.offset + self.scale * sqrt(self.dexpr)

    def conjugate(self):

        if self.scale.is_imaginary:
            return self.offset - self.scale * sqrt(self.dexpr)
        else:
            return self.expr.conjugate()

    @classmethod
    def from_expr(cls, expr, damping=None):

        offset, scale, dexpr = _decompose(expr)

        if damping is None:
            pass
        elif damping == 'over':
            if scale.is_imaginary:
                dexpr = -dexpr
                scale /= I
        elif damping == 'critical':
            # This puts constraints on variables in dexpr.  For example,
            # if have sqrt(a - b) then a = b.
            pass
        elif damping == 'under':
            if not scale.is_imaginary:
                dexpr = -dexpr
                scale *= I
        else:
            raise ValueError('Unknown damping %s' % damping)

        if dexpr == 0:
            return None

        if damping == 'critical':
            dexpr = 0

        return cls(offset, scale, dexpr, damping)

    def is_conjugate_pair(self, d):

        if self.dexpr == 0 or d.dexpr == 0:
            return self.expr == d.expr.conjugate()

        if not self.scale.is_imaginary:
            return False

        if not d.scale.is_imaginary:
            return False

        return (self.offset == d.offset and self.scale == -d.scale and
                self.dexpr == d.dexpr)
