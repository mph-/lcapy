from sympy import sqrt, S, I

Zero = S.Zero
One = S.One

def _decompose_bar(expr):

    # Look for scale * sqrt(dexpr) where dexpr = aexpr - bexpr

    dexpr = None
    scale = One

    for factor in expr.as_ordered_factors():
        if (factor.is_Pow and factor.args[1] == One / 2 and
            factor.args[0].is_Add and
            (factor.args[0].args[1].is_Mul and
             (factor.args[0].args[1].args[0] < 0) or
             (factor.args[0].args[0].is_Mul and
              factor.args[0].args[0].args[0] < 0))):
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
    return scale, offset, scale2, dexpr


class QuadraticRoot:

    def __init__(self, scale, offset, scale2, dexpr, damping=None):

        # If the expression cannot be decomposed, then
        # scale = scale2 = 1 and dexpr = 0.

        self.offset = offset
        self.scale = scale
        self.scale2 = scale2
        self.dexpr = dexpr
        self.damping = damping

    @property
    def expr(self):

        return self.scale * (self.offset + self.scale2 * sqrt(self.dexpr))

    def conjugate(self):

        if self.scale2.is_imaginary:
            return self.scale * (self.offset - self.scale2 * sqrt(self.dexpr))
        else:
            return self.expr.conjugate()

    @classmethod
    def from_expr(cls, expr, damping=None):

        scale, offset, scale2, dexpr = _decompose(expr)

        if damping in (None, 'over'):
            pass
        elif damping == 'critical':
            # This puts constraints on variables in dexpr.
            dexpr = 0
        elif damping == 'under':
            if not scale2.is_imaginary:
                dexpr = -dexpr
                scale2 *= I
        else:
            raise ValueError('Unknown damping %s' % damping)

        if dexpr == 0:
            return None

        return cls(scale, offset, scale2, dexpr, damping)

    def is_conjugate_pair(self, d):

        if self.dexpr == 0 or d.dexpr == 0:
            return self.expr == d.expr.conjugate()

        if not self.scale2.is_imaginary:
            return False

        if not d.scale2.is_imaginary:
            return False

        return (self.scale == d.scale and self.offset == d.offset and
                self.scale2 == -d.scale2 and self.dexpr == d.dexpr)
