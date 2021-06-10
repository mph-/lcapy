from sympy import Eq, Symbol
from .utils import factor_const

class Transformer(object):

    name = 'undefined'
    is_bilateral = True
    
    def __init__(self):

        self.cache = {}
        self.expr = None

    def clear_cache(self):

        self.cache = {}        

    def error(self, message=''):
        if message != '':
            raise ValueError('Could not compute %s for %s' % (self.name, self.expr))
        raise ValueError('Could not compute %s for %s: %s' % (self.name, self.expr, message))        
    
    def simplify_term(self, expr):
        return expr

    def rewrite(self, expr):
        return expr    

    def transform(self, expr, var, conjvar, evaluate=True, **assumptions):

        # Squirrel away original expression for error messages
        self.expr = expr
        
        if expr.is_Equality:
            return Eq(self.transform(expr.args[0], var, conjvar),
                      self.transform(expr.args[1], var, conjvar))        

        # Handle Lcapy symbols for the variable and conjugate variable
        try:
            var = var.expr
        except:
            pass

        try:
            conjvar = conjvar.expr
        except:
            pass        

        self.check(expr, var, conjvar, **assumptions)

        # Unilateral transforms ignore expr for t < 0 so remove Piecewise.
        if not self.is_bilateral:
            if expr.is_Piecewise and expr.args[0].args[1].has(var >= 0):
                expr = expr.args[0].args[0]
        
        if not evaluate:
            return self.noevaluate(expr, var, conjvar)

        const, expr = factor_const(expr, var)

        # TODO, apply similarity theorem.
        
        key = self.key(expr, var, conjvar, **assumptions)
        if key in self.cache:
            return const * self.cache[key]

        # The variable may have been created with different attributes,
        # say when using sym.sympify('DiracDelta(t)') since this will
        # default to assuming that t is complex.  So if the symbol has the
        # same representation, convert to the desired one.
        svar = Symbol(str(var))
        expr = expr.replace(svar, var)

        orig_expr = expr
        expr = self.rewrite(expr)        

        terms = expr.expand().as_ordered_terms()
        result = 0

        try:
            for term in terms:
                sterm = self.simplify_term(term)
                result += self.term(sterm, var, conjvar)
        except ValueError:
            self.error()

        self.cache[key] = result
        return const * result

    