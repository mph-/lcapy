import sympy as sym

class Function(object):

    def __init__(self, arg):
        self.expr = arg
    
    def __call__(self, *args):

        cls = args[0].__class__

        tweak_args = list(args)
        for m, arg in enumerate(args):
            if isinstance(arg, (Expr, Function)):
                tweak_args[m] = arg.expr

        result = self.expr(*tweak_args)

        if isinstance(args[0], Expr):
            result = cls(result)

        for m, arg in enumerate(args[1:]):
            if isinstance(arg, (Expr, Function)):
                result = result.subs(tweak_args[m], arg)

        return result

    
class Log10(Function):

        def __call__(self, arg):
            return super(Log10, self).__call__(arg, 10)
    

conjugate = Function(sym.conjugate)
    
sqrt = Function(sym.sqrt)

exp = Function(sym.exp)

log = Function(sym.log)

log10 = Log10(sym.log)

sin = Function(sym.sin)

cos = Function(sym.cos)

tan = Function(sym.tan)

atan = Function(sym.atan)

atan2 = Function(sym.atan2)

sinh = Function(sym.sinh)

cosh = Function(sym.cosh)

tanh = Function(sym.tanh)

gcd = Function(sym.gcd)

u = H = Heaviside = Function(sym.Heaviside)

delta = DiracDelta = Function(sym.DiracDelta)


from .expr import Expr
