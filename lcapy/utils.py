import sympy as sym

def factor_const(expr, var):
    """Extract constant factor from expression and return tuple
    of constant and the rest of the expression."""    

    # Perhaps use expr.as_coeff_Mul() ?

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
    
    if not expr.as_poly(var).is_linear:
        raise ValueError('Expression not a linear function of %s: %s' % (var, expr))

    scale = expr.coeff(var, 1)
    shift = expr.coeff(var, 0)        

    return scale, shift


def as_N_D(expr, var, monic_denominator=False):

    N = 1
    D = 1
    factors = expr.as_ordered_factors()
    
    for factor in factors:
        a, b = factor.as_numer_denom()
        N *= a
        if b.is_polynomial(var):
            D *= b
        else:
            N /= b
                
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

    return [term / D for term in N.expand().as_ordered_terms ()]


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


def separate_dirac_delta(expr):

    terms = expr.as_ordered_terms()
    deltas = []
    rest = 0

    for term in terms:
        if term.has(sym.DiracDelta):
            deltas.append(term)
        else:
            rest += term
            
    return rest, deltas


def remove_images(expr, var, dt, m1=0, m2=0):

    const, expr1 = factor_const(expr, var)

    result = sym.S.One
    terms = expr1.as_ordered_terms()

    if len(terms) > 1:
        for term in expr1.as_ordered_terms():
            result += remove_images(term, var, dt, m1, m2)
        return const * result
        
    if not isinstance(expr1, sym.Sum):
        return expr
    sumsym = expr1.args[1].args[0]
    foo = var - sumsym / dt
    if not expr1.args[0].has(foo):
        return expr    

    if m1 == 0 and m2 == 0:
        return const * expr1.args[0].replace(foo, var)

    return const * sym.Sum(expr1.args[0], (expr1.args[1][0], m1, m2))


