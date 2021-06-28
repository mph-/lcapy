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

    if m2 == 0 and isinstance(m1, tuple) and len(m1) == 2:
        # Perhaps should warn that this might be deprecated?
        m1, m2 = m1

    remove_all =  m1 == 0 and m2 == 0
        
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
