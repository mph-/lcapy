"""
This module implements experimental network synthesis.

Copyright 2020 Michael Hayes, UCECE
"""

from .oneport import L, C, R, G


class Synthesis(object):

    def __init__(self, expr):
        self.expr = expr
    

class CauerIForm(Synthesis):

    def __call__(self):

        # This can result in negative valued components.
        
        coeffs = self.expr.continuous_fraction_coeffs()

        def series_net(cls1, a1, cls0, a0):

            if a1 == 0 and a0 == 0:
                return None
            if a1 == 0:
                return cls0(a0)
            if a0 == 0:
                return cls1(a1)
            return cls1(a1) + cls0(a0)

        net = None
        for m, coeff in enumerate(reversed(coeffs)):

            n = len(coeffs) - m - 1
            
            parts = coeff.coeffs()
            if len(parts) > 2:
                raise ValueError('Cannot express %s as a network' % coeff)
            
            if n & 1 == 0:
                if len(parts) == 2:
                    net1 = series_net(L, parts[0], R, parts[1])
                    if net1 is None:
                        continue
                else:
                    if parts[0] == 0:
                        continue
                    net1 = R(parts[0])
                
                if net is None:
                    net = net1
                else:
                    net = net1 + net
            else:
                if len(parts) == 2:                
                    net1 = series_net(C, parts[0], G, parts[1])
                    if net1 is None:
                        continue
                else:
                    if parts[0] == 0:
                        continue
                    net1 = G(parts[0])                    

                if net is None:
                    net = net1
                else:
                    net = net1 | net
        return net
    

def synthesise(expr, form='CauerI'):

    # Perhaps convert to s-domain?

    forms = {'CauerI': CauerIForm}

    if form not in forms:
        raise ValueError('Unknown form %s, known forms: %s' % forms.keys())
    return forms[form](expr)()
