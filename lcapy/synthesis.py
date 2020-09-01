"""
This module implements experimental network synthesis.

Copyright 2020 Michael Hayes, UCECE
"""

from .oneport import L, C, R, G
from .sexpr import s, Zs


# Could make this a mixin for Zs and Ys but would need a Z flavour
# and a Y flavour.
class Synthesis(object):

    def seriesRL(self, lexpr):
        """Z = s * L + R"""
        
        if lexpr == 0:
            # Perhaps return Wire object
            return None

        var = lexpr.var
        d = lexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, None)
        if a is not None:
            net = R(a)

        a = d.pop(var, None)
        if a is not None:
            if net is None:
                net = L(a)
            else:
                net = net + L(a)            

        if d != {}:
            raise ValueError('Not series RL')
        
        return net

    def seriesRC(self, lexpr):
        """Z = R + 1 / (s * C)"""

        if lexpr == 0:
            raise ValueError('Not series RC')

        var = lexpr.var
        d = lexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, None)
        if a is not None:
            net = R(a)

        a = d.pop(1 / var, None)
        if a is not None:
            if net is None:
                net = C(1 / a)
            else:
                net = net + C(1 / a)            

        if d != {}:
            raise ValueError('Not series RC')

        return net

    def seriesLC(self, lexpr):
        """Z = s * L + 1 / (s * C)"""

        if lexpr == 0:
            raise ValueError('Not series LC')

        var = lexpr.var
        d = lexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, 0)
        if a != 0:
            raise ValueError('Not series LC')            

        a = d.pop(1 / var, None)
        if a is not None:
            net = C(1 / a)

        a = d.pop(var, None)
        if a is not None:
            if net is None:
                net = L(a)
            else:
                net = net + L(a)                            

        if d != {}:
            raise ValueError('Not series LC')

        return net

    def seriesRLC(self, lexpr):
        """Z = s * L + R + 1 / (s * C)"""

        if lexpr == 0:
            raise ValueError('Not series RLC')

        var = lexpr.var
        d = lexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, None)
        if a is not None:
            net = R(a)

        a = d.pop(1 / var, None)
        if a is not None:
            if net is None:
                net = C(1 / a)
            else:
                net = net + C(1 / a)

        a = d.pop(var, None)
        if a is not None:
            if net is None:
                net = L(a)
            else:
                net = net + L(a)                            

        if d != {}:
            raise ValueError('Not series RLC')

        return net        

    def FosterI(self, lexpr):

        expr = lexpr.partfrac(combine_conjugates=True)

        # TODO
        
    def CauerI(self, lexpr):

        # This can result in negative valued components.
        
        coeffs = lexpr.continuous_fraction_coeffs()

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

    def default(self, lexpr):

        return self.CauerI(lexpr)
    
    def network(self, lexpr, form='default'):

        if form == 'default':
            form = 'CauerI'

        if not isinstance(lexpr, Zs):
            raise ValueError('Expression needs to be Zs object')

        # Should test if a positive real function.
        
        forms = {'CauerI': self.CauerI,
                 'FosterI': self.FosterI,
                 'seriesRL': self.seriesRL,
                 'seriesRC': self.seriesRC,
                 'seriesLC': self.seriesLC,                 
                 'seriesRLC': self.seriesRLC}

        if form not in forms:
            raise ValueError('Unknown form %s, known forms: %s' % forms.keys())
        return forms[form](lexpr)

    
def network(lexpr, form='default'):        

    return Synthesis().network(lexpr, form)
