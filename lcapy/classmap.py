from .cexpr import ConstantDomainExpression
from .fexpr import FourierDomainExpression
from .sexpr import LaplaceDomainExpression
from .texpr import TimeDomainExpression
from .omegaexpr import AngularFourierDomainExpression
from .phasor import PhasorTimeDomainExpression
from .noiseomegaexpr import AngularFourierNoiseDomainExpression


classmap = {'s': LaplaceDomainExpression,
            'ivp': LaplaceDomainExpression,
            'laplace': LaplaceDomainExpression,
            't': TimeDomainExpression,
            'time': TimeDomainExpression,            
            'f': FourierDomainExpression,
            'dc': ConstantDomainExpression,
            'n': AngularFourierNoiseDomainExpression}

symbolmap = {'s': '(s)',
             'laplace': '(s)',
             't': '(t)',
             'time': '(t)',
             'f': '(f)',
             'dc': ''}

domainmap = {'constant': 'dc',
             'time': 't',
             'laplace': 's',
             'fourier': 'f',
             'angular fourier': 'omega'}

def domain_kind_to_class(kind):

    if isinstance(kind, str) and kind[0] == 'n':
        return classmap['n']
    
    try:
        return classmap[kind]
    except:
        return PhasorTimeDomainExpression


def domain_kind_quantity_to_class(kind, quantity='undefined'):

    cls = domain_kind_to_class(kind)
    # TODO: There must be a better way than creating a dummy object.
    return cls(0)._class_by_quantity(quantity)


def domain_quantity_to_class(domain, quantity='undefined'):

    kind = domainmap[domain]
    return domain_kind_quantity_to_class(kind, quantity)


def domain_kind_to_symbol(kind, name):

    try:
        return name + symbolmap[kind]
    except:
        return name + '(%s)' % kind
    

def domain_to_domain_kind(domain):
    
    try:
        return domainmap[domain]
    except:
        return None


def expr_to_domain_kind(expr):
    
    try:
        return domainmap[expr.domain]
    except:
        return None
    
