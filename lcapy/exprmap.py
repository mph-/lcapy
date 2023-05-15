from .domains import domains
from .quantities import quantities
from .exprclasses import exprclasses


def exprmap(quantity, domain):

    domaincls = domains[domain]

    if quantity == 'undefined':
        return exprclasses[domain][quantity]

    quantitycls = quantities[quantity]

    if (domaincls.is_constant_domain or domaincls.is_constant_time_domain):
        if quantitycls.is_ratio:
            domain = 'constant frequency response'

    if (domaincls.is_constant_domain or domaincls.is_constant_frequency_response_domain):
        if not quantitycls.is_ratio:
            domain = 'constant time'

    return exprclasses[domain][quantity]
