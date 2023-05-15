from .domains import domains
from .quantities import quantities
from .seqclasses import seqclasses


def seqmap(quantity, domain):

    domaincls = domains[domain]

    if quantity == 'undefined':
        return seqclasses[domain][quantity]

    return seqclasses[domain][quantity]
