"""This module provides a factory function for schematic placement algorithms.

Copyright 2021 Michael Hayes, UCECE
"""


from .schemgraphplacer import SchemGraphPlacer
from .schemlineqplacer import SchemLineqPlacer

def schemplacer(elements, nodes, method='graph', debug=False):

    if method == 'graph':
        placer = SchemGraphPlacer(elements, nodes, debug)
    elif method == 'lineq':
        placer = SchemLineqPlacer(elements, nodes, debug)
    else:
        raise ValueError('Unknown placer method %s' % method)

    return placer

    
