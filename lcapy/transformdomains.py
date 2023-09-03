"""This module provides the TransformDomains class.

Copyright 2023 Michael Hayes, UCECE

"""


from .expr import Expr
from .symbols import omega


class TransformDomains(dict):

    def __getattr__(self, attr):
        if attr not in self:
            raise AttributeError('Unknown attribute %s' % attr)
        return self[attr]

    def __getitem__(self, key):
        if key == 'w':
            key = omega
        # This allows a[omega] to work if omega used as key
        # instead of 'omega'.
        if isinstance(key, Expr):
            key = key.expr
        return super(TransformDomains, self).__getitem__(key)
