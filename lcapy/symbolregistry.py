"""This module maintains a registry of symbols.

Copyright 2022 Michael Hayes, UCECE

"""

from sympy import Symbol
from .attrdict import AttrDict
from warnings import warn

# Cannot add attributes to SymPy symbols since the Symbol
# class defines __slots__ = ().  So perhaps use facade for
# symbols that includes the symbol kind?
kinds = {}


class SymbolRegistry(AttrDict):
    """This maintains a registry of symbols so that there is only
    a single symbol of a given name.

    The symbol kinds are:
    'domain' defined by Lcapy for domain variables
    'misc' defined by Lcapy for misc. variables, such as nu
    'user' defined by the user using `symbol()` or `symbols()`
    'expr' defined by the user using `expr()`
    """

    def register(self, name, symbol, kind='user'):
        from .state import state

        if state.notify_symbol_add:
            print("Adding symbol '%s'" % name)

        self[name] = symbol
        kinds[name] = kind
        return symbol

    def add(self, name: str, kind='user', override=True,
            force=False, **assumptions):

        if name in self:
            symbol = self[name]

            if not override:
                return symbol

            if kinds[name] == 'domain' and not force:
                raise ValueError(
                    'Cannot override domain symbol %s without force=True' % name)

        if assumptions == {}:
            assumptions['positive'] = True
            # Note this implies that imag is False.   Also note that all
            # reals are considered complex (but with a zero imag part).

        elif 'positive' in assumptions:
            if not assumptions['positive']:
                assumptions.pop('positive')

        symbol = Symbol(name, **assumptions)
        return self.register(name, symbol, kind)

    def delete(self, name):

        self.pop(name)

    def add_domain(self, name: str, **assumptions):

        return self.add(name, kind='domain', **assumptions)

    def add_misc(self, name: str, **assumptions):

        return self.add(name, kind='misc', **assumptions)

    def add_user(self, name: str, **assumptions):

        return self.add(name, kind='user', **assumptions)

    def add_from_expr(self, expr):
        """Add all the symbols in an expression."""

        for symbol in expr.atoms(Symbol):
            name = str(symbol)
            self.register(name, symbol, kind='expr')

    def delete_user(self):
        """Delete user symbols."""

        for name in self.user().keys():
            self.delete(name)

    def by_kind(self, kind):

        syms = self.__class__()
        for name, sym in self.items():
            if kinds[name] == kind:
                syms[name] = sym
        return syms

    def user(self):
        """Return user defined symbols."""

        return self.by_kind('user')

    def expr(self):
        """Return expr defined symbols."""

        return self.by_kind('expr')

    def delete_expr(self):
        """Delete expr symbols."""

        for name in self.expr().keys():
            self.delete(name)

    def clear(self):

        self.delete_expr()
        self.delete_user()

    def lookup(self, name):

        new = name
        if not isinstance(name, str):
            name = str(name)

        # Replace symbol names with symbol definitions to
        # avoid problems with real or positive attributes.
        try:
            return self[name]
        except:
            from .state import state

            # This can occur when trying to substitute an expression
            # rather than a symbol, for example, when making state-space.
            if state.warn_unknown_symbol:
                warn('Lcapy does not know symbol %s' % name)
            if state.break_unknown_symbol:
                import pdb
                pdb.set_trace()

        return new

    def kind(self, name):
        """Return the symbol kind or 'unknown' if unknown."""

        try:
            return kinds[name]
        except:
            return 'unknown'
