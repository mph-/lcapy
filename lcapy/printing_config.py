import sympy as sym
from sympy.printing.str import StrPrinter
from sympy.printing.latex import LatexPrinter
from sympy.printing.pretty.pretty import PrettyPrinter


class PrintingConfig(object):

    def __init__(self):

        self._abbreviate_units = False
        self._order = None

    @property
    def abbreviate_units(self):

        return self._abbreviate_units

    @abbreviate_units.setter
    def abbreviate_units(self, val):

        self._abbreviate_units = val

        # Print abbreviated units, V not volt
        StrPrinter._default_settings['abbrev'] = val

    @property
    def order(self):

        return self._order

    @order.setter
    def order(self, val):
        """Use 'none' to disable symbol ordering when expressions are printed."""

        self._order = val

        StrPrinter._default_settings['order'] = val
        LatexPrinter._default_settings['order'] = val
        PrettyPrinter._default_settings['order'] = val
