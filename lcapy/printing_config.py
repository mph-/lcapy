import sympy as sym

class PrintingConfig(object):

    def __init__(self):

        self._abbreviate_units = False

    @property
    def abbreviate_units(self):

        return self._abbreviate_units

    @abbreviate_units.setter
    def abbreviate_units(self, val):    

        self._abbreviate_units = val
    
        # Print abbreviated units, V not volt
        sym.printing.str.StrPrinter._default_settings['abbrev'] = val
        
