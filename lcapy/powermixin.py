from .quantity import Quantity

class PowerMixin(Quantity):

    quantity = 'power'
    quantity_label = 'Power'
    quantity_units = 'W'
    is_power = True

    # TODO: the units are not power in some of the transform domains.
    # For example, in the Fourier domain V * I has units of V/Hz * A/Hz
    # which is equivalent to W / Hz^2 = J / Hz.   This is an energy spectral
    # density.
    
    
