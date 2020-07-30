"""This module performs network transformations.

Copyright 2020 Michael Hayes, UCECE

"""


def Z_wye_to_delta(Z1, Z2, Z3):
    """Perform wye to delta transformation of three impedances.

    This is equivalent to a tee-pi transform or a star-mesh transform."""

    N = Z1 * Z2 + Z1 * Z3 + Z2 * Z3
    
    Za = N / Z1
    Zb = N / Z2
    Zc = N / Z3

    return Za, Zb, Zc


def Z_delta_to_wye(Za, Zb, Zc):
    """Perform wye to delta transformation of three impedances.

    This is equivalent to a pi-tee transform or a mesh-star transform."""    

    D = Za + Zb + Zc
    
    Z1 = Zb * Zc / D
    Z2 = Za * Zc / D
    Z3 = Za * Zb / D

    return Z1, Z2, Z3


def Y_delta_to_wye(Ya, Yb, Yc):
    """Perform wye to delta transformation of three admittances.
        
    This is equivalent to a tee-pi transform or a star-mesh transform."""

    N = Ya * Yb + Ya * Yc + Yb * Yc
    
    Y1 = N / Ya
    Y2 = N / Yb
    Y3 = N / Yc

    return Y1, Y2, Y3


def Y_wye_to_delta(Y1, Y2, Y3):
    """Perform wye to delta transformation of three admittances.

    This is equivalent to a pi-tee transform or a mesh-star transform."""    

    D = Y1 + Y2 + Y3
    
    Ya = Y2 * Y3 / D
    Yb = Y1 * Y3 / D
    Yc = Y1 * Y2 / D

    return Ya, Yb, Yv




