"""This module defines the domain classes.

Copyright 2021 Michael Hayes, UCECE

"""

from .units import u as uu
from sympy import sqrt

class Domain(object):
    is_undefined_domain = False
    is_constant_domain = False
    is_time_domain = False
    is_laplace_domain = False
    is_fourier_domain = False
    is_angular_fourier_domain = False
    is_phasor_domain = False
    is_phasor_time_domain = False
    is_phasor_frequency_domain = False
    is_fourier_noise_domain = False
    is_angular_fourier_noise_domain = False    
    is_discrete_time_domain = False
    is_discrete_fourier_domain = False
    is_Z_domain = False
    is_transform_domain = False
    is_superposition_domain = False
    is_one_sided = False


class UndefinedDomain(Domain):
    domain = 'undefined'
    domain_label = 'Undefined'
    domain_units = 1
    is_undefined_domain = True


class ConstantDomain(Domain):    
    domain = 'constant'
    domain_label = 'Constant'    
    domain_units = 1
    is_constant_domain = True
    

class TimeDomain(Domain):
    domain = 'time'
    domain_label = 'Time'
    domain_units = uu.s
    is_time_domain = True

    
class LaplaceDomain(Domain):
    domain = 'laplace'        
    domain_label = 'Laplace frequency'
    domain_units = uu.rad / uu.s
    is_laplace_domain = True
    is_transform_domain = True

    
class FourierDomain(Domain):
    domain = 'fourier'    
    domain_label = 'Frequency'
    domain_units = uu.Hz
    is_fourier_domain = True
    is_transform_domain = True    

    
class AngularFourierDomain(Domain):
    domain = 'angular fourier'
    domain_label = 'Angular frequency'
    domain_units = uu.rad / uu.s
    is_angular_fourier_domain = True
    is_transform_domain = True    
    

class PhasorTimeDomain(Domain):
    domain = 'phasor'
    domain_label = ''
    domain_units = 1
    is_phasor_domain = True
    is_transform_domain = True

    
class PhasorFrequencyDomain(Domain):
    domain = 'phasor'
    domain_label = 'Angular Frequency'
    domain_units = uu.rad / uu.s
    is_phasor_domain = True
    is_transform_domain = True            


class FourierNoiseDomain(Domain):
    domain = 'fourier noise'    
    domain_label = 'Frequency'
    domain_units = uu.Hz
    is_fourier_noise_domain = True
    is_transform_domain = True
    is_one_sided = True


class AngularFourierNoiseDomain(Domain):
    domain = 'angular fourier noise'    
    domain_label = 'Angular frequency'
    domain_units = uu.rad / uu.s
    is_angular_fourier_noise_domain = True
    is_transform_domain = True
    is_one_sided = True    
    

class DiscreteTimeDomain(Domain):
    domain = 'discrete time'    
    domain_label = 'Discrete time'
    domain_units = 1
    is_discrete_time_domain = True
    is_transform_domain = False    
    

class DiscreteFourierDomain(Domain):
    domain = 'discrete fourier'
    domain_label = 'Discrete Fourier'    
    domain_units = 1
    is_discrete_fourier_domain = True
    is_transform_domain = True    
    
    
class ZDomain(Domain):
    domain = 'Z'
    domain_label = 'Z'    
    domain_units = 1
    is_Z_domain = True
    is_transform_domain = True    


class SuperpositionDomain(Domain):
    domain = 'superposition'
    domain_label = 'Superposition'    
    domain_units = 1
    is_superposition_domain = True
    is_transform_domain = True    
    
    
