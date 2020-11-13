"""This module provides the Super class.  It is the base class for
Current and Voltage.  It represents voltages and currents as a
superposition in different transform domains.

Copyright 2019--2020 Michael Hayes, UCECE

"""

from __future__ import division
from .expr import Expr, ExprDict, expr
from .sym import tsym, omega0sym, symbols_find, is_sympy, symsymbol
from .acdc import is_ac
from .printing import pprint, pretty, latex
import six

__all__ = ('Super', 'Voltage', 'Current')

class Super(ExprDict):
    """This class represents a superposition of different signal types,
    DC, AC, transient, and noise.
    
    The time-domain representation is returned with the time method,
    V.time(), or with the notation V(t).  This does not include the
    noise component.

    The Laplace representation is returned with the laplace method,
    V.laplace() or with the notation V(s).  This does not include the
    noise component.

    Noise components with different noise identifiers are stored
    separately, keyed by the noise identifier.  They are ignored by
    the laplace() and time() methods.

    The total noise can be accessed with the .n attribute.  This sums
    each of the noise components in quadrature since they are
    independent.

    For example, consider V = Voltage('cos(3 * t) + exp(-4 * t) + 5')

    str(V(t)) gives 'cos(3*t) + 5 + exp(-4*t)'

    str(V(s)) gives 's/(9*(s**2/9 + 1)) + 1/(4*(s/4 + 1)) + 5/s'

    V.dc gives 5

    V.ac gives {3: 1}

    V.s gives 1/(s + 4)
    """

    # TODO: rework this class to better show 0 result.
    
    # Where possible this class represents a signal in the time-domain.
    # It can decompose a signal into AC, DC, and transient components.
    # The 't' key is the transient component viewed in the time domain.
    # The 's' key is the transient component viewed in the Laplace domain.    
    
    def __init__(self, *args, **kwargs):
        super (Super, self).__init__()

        for arg in args:
            self.add(arg)

    def _representation(self):
        if not any(self):
            return 0
        if False:
            return self
        # It is probably less confusing for a user to display
        # using a decomposition in the transform domains.
        # We could present the result in the time-domain but this
        # hides the underlying way the signal is analysed.
        return self.decompose()

    def latex(self, **kwargs):
        """Latex"""
        return latex(self._representation(), **kwargs)
    
    def _repr_pretty_(self, p, cycle):
        """This is used by jupyter notebooks to display an expression using
        unicode."""

        p.text(pretty(self._representation()))

    def _repr_latex_(self):
        """This is used by jupyter notebooks to display an expression using
        LaTeX markup.  However, this requires matjax.  If this method
        is not defined, jupyter falls back on _repr__pretty_ which
        outputs unicode."""

        return '$$' + latex(self._representation()) + '$$'

    def pprint(self):
        """Pretty print"""

        return pprint(self._representation())

    def __getitem__(self, key):
        # This allows a[omega] to work if omega used as key
        # instead of 'omega'.
        if isinstance(key, Expr):
            key = key.expr
        return super(Super, self).__getitem__(key)

    @property
    def symbols(self):
        """Return dictionary of symbols in the expression keyed by name."""

        syms = {}
        for expr in self.values():
            syms.update(expr.symbols)
        return syms

    def ac_keys(self):
        """Return list of keys for all ac components."""

        keys = []
        for key in self.decompose().keys():
            if not isinstance(key, str) or key == 'w':
                keys.append(key)
        return keys

    def noise_keys(self):
        """Return list of keys for all noise components."""

        keys = []
        for key in self.keys():
            if isinstance(key, str) and key[0] == 'n':
                keys.append(key)
        return keys    

    def has(self, subexpr):
        """Test whether the sub-expression is contained.  For example,
         V.has(exp(t)) 
         V.has(t)

        """        
        for key,expr  in self.items():
            if expr.has(subexpr):
                return True
        return False

    def has_symbol(self, sym):
        """Test if have symbol.  For example,
        V.has_symbol('a')
        V.has_symbol(t)
        
        """                        
        return self.has(symsymbol(sym))
    
    @property
    def has_dc(self):
        """True if there is a DC component."""                
        return self == 0 or 'dc' in self.decompose()

    @property
    def has_ac(self):
        """True if there is an AC component."""        
        return self.ac_keys() != []

    @property
    def has_s_transient(self):
        """True if have transient component defined in the s-domain."""
        return 's' in self

    @property
    def has_t_transient(self):
        """True if have transient component defined in the time-domain."""
        return 't' in self

    @property
    def has_transient(self):
        """True if have transient component."""        
        return self.has_s_transient or self.has_t_transient

    @property
    def has_noisy(self):
        """True if there is a noise component."""                
        return self.noise_keys() != []

    @property
    def is_dc(self):
        """True if only has a DC component."""                
        return self == 0 or (self.has_dc and list(self.decompose().keys()) == ['dc'])

    @property
    def is_ac(self):
        """True if only has AC components."""                        
        return self.has_ac and self.ac_keys() == list(self.keys())

    @property
    def is_noisy(self):
        """True if only has noise components."""                                
        return self.has_noisy and self.noise_keys() == list(self.keys())

    @property
    def is_s_transient(self):
        """True if only has s-domain transient component."""
        return list(self.decompose().keys()) == ['s']

    @property
    def is_t_transient(self):
        """True if only has t-domain transient component."""
        return list(self.decompose().keys()) == ['t']    

    @property
    def is_transient(self):
        """True if only has transient component(s).  Note, should
        not have both t and s components."""
        return self.is_s_transient or self.is_t_transient
    
    @property
    def is_causal(self):
        return (self.is_transient and self.s.is_causal) or self == 0

    @property
    def is_superposition(self):
        return len(self.keys()) > 1

    def __call__(self, arg, **assumptions):
        """
        arg determines the returned representation
          t: time-domain representation
          s: Laplace domain representation
          f: Fourier representation
          omega: Fourier representation (angular frequency)
          jomega: Laplace domain representation with s = jomega

        For example, V(t), V(f), or V(2 * t).

        If arg is a constant, the time-domain representation
        evaluated at the argument is returned, for example,
        V(0) returns the dc value.
        """

        from .transform import call
        return call(self, arg, **assumptions)        

    def subs(self, *args, **kwargs):

        new = self.__class__()
        for kind, value in self.items():
            new[kind] = value.subs(*args, **kwargs)
        return new        
    
    def initial_value(self):
        """Determine value at t = 0. 
        See also pre_initial_value and post_initial_value"""

        return self.time().initial_value()

    def pre_initial_value(self):
        """Determine value at t = 0-.
        See also initial_value and post_initial_value"""

        return self.time().pre_initial_value()

    def post_initial_value(self):
        """Determine value at t = 0+.
        See also pre_initial_value and initial_value"""

        return self.time().post_initial_value()

    def final_value(self):
        """Determine value at t = oo."""

        return self.time().final_value()

    def transform(self, arg, **assumptions):
        """Transform into a different domain."""        

        from .transform import transform
        return transform(self, arg, **assumptions)

    def __add__(self, x):

        def _is_s_arg(x):
            return isinstance(x, sExpr) or (isinstance(x, Super) and 's' in x)

        if _is_s_arg(x):
            new = self.decompose()
        else:
            new = self.__class__(self)            

        if isinstance(x, Super):
            for value in x.values():
                new.add(value)
        else:
            new.add(x)
        return new

    def __radd__(self, x):
        return self.__add__(x)

    def __sub__(self, x):
        return -x + self

    def __rsub__(self, x):
        return -self + x

    def __neg__(self):
        new = self.__class__()
        for kind, value in self.items():
            new[kind] = -value
        return new

    def __scale__(self, x):
        new = self.__class__()
        for kind, value in self.items():
            new[kind] = value * x
        return new

    def __eq__(self, x):
        """ Test for mathematical equality as far as possible.
        This cannot be guaranteed since it depends on simplification.
        Note, SymPy comparison is for structural equality.  

        Also note, noise cannot be compared."""
        
        diff = (self - x).simplify().decompose()
        for kind, value in diff.items():
            if value != 0:
                return False
        return True

    def __ne__(self, x):
        """ Test for mathematical inequality as far as possible.
        This cannot be guaranteed since it depends on simplification.
        Note, SymPy comparison is for structural equality.  

        Also note, noise cannot be compared."""
        
        diff = (self - x).simplify().decompose()
        for kind, value in diff.items():
            if value != 0:
                return True
        return False    

    def _decompose(self, expr):
        """Decompose a t-domain expr into AC, DC, and s-domain
        transient components."""

        # Extract DC components
        dc = expr.expr.coeff(tsym, 0)
        if dc != 0:
            self['dc'] = cExpr(dc)
            expr -= dc

        if expr == 0:
            return

        # Extract AC components        
        ac = 0
        terms = expr.expr.as_ordered_terms()
        for term in terms:
            if is_ac(term, tsym):
                self.add(tExpr(term).phasor())
                ac += term

        expr -= ac
        if expr == 0:
            return

        # The remaining components are considered transient
        # so convert to Laplace representation.
        sval = expr.laplace()

        self['s'] = sval

    def decompose(self):
        """Decompose into a new representation in the transform domains."""
        
        if hasattr(self, '_decomposition'):
            return self._decomposition

        new = self.__class__()
        if 't' in self:
            new._decompose(self['t'])
        for kind, value in self.items():
            if kind != 't':
                new.add(value)
        self._decomposition = new                
        return new
    
    def select(self, kind):
        """Select a component of the signal representation by kind where:
        'super' : the entire superposition
        'time' :  the time domain representation (equivalent to self.time())
        'laplace' :  the laplace domain representation (equivalent to self.laplace())
        'ivp' :  the s-domain representation (equivalent to self.laplace())
        'dc' : the DC component
        'ac' : the AC component with angular frequency omega0
        omega : the AC component with angular frequency omega
        's' : the transient component in the s-domain
        'n' : the noise component
        't' : the time-domain transient component (this may or may not
              include the DC and AC components).

        """
        if kind == 'super':
            return self
        elif kind == 'time':
            return self.time()
        elif kind in ('ivp', 'laplace'):
            return self.laplace()

        if isinstance(kind, str) and kind[0] == 'n':
            if kind not in self:
                return self.decompose_domains['n'](0)
            return self[kind]
        
        obj = self
        if 't' in self and (kind == omega or 't' != kind):
            # The rationale here is that there may be
            # DC and AC components included in the 't' part.
            obj = self.decompose()
            
        if kind not in obj:
            if kind not in obj.decompose_domains:
                kind = 'ac'
            return obj.decompose_domains[kind](0)
        return obj[kind]

    def netval(self, kind):

        def kind_keyword(kind):
            if isinstance(kind, str) and kind[0] == 'n':
                return 'noise'
            elif kind in ('ivp', 'laplace'):
                return 's'
            elif kind in ('t', 'time'):
                return ''                
            elif not isinstance(kind, str):
                return 'ac'
            return kind

        val = self.select(kind)
        if kind in ('s', 'ivp') and (val.is_causal or val.is_dc or val.is_ac):
            # Convert to time representation so that can re-infer
            # causality, etc.
            return val.time()

        keyword = kind_keyword(kind)
        
        if 'nid' in val.assumptions:
            return keyword, val, val.nid

        if keyword == 'ac':
            return keyword, val, 0, val.omega

        return keyword, val
    
    def _kind(self, value):
        if isinstance(value, Phasor):
            # Use angular frequency for key.  This can be a nuisance
            # for numerical values since cannot use x.3 syntax
            # say for an angular frequency of 3.
            key = value.omega
            if isinstance(key, Expr):
                key = key.expr
            return key

        for kind, mtype in self.decompose_domains.items():
            if isinstance(value, mtype):
                return kind
        return None

    def _parse(self, string):
        """Parse t or s-domain expression or symbol, interpreted in time
        domain if not containing s, omega, or f.  Return most
        appropriate transform domain.

        """

        symbols = symbols_find(string)
        if 's' in symbols:
            return self.add(sExpr(string))

        if 'omega' in symbols:
            if 't' in symbols:
                # Handle cos(omega * t)
                return self.add(tExpr(string))
            return self.add(omegaExpr(string))

        if 'f' in symbols:
            # TODO, handle AC of different frequency
            return self.add(fExpr(string))

        return self.add(tExpr(string))

    def _add_noise(self, value):

        if value.nid not in self:
            self[value.nid] = value
        else:
            self[value.nid] += value
            if self[value.nid] == 0:
                self.pop(value.nid)
    
    def add(self, value):
        """Add a value into the superposition."""

        # Avoid triggering __eq__ for Super otherwise have infinite recursion
        if not isinstance(value, Super):
            try:
                val = value.expr
            except:
                val = value
            if val == 0:
                return

        if '_decomposition' in self:
            delattr(self, '_decomposition')

        if isinstance(value, Super):
            for kind, value in value.items():
                self.add(value)
            return

        if isinstance(value, six.string_types):
            return self._parse(value)

        if isinstance(value, noiseExpr):
            self._add_noise(value)
            return
        
        # DC should be real but allow const complex value.
        if isinstance(value, (int, float, complex)):
            value = self.decompose_domains['dc'](value)
        elif is_sympy(value):
            try:
                # Look for I, 5 * I, etc.
                if value.is_constant:
                    value = self.decompose_domains['dc'](value)
            except:
                pass

        # TODO, perhaps handle Fourier domain expressions in the
        # decomposition?  For now, convert to time domain.
        if isinstance(value, (fExpr, omegaExpr)):
            value = value.time()

        kind = self._kind(value)
        if kind is None:
            if value.__class__ in self.type_map:
                value = self.type_map[value.__class__](value)
            else:
                for cls1, cls2 in self.type_map.items():
                    if isinstance(value, cls1):
                        value = cls2(value)
                        break
            kind = self._kind(value)

        if kind is None:
            raise ValueError('Cannot handle value %s of type %s' %
                             (value, type(value).__name__))

        if kind not in self:
            self[kind] = value
        else:
            self[kind] += value

    @property
    def dc(self):
        """Return the DC component."""        
        return self.select('dc')

    @property
    def ac(self):
        """Return the AC components."""                
        if 't' in self.keys():
            self = self.decompose()        
        return ExprDict({k: v for k, v in self.items() if k in self.ac_keys()})

    @property
    def transient(self):
        """Return the transient component."""        
        return self.select('s').time()
    
    @property
    def s(self):
        """Return the s-domain representation of the transient component.
        This is not the full s-domain representation as returned by the
        laplace method V.laplace() or V(s).

        This attribute may be deprecated due to possible confusion."""
        return self.select('s')

    @property
    def n(self):
        result = self.decompose_domains['n'](0)
        for key in self.noise_keys():
            result += self[key]
        return result

    @property
    def noise(self):
        """Return the total noise."""
        return self.n

    @property
    def w(self):
        """Return the AC component with angular frequency omega0.
        This should be deprecated."""
        return self.select(omega0sym)    

    def time(self, **assumptions):
        """Convert to time domain."""

        result = self.time_class(0)

        # TODO, integrate noise
        for val in self.values():
            if hasattr(val, 'time'):
                result += val.time(**assumptions)
            else:
                result += val
        return result

    def transient_response(self, tvector=None):
        """Evaluate transient (impulse) response."""

        texpr = self.time()

        if tvector is None:
            return texpr

        return texpr.evaluate(tvector)
    
    def frequency_response(self, fvector=None):
        """Convert to frequency domain and evaluate response if frequency
        vector specified.

        """

        X = self.fourier()

        if fvector is None:
            return X

        return X.evaluate(fvector)

    def laplace(self, **assumptions):
        """Convert to s-domain."""                

        result = self.laplace_class(0)
        for val in self.values():
            result += val.laplace(**assumptions)
        return result

    def fourier(self, **assumptions):
        """Convert to Fourier domain."""        

        # TODO, could optimise.
        return self.time(**assumptions).fourier(**assumptions)

    def canonical(self):
        new = self.__class__()
        for kind, value in self.items():
            new[kind] = value.canonical()

        return new

    def simplify(self):
        new = self.__class__()
        for kind, value in self.items():
            new[kind] = value.simplify()

        return new

    def oneport(self):
        """Create oneport component."""
        return self.cpt()
    
from .cexpr import cExpr        
from .fexpr import fExpr    
from .sexpr import sExpr
from .texpr import tExpr
from .noiseexpr import noiseExpr
from .phasor import Phasor
from .omegaexpr import omegaExpr
from .symbols import s, omega
