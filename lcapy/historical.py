# These are historical classnames and will be deprecated.

from .cexpr import ConstantExpression as cExpr
from .cexpr import ConstantVoltage as Vc
from .cexpr import ConstantCurrent as Ic

from .fexpr import FourierDomainExpression as fExpr
from .fexpr import FourierDomainTransferFunction as Hf
from .fexpr import FourierDomainCurrent as If
from .fexpr import FourierDomainVoltage as Vf
from .fexpr import FourierDomainAdmittance as Yf
from .fexpr import FourierDomainImpedance as Zf

from .sexpr import LaplaceDomainExpression as sExpr
from .sexpr import LaplaceDomainTransferFunction as Hs
from .sexpr import LaplaceDomainCurrent as Is
from .sexpr import LaplaceDomainVoltage as Vs
from .sexpr import LaplaceDomainAdmittance as Ys
from .sexpr import LaplaceDomainImpedance as Zs

from .texpr import TimeDomainExpression as tExpr
from .texpr import TimeDomainImpulseResponse as Ht
from .texpr import TimeDomainCurrent as It
from .texpr import TimeDomainVoltage as Vt
from .texpr import TimeDomainAdmittance as Yt
from .texpr import TimeDomainImpedance as Zt

from .omegaexpr import AngularFourierDomainExpression as omegaExpr
from .omegaexpr import AngularFourierDomainTransferFunction as Homega
from .omegaexpr import AngularFourierDomainCurrent as Iomega
from .omegaexpr import AngularFourierDomainVoltage as Vomega
from .omegaexpr import AngularFourierDomainAdmittance as Yomega
from .omegaexpr import AngularFourierDomainImpedance as Zomega

from .noiseexpr import NoiseExpression as noiseExpr
