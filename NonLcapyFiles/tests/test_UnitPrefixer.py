import lcapy
from sympy import latex
from lcapy.unitPrefixer import SIUnitPrefixer
from lcapy import omega0
from lcapy import resistance

lcapy.state.show_units = True


# I know this is a terrible test, but it does for now
def test_UnitPrefixer():

    prefixer = SIUnitPrefixer()
    assert latex(prefixer.getSIPrefixedMul(10)) == "10.0"
    assert latex(prefixer.getSIPrefixedMul(100)) == "100.0"
    assert latex(prefixer.getSIPrefixedMul(1000)) == "1.0 \\text{k}"
    assert latex(prefixer.getSIPrefixedMul(100.0)) == "100.0"
    assert latex(prefixer.getSIPrefixedMul(1000.0)) == "1.0 \\text{k}"
    assert latex(prefixer.getSIPrefixedMul(100 * omega0)) == "100.0 \\omega_{0}"
    assert latex(prefixer.getSIPrefixedMul(1000 * omega0)) == "1.0 \\omega_{0} \\text{k}"
    assert latex(prefixer.getSIPrefixedMul(1000000 * omega0)) == "1.0 \\omega_{0} \\text{M}"
    assert latex(prefixer.getSIPrefixedMul(100000 * omega0)) == "0.1 \\omega_{0} \\text{M}"
    assert latex(prefixer.getSIPrefixedMul(100.0 * omega0)) == "100.0 \\omega_{0}"
    assert latex(prefixer.getSIPrefixedMul(1000.0 * omega0)) == "1.0 \\omega_{0} \\text{k}"
    assert latex(prefixer.getSIPrefixedMul(resistance(1000))) == "1.0 \\text{k} \\Omega"
