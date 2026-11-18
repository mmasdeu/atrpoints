r"""
Tests for :mod:`atrpoints.utility`.
"""
import pytest
from sage.all import PolynomialRing, QuadraticField, QQ

from atrpoints.utility import module_generators


@pytest.mark.parametrize('D, alpha', [(29, (3, -13)), (109, (7, -73)), (5, (-2, 2)), (37, (-2, -4)), (41, (2, -14))])
def test_module_generators(D, alpha):
    r"""
    [1, w] must be an O_F-basis of O_K: w integral and disc(minpoly of w over F)
    generating the relative discriminant. ``alpha`` = (s, t) means (s*sqrt(D) + t)/2.
    (D = 29 and D = 109 are the ATR fields of ``examples/atr_example.sage``; the old
    PARI-based version failed for both.)
    """
    F = QuadraticField(D, names='r')
    r = F.gen()
    a = (alpha[0] * r + alpha[1]) / QQ(2)
    K = F.extension(PolynomialRing(F, names='y').gen() ** 2 - a, 'beta')
    one, w = module_generators(K)
    assert one == 1 and w.is_integral()
    assert F.ideal(w.minpoly().discriminant()) == K.relative_discriminant()
