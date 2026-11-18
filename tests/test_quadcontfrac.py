r"""
Tests for :mod:`atrpoints.quadcontfrac` (continued fractions in real quadratic fields).
"""
import pytest
from sage.all import QuadraticField, ZZ

from atrpoints.quadcontfrac import all_qcf_up_to_length, quadratic_continued_fraction


def evaluate_cf(quotients):
    r"""
    The value a_0 - 1/(a_1 - 1/(a_2 - ...)) or a_0 + 1/(a_1 + ...) is not fixed by the
    convention, so recompute it from the convergents the way ``quotients_and_limits`` does.
    """
    p, q = [ZZ(0), ZZ(1)], [ZZ(1), ZZ(0)]
    for a in quotients:
        p.append(p[-1] * a + p[-2])
        q.append(q[-1] * a + q[-2])
    return p[-1] / q[-1]


@pytest.mark.parametrize('D', [5, 13, 29])
def test_quadratic_continued_fraction(D):
    F = QuadraticField(D, names='r')
    r = F.gen()
    w = F.ring_of_integers().ring_generators()[0]
    for x in [F(7) / 3, (2 + w) / 5, (11 * (17 + 8 * w) + 5) / 8, -r / 9]:
        cf = quadratic_continued_fraction(x)
        assert all(a.is_integral() for a in cf)
        assert evaluate_cf(cf) == x


def test_all_qcf_up_to_length():
    F = QuadraticField(29, names='r')
    w = F.ring_of_integers().ring_generators()[0]
    x = (2 + w) / 2
    # get_img_part scores a candidate expansion (the ATR code uses the imaginary part of the
    # limits); a constant score accepts every candidate above the threshold.
    cfs = all_qcf_up_to_length(x, max_length=4, get_img_part=lambda vec: 1, threshold=0)
    assert len(cfs) > 0
    assert all(evaluate_cf(cf) == x and len(cf) <= 4 for cf in cfs)
