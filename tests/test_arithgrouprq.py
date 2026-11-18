r"""
Tests for :class:`atrpoints.arithgrouprq.ArithGroupRQ` (SL_2(O_F) and Gamma_1(N)).

Run with:  sage -python -m pytest tests
See IMPLEMENTATION.md for the checks and why they matter.
"""
from collections import defaultdict

import pytest
from sage.all import QuadraticField, ZZ, libgap, lcm, matrix, prod, set_random_seed, sign

from atrpoints.arithgrouprq import ArithGroupRQ
from darmonpoints.util import get_weight_vector, multiply_out

F = QuadraticField(5, names='r')
r = F.gen()
eps = (1 + r) / 2
LEVELS = {'2': F.ideal(2), 'sqrt5': F.ideal(r), '3': F.ideal(3), 'p11': F.ideal(3 * r / 2 - ZZ(1) / 2)}
INDICES = {'2': 15, 'sqrt5': 24, '3': 80, 'p11': 120}


@pytest.fixture(scope='module')
def groups():
    return {k: ArithGroupRQ(F, level=N, use_shapiro=False) for k, N in LEVELS.items()}


def random_ambient_word(length=12):
    return [ZZ.random_element(1, 5) * (-1) ** ZZ.random_element(2) for _ in range(length)]


def random_gamma1_word(G, length=12):
    r"""
    A word in the generators of SL_2(O_F) representing an element of Gamma_1(N):
    a random word followed by the inverse of the transversal word of its coset.
    """
    wd = random_ambient_word(length)
    k = G._coset_of(multiply_out(wd, G._ambient_Ugens, G.B.one()))
    return wd + [-i for i in reversed(G._transversal_words[k])]


def random_gamma1(G, length=12):
    return multiply_out(random_gamma1_word(G, length), G._ambient_Ugens, G.B.one())


def test_level_one():
    G = ArithGroupRQ(F, use_shapiro=False)
    assert len(G.gens()) == 4
    gamma4 = matrix(F, 2, 2, [165 * eps + 79, 5 * eps - 2, -8512 * eps - 5160, -159 * eps - 76])
    assert multiply_out(G(gamma4).word_rep, G.Ugens, G.B.one()) == gamma4
    # SL_2(O_F) is perfect once the presentation is complete (see IMPLEMENTATION.md)
    assert tuple(G.abelianization().abelian_invariants()) == ()


@pytest.mark.parametrize('key', list(LEVELS))
def test_presentation(groups, key):
    G = groups[key]
    assert len(G._cosets) == G.index() == INDICES[key]
    assert all(G._is_in_group(g) for g in G.Ugens)
    assert all(multiply_out(rel, G.Ugens, G.B.one()) == 1 for rel in G.get_relation_words())


@pytest.mark.parametrize('key', list(LEVELS))
def test_word_problem(groups, key):
    G = groups[key]
    set_random_seed(1234)
    for _ in range(20):
        g = random_gamma1(G)
        assert G._is_in_group(g)
        wd = G.get_word_rep(g)  # check=True multiplies the word out
        assert multiply_out(wd, G.Ugens, G.B.one()) == g
    with pytest.raises(ValueError):
        G.get_word_rep(matrix(F, 2, 2, [1, 0, 1, 1]))


@pytest.mark.parametrize('key', list(LEVELS))
def test_word_from_ambient(groups, key):
    G = groups[key]
    set_random_seed(4321)
    for _ in range(10):
        wd = random_gamma1_word(G)
        g = multiply_out(wd, G._ambient_Ugens, G.B.one())
        assert multiply_out(G.word_from_ambient(wd), G.Ugens, G.B.one()) == g
        assert G.element_from_ambient_word(wd).quaternion_rep == g
    with pytest.raises(ValueError):
        G.word_from_ambient([4])  # S is not in Gamma_1(N)


@pytest.mark.parametrize('key', ['2', 'sqrt5', '3'])
def test_index_by_coset_enumeration(groups, key):
    r"""
    Independent check (Todd-Coxeter in GAP): the generators, written as words in the
    ambient generators by the continued fraction solver, generate a subgroup of the right
    index of the group presented by the ambient relations. This fails if the ambient
    presentation of SL_2(O_F) is incomplete.
    """
    G = groups[key]
    Fr = libgap.FreeGroup(4)
    fam = libgap.FamilyObj(Fr.One())
    Gfp = Fr / [libgap.AssocWordByLetterRep(fam, wd) for wd in G._ambient_relation_words]
    Ggens = Gfp.GeneratorsOfGroup()

    def amb(g):
        return prod((Ggens[abs(i) - 1] ** sign(i) for i in G._get_word_rep_ambient(g)), Gfp.One())

    assert ZZ(Gfp.Index(Gfp.Subgroup([amb(g) for g in G.Ugens]))) == INDICES[key]


@pytest.mark.parametrize('key', list(LEVELS))
def test_abelianization(groups, key):
    G = groups[key]
    # Tietze simplification does not change the abelianization
    Fr = libgap.FreeGroup(len(G._schreier_matrices))
    fam = libgap.FamilyObj(Fr.One())
    Qraw = Fr / [libgap.AssocWordByLetterRep(fam, wd) for wd in G._schreier_relation_words]
    raw = sorted(ZZ(o) for o in Qraw.AbelianInvariants())
    assert sorted(ZZ(o) for o in G.abelianization().abelian_invariants()) == raw
    assert sorted(G.abelian_invariants_fast()) == raw


def test_abelian_order(groups):
    G = groups['2']
    e = lcm(list(G.abelian_invariants_fast()) + [1])
    set_random_seed(99)
    for _ in range(5):
        wd = G.get_word_rep(random_gamma1(G, 8))
        k = G.abelian_order(wd)
        assert e % k == 0
        wt = get_weight_vector(wd, len(G.gens()))
        assert G._in_relation_lattice([k * o for o in wt])
        assert all(not G._in_relation_lattice([d * o for o in wt]) for d in range(1, k))


def test_find_bounding_cycle(groups):
    G = groups['2']
    n = lcm(list(G.abelian_invariants_fast()) + [1])
    set_random_seed(7)
    g = random_gamma1(G, 8)
    bc = G.find_bounding_cycle(G(g), npow=n, verbose_progress=False)
    # sum of the boundaries d(g1|g2) = g1 - g1 g2 + g2, plus g^n, must vanish
    ans = defaultdict(ZZ)
    for i, g1list, g2 in bc:
        for g1 in g1list:
            ans[tuple(G(g1).quaternion_rep.list())] += i
            ans[tuple((G(g1) * G(g2)).quaternion_rep.list())] -= i
            ans[tuple(G(g2).quaternion_rep.list())] += i
    ans[tuple((g ** n).list())] -= 1
    assert all(v == 0 for v in ans.values())


def test_save_and_load(groups, tmp_path):
    G = groups['sqrt5']
    fname = str(tmp_path / 'gamma1_sqrt5.sobj')
    G.save_data(fname)
    H = ArithGroupRQ(F, level=LEVELS['sqrt5'], data=fname, use_shapiro=False)
    assert H.Ugens == G.Ugens and H.get_relation_words() == G.get_relation_words()
    set_random_seed(5)
    g = random_gamma1(H)
    assert multiply_out(H.get_word_rep(g), H.Ugens, H.B.one()) == g
    with pytest.raises(ValueError):
        ArithGroupRQ(F, level=LEVELS['2'], data=fname, use_shapiro=False)
