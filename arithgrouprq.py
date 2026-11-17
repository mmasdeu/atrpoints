from darmonpoints.arithgroup import *
from darmonpoints.homology import *
from darmonpoints.util import *
from quadcontfrac import *
from sage.all import Matrix, MatrixSpace, identity_matrix, RR, ZZ, QQ, libgap, prod
from collections import deque


def free_reduce(word):
    r"""
    Freely reduce a word in Tietze form (cancel adjacent ``i, -i``), in linear time.
    """
    ans = []
    for i in word:
        if ans and ans[-1] == -i:
            ans.pop()
        else:
            ans.append(i)
    return ans


class ArithGroupRQ(ArithGroup_matrix_generic):
    r"""
    The group Gamma_1(N) inside SL_2(O_F), for F a real quadratic field.

    INPUT:

    - ``F`` -- a real quadratic field
    - ``level`` -- an integral ideal N of F (or a generator of it), default 1.
      For N = 1 the group is SL_2(O_F) with the presentation of [CD].
    - ``simplify`` -- (default True) whether to simplify the Reidemeister-Schreier
      presentation of Gamma_1(N) with GAP's Tietze transformations.

    For N != 1 the presentation and the word problem are obtained from those of
    SL_2(O_F) by the Reidemeister-Schreier method, using that the right cosets
    Gamma_1(N) \ SL_2(O_F) are parametrized by bottom rows mod N.
    See IMPLEMENTATION.md for details.
    """
    def __init__(self, F, level=1, simplify=True, **kwargs):
        self.F = F
        self.B = MatrixSpace(F, 2)
        self.level = F.ideal(level)
        if self.level == 0 or not self.level.is_integral():
            raise ValueError('The level must be a nonzero integral ideal')
        self.w = F.ring_of_integers().ring_generators()[0]
        self.phi =  self.w.coordinates_in_terms_of_powers()
        (self.eps,) = F.units()
        self._ambient_Ugens = [ # Notation as in [CD]
            Matrix(self.F,2,2,[1,1,0,1]), # g_1 = T_1
            Matrix(self.F,2,2,[1,self.w,0,1]), # 2 = T_w
            Matrix(self.F,2,2,[self.eps,0,0,self.eps**-1]), # 3 = U
            Matrix(self.F,2,2,[0,-1,1,0])] # 4 = S
        for g in self._ambient_Ugens:
            g.set_immutable()
        self._ambient_relation_words = [
            [4, 3, -4, 3], [4,4,4,4], [4,1,4,1,4,1,4,4], [1,2,-1,-2]] + [[i,4,4,-i,-4,-4] for i in range(1,4)]
        wd = lambda n, i : n * [i] if n > 0 else (-n) * [-i]
        Tword = lambda x : wd(self.phi(x)[0], 1) + wd(self.phi(x)[1], 2) # T_x = T_1^a T_w^b
        for theta in [self.F(1), self.w]:
            a0, b0 = self.phi(theta)
            a1, b1 = self.phi(-theta*self.eps**2)
            self._ambient_relation_words.append([3] + wd(a0, 1) + wd(b0,2) + [-3] + wd(a1, 1) + wd(b1, 2))
        # Unit analogue of (ST)^3 = S^2 (Cohn's universal GE_2 relation for the unit eps):
        # T_eps S T_{eps^-1} S T_eps S = S^2 U. Without it the group presented is
        # a proper cover of SL_2(O_F) (detected by coset enumeration, see IMPLEMENTATION.md).
        self._ambient_relation_words.append(
            Tword(self.eps) + [4] + Tword(self.eps**-1) + [4] + Tword(self.eps) + [4, -3, 4, 4])
        for r in self._ambient_relation_words:
            assert multiply_out(r, self._ambient_Ugens, self.B.one()) == 1, r

        if self.level == 1:
            self.Ugens = self._ambient_Ugens
            self._relation_words = self._ambient_relation_words
        else:
            self._compute_coset_data()
            self._compute_schreier_presentation()
            self._compute_simplified_presentation(simplify)
        self._gens = [
                self.element_class(
                    self, quaternion_rep=newg, word_rep=[i + 1], check=False
                ) for i, newg in enumerate(self.Ugens)
        ]
        super().__init__(**kwargs)
        ArithGroup_generic.__init__(self, **kwargs)
        for g in self.Ugens:
            assert self._is_in_group(g), g
        for r in self._relation_words:
            assert self(r).quaternion_rep == 1, r

    def _repr_(self):
        if self.level == 1:
            return 'SL_2(O_F) for F = %s' % self.F
        return 'Gamma_1(%s) in SL_2(O_F) for F = %s' % (self.level.gens_reduced(), self.F)

    def embed(self, x, prec=None):
        return x

    ##############################################
    #  Membership and cosets
    ##############################################

    def _is_in_group(self, gamma):
        r"""
        Whether the matrix ``gamma`` belongs to Gamma_1(N).
        """
        gamma = self.B(gamma)
        if gamma.determinant() != 1 or not all(x.is_integral() for x in gamma.list()):
            return False
        N = self.level
        return gamma[1,0] in N and gamma[1,1] - 1 in N

    def index(self):
        r"""
        The index [SL_2(O_F) : Gamma_1(N)] = Nm(N)^2 prod_{P | N} (1 - Nm(P)^-2).
        """
        N = self.level
        return ZZ(N.norm()**2 * prod(1 - QQ(P.norm())**-2 for P, _ in N.factor()))

    def _coset_key(self, c, d):
        N = self.level
        return (N.reduce(c), N.reduce(d))

    def _act_on_coset(self, key, g):
        r"""
        Right action of the matrix g on the coset Gamma_1(N) h, represented
        by the bottom row (c, d) of h mod N.
        """
        c, d = key
        return self._coset_key(c * g[0,0] + d * g[1,0], c * g[0,1] + d * g[1,1])

    def _coset_of(self, gamma):
        return self._coset_index[self._coset_key(gamma[1,0], gamma[1,1])]

    def transversal(self):
        r"""
        Matrices r_k, one for each right coset of Gamma_1(N) in SL_2(O_F),
        such that the k-th coset is Gamma_1(N) r_k (r_0 = 1).
        """
        if self.level == 1:
            return [self.B.one()]
        return self._transversal_matrices

    def _compute_coset_data(self):
        r"""
        Enumerate Gamma_1(N) \ SL_2(O_F) by BFS from the trivial coset, and build:

        - ``_cosets``, ``_coset_index``: keys (bottom rows mod N) and their indices;
        - ``_transversal_words``, ``_transversal_matrices``: Schreier transversal;
        - ``_coset_table[k][i]`` (resp. ``_coset_table_inv``): coset of k * g_i (resp. g_i^-1);
        - ``_tree_edges``: pairs (k, i) whose Schreier generator is trivial.
        """
        gens = self._ambient_Ugens
        gens_inv = [g**-1 for g in gens]
        ngens = len(gens)
        v0 = self._coset_key(0, 1)
        self._cosets = [v0]
        self._coset_index = {v0: 0}
        self._transversal_words = [[]]
        self._transversal_matrices = [self.B.one()]
        self._tree_edges = set()
        table = [[None] * ngens]
        queue = deque([0])
        while queue:
            k = queue.popleft()
            v = self._cosets[k]
            for sgn in [1, -1]:
                for i in range(ngens):
                    g = gens[i] if sgn == 1 else gens_inv[i]
                    w = self._act_on_coset(v, g)
                    kk = self._coset_index.get(w)
                    if kk is None:
                        kk = len(self._cosets)
                        self._cosets.append(w)
                        self._coset_index[w] = kk
                        self._transversal_words.append(self._transversal_words[k] + [sgn * (i + 1)])
                        self._transversal_matrices.append(self._transversal_matrices[k] * g)
                        table.append([None] * ngens)
                        queue.append(kk)
                        # The Schreier generator t_k g t_kk^-1 is freely trivial
                        self._tree_edges.add((k, i) if sgn == 1 else (kk, i))
                    if sgn == 1:
                        table[k][i] = kk
        for r in self._transversal_matrices:
            r.set_immutable()
        n = len(self._cosets)
        assert n == self.index(), 'Found %s cosets, expected %s' % (n, self.index())
        self._coset_table = table
        table_inv = [[None] * ngens for _ in range(n)]
        for k in range(n):
            for i in range(ngens):
                table_inv[table[k][i]][i] = k
        self._coset_table_inv = table_inv
        assert len(self._tree_edges) == n - 1

    ##############################################
    #  Reidemeister-Schreier
    ##############################################

    def _compute_schreier_presentation(self):
        r"""
        Compute the Schreier generators s_{k,i} = r_k g_i r_{k g_i}^-1 for the
        non-tree edges (k, i), and the Reidemeister relators t_k R t_k^-1
        rewritten in them.
        """
        gens = self._ambient_Ugens
        ngens = len(gens)
        n = len(self._cosets)
        self._schreier_index = {}
        self._schreier_matrices = []
        reps = self._transversal_matrices
        for k in range(n):
            for i in range(ngens):
                if (k, i) in self._tree_edges:
                    continue
                kk = self._coset_table[k][i]
                s = reps[k] * gens[i] * reps[kk]**-1
                assert self._is_in_group(s), s
                s.set_immutable()
                self._schreier_matrices.append(s)
                self._schreier_index[(k, i)] = len(self._schreier_matrices)
        assert len(self._schreier_matrices) == n * (ngens - 1) + 1
        rels = []
        for r in self._ambient_relation_words:
            for k in range(n):
                wd, kend = self._rewrite_schreier(r, start=k)
                assert kend == k
                wd = free_reduce(wd)
                if len(wd) > 0:
                    rels.append(wd)
        self._schreier_relation_words = rels

    def _rewrite_schreier(self, word, start=0):
        r"""
        Reidemeister rewriting: given a word in the ambient generators, return
        the corresponding word in the Schreier generators, starting at the coset
        ``start``, together with the final coset. The word represents an element
        of Gamma_1(N) iff the final coset is 0 (when start = 0).
        """
        k = start
        ans = []
        idx = self._schreier_index
        table = self._coset_table
        table_inv = self._coset_table_inv
        for l in word:
            if l > 0:
                i = l - 1
                j = idx.get((k, i))
                if j is not None:
                    ans.append(j)
                k = table[k][i]
            else:
                i = -l - 1
                k = table_inv[k][i]
                j = idx.get((k, i))
                if j is not None:
                    ans.append(-j)
        return ans, k

    def _compute_simplified_presentation(self, simplify=True):
        r"""
        Set ``Ugens``, ``_relation_words`` and the map ``_schreier_to_gens`` that
        writes every Schreier generator as a word in ``Ugens``. If ``simplify``
        is True, use GAP's Tietze transformations on the Schreier presentation.
        """
        m = len(self._schreier_matrices)
        if not simplify:
            self.Ugens = self._schreier_matrices
            self._relation_words = self._schreier_relation_words
            self._schreier_to_gens = [[j + 1] for j in range(m)]
            return
        Fr = libgap.FreeGroup(m)
        fam = libgap.FamilyObj(Fr.One())
        to_gap = lambda wd: libgap.AssocWordByLetterRep(fam, wd)
        from_gap = lambda x: [ZZ(o) for o in x.UnderlyingElement().LetterRepAssocWord().sage()]
        Q = Fr / [to_gap(wd) for wd in self._schreier_relation_words]
        iso = Q.IsomorphismSimplifiedFpGroup()
        H = iso.Range()
        verbose('Simplified presentation: %s generators -> %s generators' % (m, len(H.GeneratorsOfGroup())))
        self.Ugens = []
        for h in H.GeneratorsOfGroup():
            wd = from_gap(iso.PreImagesRepresentative(h))
            g = multiply_out(wd, self._schreier_matrices, self.B.one())
            g.set_immutable()
            self.Ugens.append(g)
        self._relation_words = [from_gap(r) for r in H.RelatorsOfFpGroup()]
        self._schreier_to_gens = [from_gap(iso.ImageElm(q)) for q in Q.GeneratorsOfGroup()]
        for s, wd in zip(self._schreier_matrices, self._schreier_to_gens):
            assert multiply_out(wd, self.Ugens, self.B.one()) == s

    ##############################################
    #  Word problem
    ##############################################

    def get_word_rep(self, gamma, check=True):
        r"""
        Write ``gamma`` (a matrix in Gamma_1(N)) as a word in Tietze form in the
        generators of this group.
        """
        if self.level == 1:
            return self._get_word_rep_ambient(gamma, check=check)
        if not self._is_in_group(gamma):
            raise ValueError('Matrix does not belong to Gamma_1(N):\n%s' % gamma)
        wd = self._get_word_rep_ambient(gamma, check=check)
        swd, kend = self._rewrite_schreier(wd)
        assert kend == 0
        ans = []
        s2g = self._schreier_to_gens
        for j in swd:
            ans.extend(s2g[j - 1] if j > 0 else [-o for o in reversed(s2g[-j - 1])])
        ans = free_reduce(ans)
        if check:
            A = multiply_out(ans, self.Ugens, self.B.one())
            assert A == gamma, A**-1 * gamma
        return ans

    def _get_word_rep_ambient(self, gamma, check=True):
        r"""
        Solve the word problem in SL_2(O_F) with respect to the generators of [CD],
        using quadratic continued fractions.
        """
        Mlist = []
        if gamma[1,0] != 0:
            for i, o in enumerate(quadratic_continued_fraction(gamma[0,0]/gamma[1,0])):
                a, b = self.phi(o)
                Mlist.extend([(0, (-1)**i * a), (1, (-1)**i * b), (3, 1)])
        gens = self._ambient_Ugens
        A = identity_matrix(self.F,2)
        for s, i in Mlist:
            A = A * gens[s]**i
        B = A**-1 * gamma
        assert B[1,0] == 0

        epsi = B[0,0]
        v0 = self.F.real_embeddings()[0]
        if v0(self.eps) < 0:
            v0 = self.F.real_embeddings()[1]
        b0, b1 = self.phi(B[0,1] / epsi)
        if v0(epsi) < 0:
            epsi = -epsi
            Mlist.append((3,2))
        assert v0(epsi) > 0 and v0(self.eps) > 0, (v0(epsi), v0(self.eps))
        j = RR(v0(epsi).log() / v0(self.eps).log())
        j = j.round()
        found = False
        for i in [j,j-1,j+1]:
            if self.eps**i == epsi:
                found = True
                break
        assert found, 'Could not find i such that eps^i = epsi'
        Mlist.extend([(2, i), (0, b0), (1, b1)])

        ans = syllables_to_tietze(Mlist)
        # Test
        if check:
            A = identity_matrix(self.F,2)
            for i in ans:
                A = A * gens[i-1] if i > 0 else A * gens[-i-1]**-1
            assert A == gamma, A**-1 * gamma
        return ans
