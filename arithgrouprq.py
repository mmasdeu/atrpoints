from darmonpoints.arithgroup import *
from darmonpoints.homology import *
from darmonpoints.util import *
from quadcontfrac import *
from sage.all import Matrix, MatrixSpace, identity_matrix, RR, ZZ, QQ, libgap, prod, vector
from collections import deque
import sys


class _LineFile:
    r"""File wrapper for tqdm when stdout is not a terminal (log files): one line per update."""
    def __init__(self, f):
        self.f = f
    def write(self, s):
        s = s.replace('\r', '').strip()
        if s:
            self.f.write(s + '\n')
    def flush(self):
        self.f.flush()


def progress_bar(iterable=None, total=None, desc='', disable=False, **kwds):
    r"""
    A tqdm progress bar. On a terminal it redraws in place; when stdout is redirected to a file
    it prints one line every 15 seconds (and at the end), so that ``tail -f`` shows the progress.
    """
    from tqdm import tqdm
    tty = sys.stdout.isatty()
    return tqdm(iterable, total=total, desc=desc, disable=disable, leave=True,
                file=sys.stdout if tty else _LineFile(sys.stdout),
                mininterval=0.5 if tty else 15, ascii=not tty,
                dynamic_ncols=tty, ncols=None if tty else 100, **kwds)


def multiply_syllables(word, genlist, z, desc=None):
    r"""
    Same as ``multiply_out(word, genlist, z)``, but one power g^a per syllable instead of one
    product per letter (words like T^100000 appear when solving the word problem at level N).
    """
    syl = tietze_to_syllables(word)
    it = progress_bar(syl, desc=desc, unit='syl') if (desc is not None and len(syl) > 2000) else syl
    ans = z
    for i, a in it:
        ans = ans * genlist[i] ** a
    return ans


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
    - ``data`` -- (default None) a file written by :meth:`save_data`; the cosets, the
      Schreier presentation and the Hermite form of the relation matrix are then loaded
      from it instead of being recomputed.

    See IMPLEMENTATION.md for details.
    """
    _DATA_ATTRS = ['_cosets', '_coset_index', '_transversal_words', '_transversal_matrices',
                   '_tree_edges', '_coset_table', '_coset_table_inv', '_schreier_index',
                   '_schreier_matrices', '_schreier_relation_words', 'Ugens', '_relation_words',
                   '_schreier_to_gens', '_rel_hnf']

    def __init__(self, F, level=1, simplify=True, data=None, **kwargs):
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

        if data is not None:
            self._load_data(data)
        elif self.level == 1:
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
        one = self.B.one()
        rels = self._relation_words
        for r in progress_bar(rels, desc='checking relations', unit='rel', disable=len(rels) < 50):
            assert multiply_syllables(r, self.Ugens, one) == 1, r

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
        for k in progress_bar(range(n), desc='Schreier generators', unit='coset'):
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
        for r, k in progress_bar([(r, k) for r in self._ambient_relation_words for k in range(n)],
                                 desc='Reidemeister relators', unit='rel'):
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
            A = multiply_syllables(ans, self.Ugens, self.B.one(), desc='word problem: check')
            assert A == gamma, A**-1 * gamma
        return ans

    def word_from_ambient(self, word, check=True):
        r"""
        Rewrite a word in the generators of SL_2(O_F) (T_1, T_w, U, S, Tietze form) that represents an
        element of Gamma_1(N) as a word in the generators of this group (Reidemeister rewriting).

        Use it when a short ambient word is known (e.g. gamma = gamma_1^k with gamma_1 a short word):
        ``get_word_rep`` goes through the quadratic continued fraction of gamma, which for matrices
        with large entries can give enormous exponents (T_1^100000, ...).
        """
        if self.level == 1:
            ans = free_reduce(list(word))
        else:
            swd, kend = self._rewrite_schreier(word)
            if kend != 0:
                raise ValueError('the word does not represent an element of Gamma_1(N)')
            ans = []
            s2g = self._schreier_to_gens
            for j in swd:
                ans.extend(s2g[j - 1] if j > 0 else [-o for o in reversed(s2g[-j - 1])])
            ans = free_reduce(ans)
        if check:
            assert multiply_syllables(ans, self.Ugens, self.B.one()) == \
                multiply_syllables(word, self._ambient_Ugens, self.B.one())
        return ans

    def element_from_ambient_word(self, word):
        r"""
        The element of this group given by a word in the generators of SL_2(O_F).
        """
        wd = self.word_from_ambient(word)
        return self.element_class(self, quaternion_rep=multiply_syllables(wd, self.Ugens, self.B.one()),
                                  word_rep=wd, check=False)

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
            A = multiply_syllables(ans, gens, identity_matrix(self.F, 2), desc='word problem: ambient check')
            assert A == gamma, A**-1 * gamma
        return ans
    ##############################################
    #  Bounding cycles (verbose version)
    ##############################################

    def find_bounding_cycle(self, gamma, npow=1, verbose_progress=True):
        r"""
        Same output as ``ArithGroup_generic.find_bounding_cycle`` (writes gamma^npow as a sum of
        elements g and boundaries del(g|g')), but

        - prints the progress (word lengths, relations used, terms produced, timings) if
          ``verbose_progress`` is True;
        - multiplies matrices directly (suffix products computed once, right to left) instead of
          building group elements, which at level N != 1 would solve the word problem for every term.
        """
        import time
        from sage.functions.generalized import sgn
        t0 = time.time()
        say = (lambda *a: print('[find_bounding_cycle %6.1fs]' % (time.time() - t0), *a, flush=True)) \
            if verbose_progress else (lambda *a: None)
        assert npow > 0, "npow should be positive"
        one = self.B(1)
        gprime = self(gamma)
        gprimeq = gprime.quaternion_rep
        gens_q = [g.quaternion_rep for g in self.gens()]
        say('word of gamma: length %s, %s syllables (%s generators)'
            % (len(gprime.word_rep), len(tietze_to_syllables(gprime.word_rep)), len(gens_q)))
        rels = self._calculate_relation(
            get_weight_vector(npow * gprime.word_rep, ngens=len(gens_q)), separated=True)
        say('relations needed: %s, total length %s, multiplicities %s'
            % (len(rels), sum(len(r) for _, r in rels), [m for m, _ in rels]))

        def power_chain(g, a, coef):
            # A 2-chain C(a) with boundary a[g] - [g^a], with O(log a) terms (binary splitting):
            #   C(2k) = 2 C(k) + (g^k | g^k),   C(2k+1) = C(2k) + (g^2k | g).
            # It replaces sum_{j=1}^{a-1} (g^j | g) (same boundary, a terms); the two chains differ by a
            # 2-cycle of the cyclic group <g>, i.e. only the choice of bounding cycle changes.
            out = []
            while a > 1:
                if a % 2:
                    out.append((coef, [g ** (a - 1)], g))
                    a -= 1
                else:
                    h = g ** (a // 2)
                    out.append((coef, [h], h))
                    coef *= 2
                    a //= 2
            return out

        def word_terms(syl, mult, bar):
            # terms for the word with syllables syl, counted with multiplicity mult
            n = len(syl)
            suffix = [None] * (n + 1)
            suffix[n] = one
            for k in range(n - 1, -1, -1):
                i, a = syl[k]
                suffix[k] = gens_q[i] ** a * suffix[k + 1]
            out = []
            for jj, (i, a) in enumerate(syl, start=1):
                gq = gens_q[i]
                gaq = gq ** a
                if jj < n:
                    out.append((-mult, [gaq], suffix[jj]))
                if a < 0:
                    out.append((mult, [gaq ** -1], gaq))
                if abs(a) > 1:
                    out.extend(power_chain(gq, abs(a), -sgn(a) * mult))
                bar.update(1)
                bar.set_postfix(terms=nterms[0] + len(out), refresh=False)
            nterms[0] += len(out)
            return out

        words = [(tietze_to_syllables(gprime.word_rep), npow)]
        for m, rel in rels:
            assert m > 0
            words.append((tietze_to_syllables(rel), m))
        nterms = [0]
        ans = power_chain(gprimeq, npow, -1)
        with progress_bar(total=sum(len(w) for w, _ in words), desc='bounding cycle', unit='syl',
                          disable=not verbose_progress) as bar:
            for syl, m in words:
                ans += word_terms(syl, m, bar)
        # fix the boundaries del(1|1) that were miscounted
        check = 0 if gprimeq != 1 else -npow
        for i, g1list, g2 in progress_bar(ans, desc='del(1|1) check', unit='term',
                                          disable=not verbose_progress or len(ans) < 2000):
            for g1 in g1list:
                if g1 == 1:
                    check += i
                if g1 * g2 == 1:
                    check -= i
                if g2 == 1:
                    check += i
        ans.append((-check, [one], one))
        say('done: %s terms (%s pairs (A, B))' % (len(ans), sum(len(o[1]) for o in ans)))
        return ans

    ##############################################
    #  Relations (fast lift via a cached Hermite form)
    ##############################################

    def _relation_hnf(self):
        r"""
        Column Hermite form of M^T (M = relation matrix, rows = weight vectors of the relation words),
        computed with PARI's ``mathnf(M^T, 1)``: M^T * U = [0 | H] with U unimodular and H upper
        triangular of full column rank. Returns (H, U2) with U2 the last columns of U, so that
        M^T * U2 = H. Computed once and cached (and saved by :meth:`save_data`).
        """
        if getattr(self, '_rel_hnf', None) is None:
            import time
            from sage.libs.pari import pari
            t0 = time.time()
            M = Matrix(ZZ, self.get_relation_matrix())
            Hp, Up = pari(M.transpose()).mathnf(1)
            H, U = Matrix(ZZ, Hp), Matrix(ZZ, Up)
            U2 = U[:, U.ncols() - H.ncols():]
            assert M.transpose() * U2 == H
            self._rel_hnf = (H, U2)
            print('[relations] Hermite form of the %s x %s relation matrix (rank %s) in %.1fs'
                  % (M.nrows(), M.ncols(), H.ncols(), time.time() - t0), flush=True)
        return self._rel_hnf

    def _solve_relations(self, wt):
        r"""
        Integer vector x with x * M = wt (M = relation matrix), or None if there is none.
        """
        H, U2 = self._relation_hnf()
        rem = vector(ZZ, wt)
        k = H.ncols()
        y = [ZZ(0)] * k
        for j in range(k - 1, -1, -1):
            col = H.column(j)
            p = max(i for i in range(H.nrows()) if col[i] != 0)   # pivot row of column j
            q, r = rem[p].quo_rem(col[p])
            if r != 0:
                return None
            y[j] = q
            if q != 0:
                rem -= q * col
        if rem != 0:
            return None
        return U2 * vector(ZZ, y)

    def _solve_relations_sparse(self, wt):
        r"""
        Integer vector x with x * M = wt (M = relation matrix), with small entries, or None.

        The system M^T x = wt (one equation per generator, one unknown per relation) is very sparse.
        It is solved by elimination with pivots of coefficient +-1, chosen by the Markowitz rule
        (fewest fill-in), with the unknowns that are never pivots set to 0: this is the abelian
        version of Tietze transformations and keeps the coefficients small. The equations left
        without a unit pivot (the torsion part, here a few equations) are solved with PARI's
        ``matsolvemod``. A Hermite-form solution (:meth:`_solve_relations`) is valid as well but
        its entries can have hundreds of digits, which makes the bounding cycle useless numerically.
        """
        from sage.libs.pari import pari
        M = self.get_relation_matrix()
        nrel, ngen = M.nrows(), M.ncols()
        eqs = [dict() for _ in range(ngen)]
        for (i, c), v in M.dict().items():
            eqs[c][i] = ZZ(v)
        rhs = [ZZ(o) for o in wt]
        occ = [set() for _ in range(nrel)]        # unknown -> equations containing it
        for c, e in enumerate(eqs):
            for i in e:
                occ[i].add(c)
        active = set(range(ngen))
        pivots = []
        while True:
            best = None
            for c in active:
                e = eqs[c]
                if not e:
                    continue
                le = len(e) - 1
                for i, v in e.items():
                    if v == 1 or v == -1:
                        cost = le * (len(occ[i]) - 1)
                        if best is None or cost < best[0]:
                            best = (cost, c, i)
                            if cost == 0:
                                break
                if best is not None and best[0] == 0:
                    break
            if best is None:
                break
            _, c, i = best
            e, ci = eqs[c], eqs[c][i]
            active.discard(c)
            for c2 in list(occ[i]):
                if c2 == c:
                    continue
                e2 = eqs[c2]
                f = e2[i] * ci                     # ci = +-1, so e2[i] / ci = e2[i] * ci
                for j, v in e.items():
                    w = e2.get(j, 0) - f * v
                    if w == 0:
                        if j in e2:
                            del e2[j]
                            occ[j].discard(c2)
                    else:
                        if j not in e2:
                            occ[j].add(c2)
                        e2[j] = w
                rhs[c2] -= f * rhs[c]
            pivots.append((c, i))
        x = {}
        rest = [c for c in active if eqs[c] or rhs[c] != 0]
        if rest:
            if any(not eqs[c] for c in rest):
                return None                        # an equation 0 = rhs != 0
            unk = sorted(set(i for c in rest for i in eqs[c]))
            pos = {i: k for k, i in enumerate(unk)}
            A = Matrix(ZZ, len(rest), len(unk))
            for r_, c in enumerate(rest):
                for i, v in eqs[c].items():
                    A[r_, pos[i]] = v
            sol = pari(A).matsolvemod(0, pari(vector(ZZ, [rhs[c] for c in rest])).Col())
            if sol == 0:
                return None
            for k, i in enumerate(unk):
                if ZZ(sol[k]) != 0:
                    x[i] = ZZ(sol[k])
        for c, i in reversed(pivots):
            e = eqs[c]
            val = (rhs[c] - sum(v * x.get(j, 0) for j, v in e.items() if j != i)) * e[i]
            if val != 0:
                x[i] = val
        ans = vector(ZZ, nrel)
        for i, v in x.items():
            ans[i] = v
        if ans * M != vector(ZZ, wt):
            return None
        return ans

    def _calculate_relation(self, wt, separated=False):
        r"""
        Same as ``ArithGroup_generic._calculate_relation``: write the weight vector ``wt`` as an
        integer combination of the weight vectors of the relation words and return the
        corresponding relation words. Uses the cached Hermite form instead of a new lift each time.
        """
        relwords = self.get_relation_words()
        if len(relwords) == 0:
            return []
        lin = self._solve_relations_sparse(wt)
        if lin is None:
            raise ValueError('weight vector is not a combination of the relations')
        assert lin * Matrix(ZZ, self.get_relation_matrix()) == vector(ZZ, wt)
        ans = []
        for i, lam in enumerate(lin):
            if lam == 0:
                continue
            relation = relwords[i]
            if separated:
                ans.append((ZZ(-lam), relation) if lam < 0 else (ZZ(lam), [-j for j in reversed(relation)]))
            else:
                ans += ZZ(-lam) * relation if lam < 0 else ZZ(lam) * [-j for j in reversed(relation)]
        return ans

    ##############################################
    #  Saving / loading the computed data
    ##############################################

    def save_data(self, fname):
        r"""
        Save the cosets, the presentation and (if computed) the Hermite form of the relation
        matrix to ``fname``, so that ``ArithGroupRQ(F, level, data=fname)`` skips recomputing them.
        """
        from sage.misc.persist import save
        d = {a: getattr(self, a) for a in self._DATA_ATTRS if hasattr(self, a)}
        d['_level_gens'] = self.level.gens_reduced()
        d['_field_poly'] = self.F.defining_polynomial()
        save(d, fname)

    def _load_data(self, fname):
        from sage.misc.persist import load
        d = load(fname) if isinstance(fname, str) else fname
        if d['_field_poly'] != self.F.defining_polynomial() or self.F.ideal(d['_level_gens']) != self.level:
            raise ValueError('data file %s is for a different field or level' % fname)
        for a in self._DATA_ATTRS:
            if a in d:
                setattr(self, a, d[a])
        print('[ArithGroupRQ] loaded %s (%s generators, %s relations)'
              % (fname, len(self.Ugens), len(self._relation_words)), flush=True)

    ##############################################
    #  Abelianization (fast, via Smith / Hermite forms)
    ##############################################

    def abelian_invariants_fast(self):
        r"""
        Invariants of the abelianization (0 = a factor Z), from the elementary divisors of the
        relation matrix. Same as ``self.abelianization().abelian_invariants()``, but without
        building darmonpoints' FGP-module (which takes minutes for hundreds of generators).
        """
        M = Matrix(ZZ, self.get_relation_matrix())
        # (elementary_divisors pads with zeros up to the number of rows, so keep only the nonzero ones)
        ed = [d for d in M.elementary_divisors() if d not in (0, 1)]
        ed += [ZZ(0)] * (M.ncols() - M.rank())
        return tuple(ed)

    def _in_relation_lattice(self, wt):
        return self._solve_relations_sparse(wt) is not None

    def abelian_order(self, x):
        r"""
        Order of the image of ``x`` (an element, a matrix or a word) in the abelianization
        (``Infinity`` if it has infinite order).
        """
        from sage.rings.infinity import Infinity
        from sage.arith.misc import divisors
        from sage.arith.functions import lcm
        if not isinstance(x, list):
            x = x.word_rep if hasattr(x, 'word_rep') else self(x).word_rep
        wt = vector(ZZ, get_weight_vector(x, len(self.gens())))
        inv = self.abelian_invariants_fast()
        e = lcm([ZZ(1)] + [d for d in inv if d != 0])
        for k in divisors(e):
            if self._in_relation_lattice(k * wt):
                return k
        return Infinity
