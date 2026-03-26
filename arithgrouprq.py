from darmonpoints.arithgroup import *
from darmonpoints.homology import *
from darmonpoints.util import *
from quadcontfrac import *
from sage.all import Matrix, MatrixSpace, identity_matrix, RR, ZZ

class ArithGroupRQ(ArithGroup_matrix_generic):
    def __init__(self, F, **kwargs):
        self.F = F
        self.B = MatrixSpace(F, 2)
        self.w = F.ring_of_integers().ring_generators()[0]
        self.phi =  self.w.coordinates_in_terms_of_powers()
        (self.eps,) = F.units()
        self.Ugens = [ # Notation as in [CD]
            Matrix(self.F,2,2,[1,1,0,1]), # g_1 = T_1
            Matrix(self.F,2,2,[1,self.w,0,1]), # 2 = T_w
            Matrix(self.F,2,2,[self.eps,0,0,self.eps**-1]), # 3 = U
            Matrix(self.F,2,2,[0,-1,1,0])] # 4 = S
        self._gens = [
                self.element_class(
                    self, quaternion_rep=newg, word_rep=[i + 1], check=False
                ) for i, newg in enumerate(self.Ugens)
        ]
        self._relation_words = [
            [4, 3, -4, 3], [4,4,4,4], [4,1,4,1,4,1,4,4], [1,2,-1,-2]] + [[i,4,4,-i,-4,-4] for i in range(1,4)]
        for theta in [self.F(1), self.w]:
            a0, b0 = self.phi(theta)
            a1, b1 = self.phi(-theta*self.eps**2)
            wd = lambda n, i : n * [i] if n > 0 else (-n) * [-i]
            self._relation_words.append([3] + wd(a0, 1) + wd(b0,2) + [-3] + wd(a1, 1) + wd(b1, 2))
        super().__init__(**kwargs)
        ArithGroup_generic.__init__(self, **kwargs)
        for r in self._relation_words:
            assert self(r).quaternion_rep == 1, r

    def embed(self, x, prec=None):
        return x
    
    def get_word_rep(self, gamma, check=True):
        Mlist = []
        if gamma[1,0] != 0:
            for i, o in enumerate(quadratic_continued_fraction(gamma[0,0]/gamma[1,0])):
                a, b = self.phi(o)
                Mlist.extend([(0, (-1)**i * a), (1, (-1)**i * b), (3, 1)])
        gens = self.Ugens
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

