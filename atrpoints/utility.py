r"""
Auxiliary functions for the ATR point computations with level: optimal embeddings
(tau_0 and gamma_tau), decomposition of the limits of integration into elementary
matrices, continued fractions, and the class :class:`Limit`.
"""
import resource

from sage.all import *
from itertools import product  # after sage.all: itertools.product, not Sage's product

from .quadcontfrac import *

prec_bits = 2000
TIME_LIMIT = 10
LOW_PREC_BITS = 300
LOWER_PREC_BITS = 100


def get_memory_usage():
    r"""
    Peak resident memory of this process, in MB (replaces Sage's removed ``get_memory_usage``).
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


class Limit(SageObject):
    def __init__(self,x0,x1,y0,y1,conj=1,scale = 1):
        self.x0,self.x1,self.y0,self.y1=x0,x1,y0,y1
        self.x0app,self.x1app,self.y0app=ComplexField()(x0),ComplexField()(x1),ComplexField()(y0)
        if y1 is Infinity:
            self.y1app=Infinity
        else:
            self.y1app=ComplexField()(y1)
        self.conj=conj
        self.scale = scale

    def volume(self):
        v1 = RR(acosh(1+(self.x0-self.y0).abs()/(2*self.x0.imag()*self.y0.imag())))
        v2 = RR(acosh(1+(self.x1-self.y1).abs()/(2*self.x1.imag()*self.y1.imag())))
        return v1*v2

    def unpack(self):
        return self.x0,self.x1,self.y0,self.y1

    def imaginary_part_approx(self):
        x0app,x1app,y0app,y1app=self.x0app,self.x1app,self.y0app,self.y1app
        rr=min([(z1.imag()*z2.imag()).sqrt() for z1,z2 in [(x0app,x1app),(x0app,y1app),(y0app,x1app),(y0app,y1app)] if z2 is not Infinity])
        return rr
    def __str__(self):
        return str(self.unpack())
    def x(self):
        return (x0,x1)
    def y(self):
        return (y0,y1)
    def sw(self):
        return self.x0
    def se(self):
        return self.x1
    def nw(self):
        return self.y0
    def ne(self):
        return self.y1


def tau0_and_cusp_withlevel(base,ext,v0,N,find_all = False):
    F = base
    ww=module_generators(ext)[1]
    assert ww.is_integral()
    wwminpoly = ww.minpoly()
    # phi = wwminpoly.base_ring().embeddings(F)[0]
    # K = NumberField(ww.minpoly().map_coefficients(phi),names ='w')

    K = NumberField(ww.minpoly(),names ='w')
    eps = find_the_unit_of(F,K)
    print('The unit of K that we are using is eps=',eps)
    w = K.gen(0)
    idealN = F.ideal(N)

    #these are elements in K whose norm over Q is the same as the norm over Q of N. We actually want an element in K whose norm over F coincides (up to units) with N.
    if N == 1 or N == -1:
        candidates=K.elements_of_norm(1)
    else:
        candidates=K.elements_of_norm(N.norm())
    print('Found %s candidates'%len(candidates))

    found = False
    all_B = []
    #repassar aixo perque abans no donava be!!!
    for elt in candidates:
        if elt.norm(F) in idealN:
            #the element of norm N is of the form NN=a+cw
            a,c=elt.vector().list()
            #the columns of this matrix give a base of O_K as O_F-module: O_K=<a+cw,b+dw>; let't call this the basis BB
            try:
                B=matrix_from_left_column(F,a,c)
                found = True
                if find_all == True:
                    all_B.append(B)
                else:
                    break
            except RuntimeError: pass

    if not found:
        print('Error: element in K of norm N not found; are you sure that K has elements of that norm?')
        raise RuntimeError


    #this is the matrix of "multiplication by w" in the basis AA=<1,w>
    M=matrix(F,[[0,-(w.minpoly()[0])],[1,-(w.minpoly()[1])]])
    tau0vec = []
    gtauvec = []
    u_F=F.units()[0]

    for B in all_B:
        #the optimal embedding is given by the matrix of "multiplication by w" in the basis BB; we can get it by the usual formula:
        W=B**-1*M*B

        #print('the optimal embedding sends w to the matrix W given by:')
        #print(W)
        # let's check it's correct print('checking...')
        assert W[1,0] in idealN
        print('W belongs to Gamma_0(N)...')
        assert W.minpoly() == w.minpoly()
        print('and it has the same minimal polynomial as w!')
        #now we find a suitable unit in K...
        e1,e2=eps.vector()
        gamma_tau=e1+e2*W
        print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        print('W=')
        print(W)
        print('gamma_tau=')
        print(gamma_tau)
        print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        assert gamma_tau.minpoly() == eps.minpoly()
        print('which we checked has the same minimal polynomial')
        #it belongs to Gamma_0(N)
        if not is_in_Gamma_0(gamma_tau,N):
            print('apparently you do not get an embedding of level N; try to use the conjugate of alpha instead')

        # gamma_tau lives in Gamma_0(N), but maybe not in Gamma_1(N): replace eps by
        # the smallest power eps^k whose image has top-left entry = +-1 mod N
        # (then -I*gamma_tau^k, which acts in the same way, lies in Gamma_1(N)).
        # Note that the ATR point gets multiplied by k.
        gamma_tau, k = power_in_Gamma_1(gamma_tau, N)
        print('Using eps^%s, whose image lies in Gamma_1(N)' % k)
        assert (gamma_tau[0,0] - 1) in idealN and gamma_tau[1,0] in idealN
        assert gamma_tau.det() == 1

        aa,bb,cc,dd = W.apply_map(v0).list()
        if cc > 0:
            tau0 = ((aa-dd)+sqrt((dd-aa)**2+4*bb*cc))/(2*cc)
            print('tau0=',tau0)
        else:
            tau0 = ((aa-dd)-sqrt((dd-aa)**2+4*bb*cc))/(2*cc)
            print('tau0=', tau0)
        assert tau0.imag() > 0
        # print('tau=',tau0)
        tau0vec.append(tau0)
        gtauvec.append(gamma_tau)
    if find_all == True:
        return tau0vec ,gtauvec
    else:
        return tau0vec[0],gtauvec[0]

def tau0_and_cusp(base, ext, e, v0, embedding, extra_parameters = None):
    if extra_parameters is None:
        extra_parameters = [0,0,0]
    w=module_generators(ext)[1]
    assert w.is_integral()

    eps = find_the_unit_of(base,ext)
    #now we compute the optimal embedding K-->M_2(F)
    #we need two base change matrices
    w0,w1 = w.vector()
    BW=Matrix(base,[[1,w0],[0,w1]])
    WB=BW**-1
    #first, the image of w, the generator of O_K as O_F-module
    w2_base_W=WB*(w**2).vector()
    #Construct the matrix of multiplication by w in the basis [1,w]
    W=Matrix(base,[[0,w2_base_W[0]],[1,w2_base_W[1]]])

    if embedding == 1:
        S = Matrix(base,2,2,[-e,0,0,1])
        W=S**(-1)*W*S
        #W = other_embedding(base,ext,w)

    # Modify the matrix hoping for better imaginary parts
    A = Matrix(base,2,2,[e,0,0,1])**(ZZ(extra_parameters[0]))*Matrix(base,2,2,[0,-1,1,0])**ZZ(extra_parameters[1])*Matrix(base,2,2,[e,0,0,1])**(ZZ(extra_parameters[2]))
    # a,b,c,d=W.list()
    # # If A = [x,y,z,t] and W = [a,b,c,d]
    # # Want to minimize (b*z + d*t)*z - (a*z + c*t)*t
    # # We may choose t=2b and z=(a-d)
    # if a!=d:
    #     A = matrix_from_bottom_row(base,a-d,2*b)
    # else:
    #     A = matrix_from_bottom_row(base,1,2*b)
    W = A*W*A**(-1)

    assert (w.minpoly()).subs(W)==0
    #now, since beta=(w-w[0])/w[1] we can compute the image of genK
    Beta = (W-w0)/w1
    #now we compute the fixed point tau_0
    aa,bb,cc,dd = Beta.apply_map(v0).list()
    if cc > 0:
        tau0 = ((aa-dd)+sqrt((dd-aa)**2+4*bb*cc))/(2*cc)
    else:
        tau0 = ((aa-dd)-sqrt((dd-aa)**2+4*bb*cc))/(2*cc)
    assert tau0.imag() > 0
    # print('tau=',tau0)

    #now we can compute the matrix gamma_tau (notation as in Darmon-Logan) and its image under v0
    gtau = eps.vector()[0]+eps.vector()[1]*Beta
    #and gamma_tau\cdot\infty and the continued fraction
    print('tau0,gtau')
    print(tau0)
    print(gtau)

    return tau0,gtau

def act_single(M,Z,v):
    try:
        M0=M.apply_morphism(v) #Matrix(QQbar,2,2,[v(M[ii,jj]) for ii in range(2) for jj in range(2)])
        if v(M.determinant())>0:
            conj=False
        else:
            conj=True
        if Z is Infinity:
            if(M0[1,0]==0):
                return Infinity
            else:
                tmp=M0[0,0]/M0[1,0]
                if conj:
                    tmp=tmp.conjugate()
                return tmp
        try:
            tmp=(M0[0,0]*Z+M0[0,1])/(M0[1,0]*Z+M0[1,1])
            if conj:
                tmp=tmp.conjugate()
            return tmp
        except ValueError:
            print('!!')
            return Infinity
    except RuntimeError:
        print('!! Maybe too much depth in recursion?')
        raise RuntimeError

def get_img_part(F,tau_0,v0,v0approx,ci,cip1,return_lims=False):
    assert False
    GA=Matrix(ZZ,2,2,[1,-1,0,1])
    GB=Matrix(ZZ,2,2,[1,-1,1,0])
    # Once the chain is found, get the limits from it
    M=Matrix(F,[ci[0],cip1[0],ci[1],cip1[1]])
    if M.determinant()==-1:
        M.set_col_to_multiple_of_col(1,1,-1)
    assert M.determinant() == 1
    M = M**(-1)
    if return_lims:
        x,y =act_single(GA*M,tau_0,v0), act_single(GB*M,tau_0,v0)
        x0,y0=ComplexField()(x),ComplexField()(y)
        return x,y,RealField()(min([x0.imag(),y0.imag(),x0.imag()/x0.norm(),y0.imag()/y0.norm()])).abs().sqrt()
    else:
        x,y =act_single(GA*M,tau_0,v0approx), act_single(GB*M,tau_0,v0approx)
        x0,y0=ComplexField()(x),ComplexField()(y)
        return RealField()(min([x0.imag(),y0.imag(),x0.imag()/x0.norm(),y0.imag()/y0.norm()])).abs().sqrt()

def get_img_part_vec(F,tau_0,v0,v0approx,vec):
    assert False
    img_part = 10
    chain=[[1,0]]
    p = [0, 1]
    q = [1, 0]
    for i in range(len(vec)):
        p.append(p[i+1]*vec[i]+p[i])
        q.append(q[i+1]*vec[i]+q[i])
        chain.append([p[i+2],q[i+2]])
        new_img_part = get_img_part(F, tau_0, v0,v0approx, chain[i],chain[i+1],return_lims = False)
        img_part = min([img_part,new_img_part])
    return img_part


#atencio: aquesta rutina esta tunejada per a que retorni una unitat que ens va be per al cos amb alpha=r+1
#la normal és descomentant la penultima linia, on hi ha el que retornavem abans...
#tambe hi ha uns prints que sobren...
def find_the_unit_of(F,K):
    found=False
    for eps in K.units():
        is_square,root=eps.norm(F).is_square(root=True)
        deg=eps.minpoly().degree()
        if deg==2:
            unit_not_in_F=eps
        if is_square and deg==2:
            print('unit=',(ZZ(1)/root)*eps)
            return (ZZ(1)/root)*eps
    # Not found so far: the norm of unit_not_in_F is a unit of F which is not a
    # square, and unit_not_in_F**2/norm has norm 1.
    norm=unit_not_in_F.norm(F)
    return unit_not_in_F**2/norm


def power_in_Gamma_1(gamma, N):
    r"""
    Return ``(g, k)``, where ``k >= 1`` is the smallest exponent such that the
    top-left entry of ``gamma**k`` is congruent to `\pm 1` modulo ``N``, and
    ``g = \pm gamma**k`` is the corresponding element of `\Gamma_1(N)`.

    ``gamma`` must lie in `\Gamma_0(N)`. Then the top-left entry of ``gamma**k``
    is `a^k` modulo ``N``, where `a` (the top-left entry of ``gamma``) is a unit
    modulo ``N``; so ``k`` divides the order of `a` in `(O_F/N)^\times`.

    EXAMPLES::

        sage: from atrpoints.utility import power_in_Gamma_1
        sage: F.<r> = QuadraticField(5)
        sage: g = matrix(F, 2, 2, [3, 2, 7, 5])     # 3^3 = -1 mod 7
        sage: h, k = power_in_Gamma_1(g, F.ideal(7))
        sage: k, h == -g**3
        (3, True)
        sage: h
        [-181 -126]
        [-441 -307]
    """
    F = gamma.base_ring()
    idealN = F.ideal(N)
    assert gamma.det() == 1 and gamma[1,0] in idealN, 'gamma must lie in Gamma_0(N)'
    bound = idealN.euler_phi()
    g = gamma
    for k in range(1, bound + 1):
        if (g[0,0] - 1) in idealN:
            return g, k
        if (g[0,0] + 1) in idealN:
            return -g, k
        g = g * gamma
    raise RuntimeError('no power of gamma has top-left entry +-1 mod N')
    

def module_generators(K):
    r"""
    Return ``[1, w]`` with O_K = O_F + O_F w, for K/F a relative quadratic extension
    (O_K is a free O_F-module containing 1, e.g. when F has class number one).

    Write the elements of O_K as a + c*beta (beta = ``K.gen()``, a, c in F). The
    coefficients c form a fractional ideal J = (j) of F; take x in O_K with
    coefficient j, and reduce its F-part modulo O_F (as an HNF basis would).

    (This used to be computed with PARI's ``rnfhnfbasis`` on ``K.pari_rnf()``, but
    current Sage builds that rnf from a different relative polynomial, so the basis
    was expressed in the wrong generator.)
    """
    F = K.base_field()
    OF_basis = F.ring_of_integers().basis()
    to_OF = Matrix(QQ, [list(o.vector()) for o in OF_basis]).inverse()  # F-coords -> O_F-coords
    def coords(t):
        # coordinates of t in F with respect to the integral basis of F
        return vector(QQ, list(t.vector())) * to_OF
    # Absolute coordinates (c-part first) of a Z-basis of O_K
    rows = []
    for x in K.maximal_order().basis():
        a, c = K(x).vector()
        rows.append(list(coords(c)) + list(coords(a)))
    L = span([vector(QQ, r) for r in rows], ZZ)   # O_K as a lattice, echelon basis
    E = L.basis()
    assert len(E) == 4 and all(e[0] == 0 and e[1] == 0 for e in E[2:])
    J = F.ideal([sum(e[k] * OF_basis[k] for k in range(2)) for e in E[:2]])
    if not J.is_principal():
        raise ValueError('O_K is not of the form O_F + O_F w')
    jgen = J.gens_reduced()[0]
    s = Matrix(QQ, [list(e[:2]) for e in E[:2]]).solve_left(coords(jgen))
    assert all(t in ZZ for t in s)
    xvec = s[0] * E[0] + s[1] * E[1]
    a = sum(xvec[2 + k] * OF_basis[k] for k in range(2))
    ca = coords(a)
    a = sum((ca[k] - ca[k].floor()) * OF_basis[k] for k in range(2))  # reduce a modulo O_F
    w = a + jgen * K.gen()
    assert w.is_integral()
    assert F.ideal(w.minpoly().discriminant()) == K.relative_discriminant()
    return [K(1), w]

def find_embeddings(F,K):
    embF=F.real_embeddings()
    embK=K.real_embeddings()
    if abs(embK[0](F.gen())-embF[0](F.gen()))<ZZ(2)**(-30):
        v1=embF[0]
        v0=embF[1]
        swapped_embeddings=True
    else:
        v0=embF[0]
        v1=embF[1]
        swapped_embeddings=False
    v00=v0(F.gen())
    v10=v1(F.gen())
    v00exact=QQbar.polynomial_root(F.polynomial(),RIF(v00-ZZ(1)/10,v00+ZZ(1)/10))
    v10exact=QQbar.polynomial_root(F.polynomial(),RIF(v10-ZZ(1)/10,v10+ZZ(1)/10))
    v0=F.hom([v00exact],QQbar)
    v1=F.hom([v10exact],QQbar)
    v0approx=F.hom([RealField(prec_bits)(v00exact)],RealField(prec_bits),check=False)
    v1approx=F.hom([RealField(prec_bits)(v10exact)],RealField(prec_bits),check=False)
    return v0,v1,v0approx,v1approx,swapped_embeddings

def find_fundamental_unit(F,v0,v0approx):
    e=F.units()[0]
    if v0(e)<0:
        e=-e

    if v0approx(e**2)<1:
        e = ZZ(1)/e
    return e


#rutina auxiliar
#fa la divisio entera entre a i b; retorna el quocient i el residu
#com més gran sigui radi més probabilitat hi ha que trobi el quocient
#si no el troba retorna l'string 'not found'
def division(a,b,radi=100):
    quo=a/b
    w = quo.parent().ring_of_integers().ring_generators()[0]
    # a = quo.numerator()
    # b = quo.denominator()
    compab=quo.vector()
    comp=[compab[0].n().floor(), compab[1].n().floor()]
    #la primera aproximacio es la obvia, fer part entera de les components de a/b
    first_aprox=comp[0]+comp[1]*w
    if ((quo-first_aprox).norm())**2<1: 
        return [first_aprox,a-b*first_aprox]
    #si no funciona la first_aprox, mirem d'aproximar a/b per elements del lattice <1,w>
    #anem buscant punts com a molt a distancia r (en norma L1) i fem creixer r fins al parametre radi 
    for r in range(0,radi,5):
        for i in range(-r,r):
            for j in range(-r,r):
                if i**2>=(r-5)**2 or j**2>=(r-5)**2:
                    aprox=first_aprox+i+j*w
                    if ((quo-aprox).norm())**2<1:
                        return [aprox,a-b*aprox]
    raise RuntimeError('Not found...')

#rutina auxiliar que, donat un vector v, calcula el quocient de la fraccio continua associat a v[0],v[1],...,v[n].
#la faig servir per a comprovar que la fraccio continua que calculem esta ben feta
def check_cf(v):
    if len(v)==2:
            return v[0]+ZZ(1)/v[1]
    else:
            w=[]
            for i in range(1,len(v)):
                    w.append(v[i])
            return v[0]+ZZ(1)/check_cf(w)

#aquesta es la rutina que calcula l'expressio en fraccions continues de a/b fent servir nomes divisions de longitud 1
def cf(quo):
    a = quo.numerator()
    b = quo.denominator()
    w = quo.parent().ring_of_integers().ring_generators()[0]
    v=[]
    knm1=a
    kn=b
    while True:
        an,knp1=division(knm1,kn)
        v.append(an)
        if knp1==0:
            break
        knm1=kn
        kn=knp1
    # # aqui comprovem que a/b coincideix amb la fraccio continua associada a v
    # if check_cf(v)!=a/b:
    #   raise RuntimeError, 'Error: no ha calculat be la fraccio continua!'
    return v

#this is the main function: given a matrix M in Gamma_1(N), it returns its decomposition into 6 elementary matrices
#(actually, it does not use that M belongs to Gamma_1(N)...does this mean that we can decompose any matrix in SL_2 as a product of six elementray matrices???)
def el_mat_dec(M,level):
    F=M.base_ring()
    u=F.units()[0]
    E_lambda,uu=find_lambda(M,u)
    #we know that E_lambda*M is a matrix [a,b,c,d]
    # such that c=uu+ta for some unit uu; now we find uu and t
    MM=E_lambda*M
    a,b,c,d = MM.list()
    e = 1
    # e = is_represented_by_unit(F,level,a,F.units()[0])
    t=(c-uu)/a
    assert t.is_integral()
    E1i=Matrix(F,2,2,[1,0,uu*(e-a)/e,1])
    E2i=Matrix(F,2,2,[1,-ZZ(1)/uu,0,1])
    E3i=Matrix(F,2,2,[1,0,t*(e-a),1])
    E33i = Matrix(F,2,2,[1,0,c/e,ZZ(1)/e])
    T=(E33i*E3i*E2i*E1i)**(-1)
    E_x=T*MM
    # assert E_x[1][0] == 0 and E_x[0][0] == 1
    dec = [E_lambda**(-1), E33i*E3i, E2i, E1i, E_x]
    assert prod(dec) == M
    return dec

#given a matrix M=[a,b,c,d] we run over all lambda in O until we find some satisfying that the class of c mod a+lambda*c is represented by a unit; observe that we are not requiring that a+lambda*c is a prime, nor that the units generate all the quotient ring: if the class of c (mod a+lambda*c) can be represented by a unit, that's fine
#Rmk: maybe here's room for improvement: the loop is pretty pedestrian, and moreover it's not clear at all that we're running over the best lambda'a possible...but apparently it works so far...
def find_lambda(M,u):
    F=M.base_ring()
    n=0
    a,b,c,d=M.list()
    uu = dict()
    for i in range(-2,2):
        uu[i]=u**i
    while True:
        n += 1
        for l,i,sign in product(F.elements_of_norm(n),list(range(-2,2)),[-1,1]):
            lmd=sign*l*uu[i]
            rep_unit = is_represented_by_unit(F,a+lmd*c,c,u)
            if rep_unit != False:
                return Matrix(F,2,2,[1,lmd,0,1]),rep_unit


#we want to know if the class of c mod a is represented by a unit in F
#returns the unit representing c mod a it it is, and returns False otherwise
def is_represented_by_unit(F,a,c,u):
    A = F.ideal(a)
    c_coords = A.coordinates(c)
    one_coords = A.coordinates(1)
    if (c_coords-one_coords).denominator() == 1:
        return 1
    uinv = ZZ(1)/u
    n=1
    unit1 = 1
    unit2 = 1
    while True:
        unit1 *= u
        unit2 *= uinv
        u1_coords = A.coordinates(unit1)
        u2_coords = A.coordinates(unit2)
        if (c_coords - u1_coords).denominator() == 1:
            return unit1
        if (c_coords - u2_coords).denominator() == 1:
            return unit2
        if (c_coords + u1_coords).denominator() == 1:
            return -unit1
        if (c_coords + u2_coords).denominator() == 1:
            return -unit2
        if (u1_coords - one_coords).denominator() == 1:
            return False
        n+=1

#given a 2x2 matrix M, it checks whether it belongs to Gamma_1(N)
def is_in_Gamma_1(elt,N):
    F=elt.base_ring()
    II = F.ideal(N)
    return II.coordinates(F(elt[0,0])-1).denominator() == 1 and II.coordinates(F(elt[1,0])).denominator() == 1 and elt.det() == 1

def is_in_Gamma_0(elt,N):
    F=elt.base_ring()
    return F.ideal(N).coordinates(F(elt[1,0])).denominator() == 1 and elt.det() == 1

def quotients_and_limits_withlevel(F,tau0,v0,v0approx,v1,gtau,level):
    r''' We find the limits of the integrals to be computed. These are 3-component limits.
    That is, it converts the integral \int^tau_0\int_\infty^cusp into a sum of integrals
    \int^tau_1\int_infty^0+int^\tau_2\int_\infty^0+int^\tau_3\int_\infty^0...then we use
    (a variation of) the trick of page 21 in Darmon-Logan to convert this into a sum of
    integrals of the form \int_x1^y1\int_0^\infty+\int_x2^y2\int_0^\infty+\int_x3^y3\int_0^\infty+...
    and we obtain [[x1,y1],[x2,y2],[x3,y3],...]'''
    I = QQbar(-1).sqrt()
    N0 = v0(level)
    N1 = v1(level)
    sqN1 = N1.sqrt()

    # Note: we assume that gtau is in Gamma1 !
    V = []
    tau = tau0
    # assert is_in_Gamma_1(gtau,level)
    for E in el_mat_dec(gtau,level):
        a,b,c,d = E.list()
        # assert a == 1 and d == 1
        newtau = act_single(E.inverse(),tau,v0)
        if b != 0: # upper-triangular
            assert c == 0
        else: # lower-triangular
            assert b == 0
            if v0approx(E.determinant()) > 0:
                conj = 1
            else:
                assert 0
                conj = -1
            V.extend([Limit(tau,I/sqN1,newtau,Infinity, conj = conj),Limit(-ZZ(1)/(N0*newtau),I/sqN1,-ZZ(1)/(N0*tau),Infinity, conj = conj)])
        tau = newtau
    return V,min([l.imaginary_part_approx() for l in V])


# this function gives a sequence of cusps connecting \infty to the cusp x; the output is a vector in which each component is of the form [a,b], being a/b the corresponding cusp. The first component is always [1,0], the cusp \infty
def quotients_and_limits(F,tau0,v0,v0approx,v1,cusp_cf = None,max_length = -1,threshold = None, assume_E1 = False, level = 1):
    r''' We find the limits of the integrals to be computed. These are 3-component limits.
    That is, it converts the integral \int^tau_0\int_\infty^cusp into a sum of integrals
    \int^tau_1\int_infty^0+int^\tau_2\int_\infty^0+int^\tau_3\int_\infty^0...then we use
    (a variation of) the trick of page 21 in Darmon-Logan to convert this into a sum of
    integrals of the form \int_x1^y1\int_0^\infty+\int_x2^y2\int_0^\infty+\int_x3^y3\int_0^\infty+...
    and we obtain [[x1,y1],[x2,y2],[x3,y3],...]'''
    print('level =',level)
    # level = -level
    I = QQbar(-1).sqrt()
    opt_V = []
    opt_imag_part = -1
    N0 = v0(level)
    N1 = v1(level)
    sqN1 = N1.sqrt()
    idealN = F.ideal(N)
    cusp = gtau[0,0]/gtau[1,0]

    # First, we find the chains
    if cusp_cf is None:
        if assume_E1 == True:
            print('Assuming E1!')
            all_vecs=[cf(cusp)]
        else:
            if max_length==-1:
                all_vecs=[quadratic_continued_fraction(cusp,optimize_for_length=False)]
            elif  max_length==0:
                all_vecs=quadratic_continued_fraction(cusp,optimize_for_length=True)
            else:
                all_vecs=all_qcf_up_to_length(cusp,max_length=max_length, get_img_part = lambda vec: get_img_part_vec(F, tau0, v0, v0approx, vec),threshold = threshold)
    else:
        all_vecs=[cusp_cf]

    # print('Found %s candidates...'%len(all_vecs))
    all_chains=[]
    for vec in all_vecs:
        V = []
        img_part = 100
        is_good = True
        chain=[[1,0]]
        p = [0, 1]
        q = [1, 0]
        for i in range(len(vec)):
            p.append(p[i+1]*vec[i]+p[i])
            q.append(q[i+1]*vec[i]+q[i])
            chain.append([p[i+2],q[i+2]])
            x, y, new_img_part = get_img_part(F, tau0,v0,v0approx, chain[i],chain[i+1],return_lims = True)
            img_part = min([img_part,new_img_part])
            if img_part < opt_imag_part:
                is_good = False
                break
            if idealN.coordinates(chain[i][1]).denominator().abs() == 1:
                V.extend([Limit(x,I/sqN1,y,Infinity),Limit(-ZZ(1)/(N0*y),I/sqN1,-ZZ(1)/(N0*x),Infinity)])
            else:
                assert 0
                # V.extend([Limit(y,I/sqN1,x,Infinity),Limit(-1/(N0*x),I/sqN1,-1/(N0*y),Infinity)])

            # V.extend([Limit(-1/y,I,-1/x,Infinity),Limit(x,I,y,Infinity)])
        if is_good:
            opt_imag_part = img_part
            opt_V = V
    if opt_imag_part == -1:
        raise RuntimeError
    return opt_V, opt_imag_part


def bezout(F,a,b):
    if b == 0:
        return 1,0
    q=quadratic_continued_fraction(F(a/b))
    u=[]
    v=[]
    u.extend([1,0])
    v.extend([0,1])
    #print(len(q))
    for i in range(0,len(q)):
        u.append(u[i]-q[i]*u[i+1])
        v.append(v[i]-q[i]*v[i+1])
    return u[len(u)-2],v[len(v)-2]


def matrix_from_left_column(F,a,c):
    d,b=bezout(F,a,c)
    g=b*c+a*d
    if g.norm().abs() != 1:
        print('Not coprime inputs?')
        print('g = ',g)
        print('norm g =',g.norm().abs())
        raise RuntimeError

    Mres=Matrix(F.ring_of_integers(),2,2,[a,-b/g,c,d/g])
    assert Mres.det()==1
    return Mres

def matrix_from_bottom_row(F,c,d):
    a,mb=bezout(F,d,c)
    g=mb*c+a*d
    assert g.norm().abs()==1
    return Matrix(F.ring_of_integers(),2,2,[a/g,-mb/g,c,d])

def other_embedding(F,K, w):
    a0,b0=K.class_group().gens()[0].ideal().gens_reduced()
    vecw=w.vector()
    T=Matrix(F,2,2,[1,vecw[0],0,vecw[1]]).inverse()
    veca0=a0.vector()
    vecb0=b0.vector()
    veca=T*Matrix(F,2,1,veca0)
    vecb=T*Matrix(F,2,1,vecb0)
    S=Matrix(F,2,2,[veca[0,0],vecb[0,0],veca[1,0],vecb[1,0]])
    vec=T*Matrix(F,2,1,(w**2).vector())
    W=S**(-1)*Matrix(F,2,2,[0,vec[0,0],1,vec[1,0]])*S
    if not all([W[ii,jj].is_integral() for ii in range(2) for jj in range(2)]):
        a=F(a0)
        Ma=Matrix(F,2,2,[a,0,0,1])
        W=Ma*W*Ma**-1
    assert all([W[ii,jj].is_integral() for ii in range(2) for jj in range(2)])
    return W


def find_alpha(E, x_coord):
    F = E.base_ring()
    ff=E.defining_polynomial()
    Fx=PolynomialRing(F,names='Y')
    Y=Fx.gen()
    R=ff.parent()
    phi=R.hom([x_coord,Y,1])
    alpha=phi(ff).discriminant()
    alpha*=ZZ(alpha.denominator())**2
    vv=alpha.factor()
    if vv.unit().is_square():
        alpha=F(1)
    elif (-vv.unit()).is_square():
        alpha=F(-1)
    elif (vv.unit()/F.units()[0]).is_square():
        alpha=F.units()[0]
    elif (-vv.unit()/F.units()[0]).is_square():
        alpha=-F.units()[0]
    else:
        raise RuntimeError("Don't know what to do here...")
    for xx in vv:
        if xx[1]%2==1:
            alpha*=xx[0]
    return alpha
