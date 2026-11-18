# The file that runs the tests with LEVEL

from sage.parallel.all import fork
from itertools import product
from operator import itemgetter, attrgetter
from darmonpoints.util import *
from darmonpoints.arithgroup_generic import *
from tqdm import tqdm
sys.setrecursionlimit(10**7)

prec = 1000
R = 20
F.<r>=QuadraticField(5)
w = (1+r)/2
eps = (1+r)/2
gamma1 = matrix(F,2,2,[4*eps+2, -2*eps-5, -2*eps-1, 2*eps+1])
gamma2 = matrix(F, 2, 2, [13*eps+9, 4*eps+1, -32*eps-18, -7*eps-6])
gamma3 = matrix(F,2,2,[-47*eps-17, 9*eps-6, -520*eps-308, 53*eps+20])
gamma4 = matrix(F,2,2,[165*eps+79, 5*eps-2, -8512*eps-5160, -159*eps - 76])
assert all(gamma.determinant() == 1 for gamma in [gamma1, gamma2, gamma3, gamma4])

gamma = gamma1


from atrpoints.arithgrouprq import *
G = ArithGroupRQ(F, use_shapiro=False)

emb = lambda x : x
H = OneChains(G, ZZ**1, emb=emb)

Gab = G.abelianization()
g2wordrep = G(gamma**2).word_rep
gab = Gab(G(gamma**2))

# comm_list = G.decompose_into_commutators(G(gamma),2)


id=matrix(F,2,2,[1,0,0,1])


divisors(4*norm(eps^2-1))

def min_exp(gamma):
    for d in divisors(4*norm(eps^2-1)):
        if vector(Gab(G(gamma**d)))==0:
            return d

a, b, c, d = gamma.list()
deltaaa=gamma.trace()^2-4
t=QQ['t'].gen()
K.<h>=F.extension(t^2-(21*eps-11))
tau_CD= (a-d+K(deltaaa).sqrt())/(2*c)


OF = F.ring_of_integers()
disc = F.absolute_discriminant()

v0, v1 = F.real_embeddings(prec)
found = False
for w0 in K.complex_embeddings(prec):
    if w0(h).imag().abs() > 10**-10:
        found = True
        break
assert found, "Couldn't find a complex embedding of K"
v0, v1  = F.real_embeddings(prec)
if (w0(K(F.gen())) - v0(F.gen())).abs() > 10**-10:
    v0, v1 = v1, v0
assert (w0(K(F.gen())) - v0(F.gen())).abs() < 10**-10, "Couldn't find a complex embedding of K that extends the chosen real embedding of F"
real = [v0, v1]
CC = ComplexField(prec)
pi = CC(pi)
zfm1 = K.zeta_function(prec)(-1)
rtdisc = CC(disc).sqrt()
twoIpi = CC.gen() * 2 * pi
delta=r
RF= K.regulator()

etha=eps^2
f=1

def reduction(w1,w2):
    T, w = abs(w2/w1), w1
    if(real[0](w2)*real[1](w2)>0):
        T, w = abs(w1/w2), w2
    if T < 0:
        T = -T
    if T[1]<0 or (floor(T)<T[0]<ceil(T) and T[1] < min(a-floor(a),ceil(a)-a)):
        return 0
    return [w,-floor(T-2*T[0]*r)+T]

def reduction_integers(g):
    if disc%4 == 1:
        return reduction(1,(1+g)/2)
    if (disc//4)%4 == 2 or (disc//4)%4 == 3:
        return reduction(1,g)
    
tau=reduction_integers(r)[1]
bk_liste = continued_fraction_list(real[1](tau), type = 'hj', nterms=20)
                                   
A_LAMBDA_tau = []
A0 = F(1)
Ak = A0
Ak_moins_1 = tau
k = 0

while Ak_moins_1 != etha^(-1):
    A_LAMBDA_tau.append(Ak) 
    Ak_plus_1 = bk_liste[k]*Ak - Ak_moins_1
    Ak_moins_1 = Ak; Ak = Ak_plus_1; k += 1

S2=[]; n=0; 




while abs(real[0](etha^n))+abs(real[1](etha^n))<R or abs(real[0](etha^n*A_LAMBDA_tau[1]))+abs(real[1](etha^n*A_LAMBDA_tau[1]))<R:
    for j in range(0,len(A_LAMBDA_tau)-1):
        a=1; b=0
        while abs(real[0](etha^n*a*(A_LAMBDA_tau[j])))+abs(real[1](etha^n*a*(A_LAMBDA_tau[j])))<R:
            T=a*A_LAMBDA_tau[j]*etha^n+b*A_LAMBDA_tau[j+1]*etha^n
            while abs(real[0](T))+abs(real[1](T))<R:
                S2.append(T)
                b=b+1
                T=a*A_LAMBDA_tau[j]*etha^n+b*A_LAMBDA_tau[j+1]*etha^n
            a=a+1; b=0
    n=n+1
n=-1
while abs(real[0](etha^n))+abs(real[1](etha^n))<R or abs(real[0](etha^n*A_LAMBDA_tau[1]))+abs(real[1](etha^n*A_LAMBDA_tau[1]))<R:
    for j in range(0,len(A_LAMBDA_tau)-1):
        a=1; b=0
        while abs(real[0](etha^n*a*(A_LAMBDA_tau[j])))+abs(real[1](etha^n*a*(A_LAMBDA_tau[j])))<R:
            T=a*A_LAMBDA_tau[j]*etha^n+b*A_LAMBDA_tau[j+1]*etha^n
            while abs(real[0](T))+abs(real[1](T))<R:
                S2.append(T)
                b=b+1
                T=a*A_LAMBDA_tau[j]*etha^n+b*A_LAMBDA_tau[j+1]*etha^n
            a=a+1; b=0
    n=n-1
n=0

def aire(A,B):
    return(-sign(A[1][0]*B[1][0]*((A*B)[1][0])))

def matrix_real(A,i):
    return(matrix([[real[i](A[0][0]),real[i](A[0][1])],[real[i](A[1][0]),real[i](A[1][1])]]))

def moeb_real(tau,A,i):
    a, b, c, d = A.apply_morphism(real[i]).list()
    ans = (a * tau + b ) / (c * tau + d)
    return ans

def int_exp(P,A,B,mu):
    mu_delta0 = real[0](mu/delta)
    mu_delta1 = real[1](mu/delta)
    z0 , z1 = P
    exp1 = (twoIpi*mu_delta0*moeb_real(z0,A,0)).exp()-(twoIpi*mu_delta0*z0).exp()
    exp2 = (twoIpi*mu_delta1*moeb_real(z1,A*B,1)).exp()-(twoIpi*mu_delta1*moeb_real(z1,A,1)).exp()

    return exp1*exp2/ CC(norm(mu))

def non_eis_int(z1,z2):
    if (z1*z2).imag()==0:
        if z1.imag()!=0:
            return I/2 * (z2.real()-z1.real())*(1/z1.imag())
        if z2.imag()!=0:
            return I/2 * (z2.real()-z1.real())*(-1/z2.imag())
        return 0
    return I/2 * (z2.real()-z1.real())*(1/z1.imag()-1/z2.imag())

@cached_function
def sigma_1(mu):
    return sum([ZZ(abs(norm(A))) for A in divisors(F.ideal(mu))])

def kappa_square(P,A,B):
    suma = sum(sigma_1(mu)* int_exp(P,A,B,mu) for mu in S2)
    z0, z1 = P
    const = (twoIpi)^2*zfm1*(moeb_real(z0,A,0)-z0)*(moeb_real(z1,A * B,1)-moeb_real(z1,A,1))/ rtdisc
    extr_term = 2*RF * (non_eis_int(z0,moeb_real(z0,A,0))-non_eis_int(moeb_real(z1,A,1),moeb_real(z1,A*B,1)))
    return const + 4* rtdisc * suma + extr_term

def kappa(P, A, B):
    k2 = kappa_square(P, A, B)
    C = RF * pi * CC.gen()
    A0 = A.apply_morphism(real[0])
    B0 = B.apply_morphism(real[0])
    A1 = A.apply_morphism(real[1])
    B1 = B.apply_morphism(real[1])
    return k2 + C * (aire(A1, B1) - aire(A0, B0))

def rho_abelian(G, gamma,P):
    bc = G.find_bounding_cycle(G(gamma), npow=2)
    bc = sum(([(n, a, b) for a in alist] for n, alist, b in bc), [])
    ans = 0
    for n, A, B in tqdm(bc):
        ans += -n * kappa(P, A, B)
    return ans

def rho(G, gamma,tau):
    n = min_exp(gamma)
    P = [w0(tau), CC.gen()]
    return rho_abelian(G, gamma**n, P) / (n * twoIpi**2)

rho_gamma_tau = rho(G, gamma, tau_CD)
print(rho_gamma_tau)

rho_good = -.36666666666666666 - .27784944302 * CC.gen()
Jgood = CC(-5.7303667736837191510384895962270519334997028174)
J = (rho_gamma_tau * twoIpi).exp()
for z, _ in (x^15-1).roots(CC):
    print(algdep(J*z,16).degree())