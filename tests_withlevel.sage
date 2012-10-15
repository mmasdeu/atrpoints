'''
The file that runs the tests with LEVEL
'''

from sage.parallel.all import fork
from itertools import product
from operator import itemgetter, attrgetter
load 'atrpoints_withlevel.sage'
sys.setrecursionlimit(10**7)


F.<r>=QuadraticField(5)
w=(1+r)/2
# Curve used in J.Gartner's thesis (note typo in J.G. thesis!)
E=EllipticCurve([1,-(w+1),w,-(30*w+45),-(111*w+117)])

N=E.conductor().gen(1)


''' How to run a test
=========================================================

 1) Initialize TestDL(E,prec,working_prec), where:

 - E is the elliptic curve to work with.
 - prec is the desired output precision for the ATR point.
 - working_prec is the precision used to make the calculations
   (which should be larger than prec). If one wants to work
   with machine doubles (much faster) the working_prec should be 53.

 2) Run the method check_conjecture(x_coord, max_length), where:

 - x_coord is the x-coordinate of a known point of infinite order defined
   over an atr field.
 - max_length is the maximum length allowed for the continued fraction expansion.
   One should experiment with this value: too large takes too long to find, and too
   small will either not find any expansion, or will not find the optimal one.
'''


############################################
# A working example
############################################
alpha=-r + 1
R.<z>=PolynomialRing(F)
K.<gK>=F.extension(z^2-alpha)
ww=module_generators(K)[1]
K.<w>=NumberField(ww.minpoly())
KK.<gKK>=K.absolute_field()
x_coord=1/1210*(-1090*gKK^3 + 1943*gKK^2 + 316*gKK - 3422) 
EK=E.base_extend(F.embeddings(KK)[1])
P=EK.lift_x(x_coord)
# Note that x_coord does not belong to F, but K.
# However, 2*P has the x-coordinate in F.
alpha=-r + 1
R.<y>=PolynomialRing(F)
K.<beta>=F.extension(y^2-alpha)
x_coord_2=18883/2420*alpha - 16127/2420

# Note that we use the x-coordinate of 2*P
T = TestDL(E, prec = 40, working_prec=500,level = N)
T.check_conjecture(x_coord = x_coord_2, alpha = alpha, max_length = 8)
T.check_lindeps(9)

# RESULTS:

#sage: T.check_lindeps(9)
#Verification of the imaginary part of the point
#[-1, -2, -2]~
#Verification of the real part of the point (should be torsion)
#[1, -12]~

#sage: T.Jtau
#-4.8289548170777950654568971974198568649337686191168718743971233209322236 + 4.5346965323332380780806864274681858939853704738334706173669141904269964*I

#sage: T.max_norm
#[179710, 0]

#sage: T.lambda_0_minus_K
#-1.199422384661861113660086*I

#sage: T.lambda_1_minus_K
#1.5721667861326515211451*I

#sage: T.lambda_0_plus_K 
#-0.381521644169657283699700042

#sage: T.lambda_1_plus_K 
#1.054757485987466300211

#unit= (r - 1)*w - 3/2*r + 5/2

#gtau=Matrix(F,2,2,[3/2*r - 5/2,-r + 1,8*r - 14, -9/2*r + 15/2])

#Min Imag part after prepare_limits:
#0.0115563456513215
#Updated max_norm to 179710.000000000

#tau0= 0.4392914189919321? + 0.3534081297534779?*I


'''
sage: el_mat_dec(gtau,N)   
[
[           1 -1/2*r + 1/2]  [               1                0]
[           0            1], [-36692*r + 82047                1],

[                   1 75025/2*r + 167761/2]
[                   0                    1],

[                      1                       0]
[3356817/2*r - 7506071/2                       1],

[               1 -30150*r - 67418]
[               0      1/2*r + 3/2]
]
'''

#################################
