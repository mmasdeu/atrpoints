'''
The file that runs the tests
'''

from sage.parallel.all import fork
from itertools import product
from operator import itemgetter, attrgetter
from atrpoints import *
sys.setrecursionlimit(10**7)

F.<r>=QuadraticField(29)
w=(1+r)/2
eps=2+w
E29=EllipticCurve(F,[1,0,eps^2,0,0])
x_coord_29=[(-1+w)+3,(2+w)/2,(11*(17+8*w)+5)/8,(2*(19+9*w)+1)/5,(-4*(4+3*w)-11)/15,-1/9,7*(3*w)/9+5,-1/4,(43*(1+w)+51)/10,(98*(7+5*w)+387)/13,(-3*(-5+5*w)-13)/5]

F.<r>=QuadraticField(37)
w=(1+r)/2
E37=EllipticCurve(F,[0,2,1,-(19+8*w),28+11*w])
x_coord_37=[-2*(w-3)/3-13/3,(w+1)/7-3/7,-2*(15*w+38)/165-104/165,(2*w+5)/8-5/8,115*(w+2)/588-80/147,-(4*w+10)/8-3/4,-196*(5*w-15)/675-20/9]

F.<r>=QuadraticField(41)
w=(1+r)/2
E41=EllipticCurve(F,[1,0,0,-27-10*w,0])
x_coord_41=[-1/4,(-3*(-181-67*w)-1481)/268,((697+258*w)-9)/43,(-(-697-258*w)-1729)/258,(-7102*(389+144*w)-1271153)/9884736,(29*(1+w)+49)/4]

F.<r> = QuadraticField(109)
w=(1+r)/2
E109 = EllipticCurve(F,[w,-(1+w),0,-(58*w+245),-(630*w+2944)])

F.<r> = QuadraticField(509)
w=(1+r)/2
E509 = EllipticCurve(F,[-1, 2+2*w, -w, 162+3*w, 71+34*w])


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

T = TestDL(E29, prec = 20, working_prec=53)
T.check_conjecture(x_coord_29[2], max_length = 8)

a = E109.base_ring().gen()
w = (1+a)/2
T = TestDL(E109, prec = 14, working_prec=53)
T.check_conjecture(3*w + 11, max_length = 8)
