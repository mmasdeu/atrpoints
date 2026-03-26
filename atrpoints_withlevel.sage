prec_bits=2000
TIME_LIMIT=10
LOW_PREC_BITS=300
LOWER_PREC_BITS=100
from itertools import islice,product
from sage.misc.cachefunc import disk_cached_function
import pickle
from sage.parallel.all import fork
import gc
from array import array
from bisect import bisect_left
NCPUS=sage.parallel.ncpus.ncpus()
from sage.parallel.all import fork
sys.setrecursionlimit(10**7)
from quadcontfrac import *
load('atr_cython_withlevel.spyx')
load('utility.sage')

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

#this is a Hilbert Modular Form
class HilbertModularForm(SageObject):
#we use an elliptic curve to initialize a HMF f; we have to pass the level of f, and we compute the Fourier coefficients of ideals of norm up to max_norm
    def __init__(self,E,max_norm_aprs=1,level=1,n_cpus=None,save_to_disk='',expand_factor=1,prec=53,large_prime_data = None):
        if type(E) is str:
            self.E,self.max_norm_aprs,self.N,self.F,self.D,self.aprs = load(E)
            self.init_embeddings()
        else:
            # assert(level==1)
            self.N=level
            self.E=E
            self.F=self.E.base_field()
            self.D=self.F.disc()
            assert(2**valuation(expand_factor,2)==expand_factor)
            self.aprs=dict({ZZ(1):[[(1,1,1)]]})
            self.init_embeddings()
            if n_cpus is None:
                n_cpus=sage.parallel.ncpus.ncpus()
            print('n_cpus=%s'%n_cpus)
            assert(2**valuation(n_cpus,2)==n_cpus)
            M=n_cpus*expand_factor
            congN=2*M
            pari.init_primes(2*max_norm_aprs)
            self.small_norm,rem=ZZ(max_norm_aprs).sqrtrem()
            if rem != 0:
                self.small_norm+=1

            width=ceil((max_norm_aprs-self.small_norm)/QQ(M))
            Ainv=self.E.a_invariants()
            e2=self.e2
            v00,v10=self.v00,self.v10 #F.real_embeddings(prec)
            if prec==53:
                v0=lambda x:RDF(v00(x))
                v1=lambda x:RDF(v10(x))
            else:
                v0=v00
                v1=v10

            self.v0=v0
       
            self.v1=v1
            if v0(e2).abs()<1:
                e2=1/e2
            print('M=',M)
            # inputs=[(self,max_norm_aprs,congN,1,[2])]+[(self,max_norm_aprs,congN,a%congN,[]) for a in range(3,congN,2)]
            print('Computing data for small primes... (up to %s)'%self.small_norm)
            level = self.N
            newdict=dict(compute_coefficients(self.F,2,self.small_norm,max_norm_aprs,0,Ainv,level,v0,v1,e2,prec,single = False))
            newdict[1]=[[(1,1,1,1)]]
            newvec=[None,[(1,1,1,1)]]
            for nn in range(2,self.small_norm):
                ff=factor(nn)
                if len(ff)==1:
                    try:
                        newvec.append(newdict[ff[0][0]][ff[0][1]-1])
                    except KeyError:
                        newvec.append([])
                else:
                    if any([r%2==1 and kronecker_symbol(self.D,p)==-1 for p,r in ff]):
                        newvec.append([])
                    else:
                        newvec.append(get_all_combinations_small(self,ff,newdict))

            self.aprs_small=newdict
            large_primes=filter(lambda xx:kronecker_symbol(self.F.discriminant(),xx)!=-1,prime_range(self.small_norm,max_norm_aprs))
            if large_prime_data is None:
                print('Starting the massive computation now...')
                #large_primes=[xx for xx in prime_range(self.small_norm,max_norm_aprs) if kronecker_symbol(self.F.discriminant(),xx)!=-1]
                large_prime_data=[None]*len(large_primes)
                inputs=[]
                for ii in range(self.small_norm,max_norm_aprs,width):
                    pp=next_prime(ii-1)
                    while kronecker_symbol(self.F.discriminant(),pp)==-1:
                        pp=next_prime(pp)
                    inputs.append((self.F,ii,min([ii+width,max_norm_aprs]),max_norm_aprs,bisect_left(large_primes,pp), Ainv,level,v0,v1,e2,prec,True))
                    # print('[%s,%s)'%(ii,min([ii+width,max_norm_aprs])))
                print( 'len(inputs)=',len(inputs))
                I=compute_coefficients(inputs)
                gc.disable()
                kk=0
                t=walltime()
                for res in I:
                    for data in res[1]:
                        large_prime_data[data[0]]=data[1]
                        # This line should be removed for efficiency, after debugging
                        # assert data[0]==bisect_left(large_primes,data[2])
                    print('Finished %s-th input out of %s. (%s) -- mem = %s'%(kk,len(inputs),walltime(t),get_memory_usage()))
                    t=walltime()
                    kk+=1
                    if kk%10==0:
                        gc.collect()
                gc.enable()
            self.aprs_large=large_prime_data
            self.ans_small=newvec
            self.large_primes=large_primes
            self.max_norm_aprs=max_norm_aprs

    def init_embeddings(self):
        v0,v1=self.F.real_embeddings()
        self.e=self.F.units()[0]
        if not v0(self.e).abs()>1:
            v0,v1=v1,v0
        assert(v0(self.e).abs()>1 and v1(self.e).abs()<1)
        self.e2=self.e**2

        v00=v0(self.F.gen())
        v10=v1(self.F.gen())
        self.v00=v0
        self.v10=v1
        v00exact=QQbar.polynomial_root(self.F.polynomial(),RIF(v00-1/10,v00+1/10))
        v10exact=QQbar.polynomial_root(self.F.polynomial(),RIF(v10-1/10,v10+1/10))
        self.v0approx=self.F.hom([RealField(prec_bits)(v00exact)],RealField(prec_bits),check=False)
        self.v1approx=self.F.hom([RealField(prec_bits)(v10exact)],RealField(prec_bits),check=False)

    def save_form(self,filename):
        save([self.E,self.max_norm_aprs,self.N,self.F,self.D,self.aprs],filename)

    def load_coefficients(self,filename):
        self.aprs=load(filename)

    def equalize(self,x):
        x0,x1=self.v0approx(x),self.v1approx(x)
        if not (x0>0 and x1>0):
            x=-x
            x0,x1=self.v0approx(x),self.v1approx(x)
        e0=self.v0approx(self.e2).abs()
        while x0*e0<x1:
            x*=self.e2
            x0,x1=self.v0approx(x).abs(),self.v1approx(x).abs()
        while x0>x1*e0:
            x/=self.e2
            x0,x1=self.v0approx(x).abs(),self.v1approx(x).abs()
        return x



#this class represents an ATR extension of F. Thus, we pass the field F and an element alpha belonging to F such that the ATR extension is K=F(sqrt(alpha))
class ATRField(SageObject):
    def __init__(self,F,alpha,C = None, embedding=0, extra_parameters = None,level = 1):
        self.minW=1
        self.cached_ints=dict()
        self.cached_ints_N=0
        #base is the base field F
        self.base=F
        self.I=QQbar(-1).sqrt()
        self.genF=self.base.gen()
        self.alpha=alpha
        self.level = level
        if not self.is_ATR():
            raise ValueError,'Error: the extension is not ATR'
        Q.<y>=PolynomialRing(F)
        #ext is the field extension K, with generator beta
        self.ext=F.extension(y^2-alpha,'beta')

       	#we create the embeddings v0,v1 of F, so that v0 extends to a complex embedding of K and v1 extends to a real embedding of K
        self.v0,self.v1,self.v0approx,self.v1approx,self.swapped_embeddings = find_embeddings(F,self.ext)
        #e is the fundamental unit, chosen such that v0(e)>0 and v1(e)<0 and v0(e)>1.
        self.e = find_fundamental_unit(F,self.v0,self.v0approx)
        self.e2=self.e**2

        self.e0=self.v0(self.e)
        self.e0approx=self.v0approx(self.e)
        self.epssqrt=((RealField(prec_bits)(2**LOW_PREC_BITS*self.e0approx.sqrt())).round())/2**LOW_PREC_BITS
        self.e1=self.v1(self.e)
        #this is a unit of K (we have to check that it is not the fundamental unit of F)

        self.embedding=embedding
        if self.level == 1:
            self.tau_0,self.gtau = tau0_and_cusp(self.base, self.ext, self.e, self.v0, embedding, extra_parameters = extra_parameters)
        else:
            self.tau_0,self.gtau = tau0_and_cusp_withlevel(self.base, self.ext, self.v0,self.level,find_all = True)


        r=self.base.ring_of_integers().ring_generators()[0]
        self._changebasismatrix=Matrix(QQbar,2,2,[1,self.v0(r),1,self.v1(r)]).inverse()
        self._changebasismatrixapprox=Matrix(RealField(prec_bits),2,2,[self._changebasismatrix[ii,jj] for ii in range(2) for jj in range(2)])

        d_temp=self.base.different().gens_reduced()[0]
        #we have to take a totally positive generator
        if (self.v0(d_temp)>0) and (self.v1(d_temp)>0):
            self.d=d_temp
        if (self.v0(d_temp)>0) and (self.v1(d_temp)<0):
            self.d=d_temp*self.e
        if (self.v0(d_temp)<0) and (self.v1(d_temp)>0):
            self.d=-(d_temp)*self.e
        if (self.v0(d_temp)<0) and (self.v1(d_temp)<0):
            self.d=-d_temp
        self.d0=self.v0(self.d)
        self.d1=self.v1(self.d)

        self.estimated_C=self.estimate_C(1/10).ceil()
        if C is None:
            self.C=self.estimated_C
        else:
            self.C=C
        self.ddelta=self.epssqrt/(self.C*(1+self.e0approx))
        print('Estimated delta = ',RR(self.ddelta))



    def how_many_coeffs(self,r):
        return ceil(RealField()((self.base.discriminant()*r**2)/(30*self.ddelta**2)))

    #checks if the field F(sqrt{alpha}) is ATR
    def is_ATR(self,alpha=None):
        F=self.base
        if alpha is None:
            alpha=self.alpha
        Q.<y>=PolynomialRing(F)
        y=Q.gen()
        if not (y^2-alpha).is_irreducible():
            return False
        K=NumberField(y^2-alpha,'beta')
        if len(K.real_embeddings())==2:
            return True
        return False

    def estimate_number_limits(self, original_volume):
        return ZZ((original_volume/(acosh(RR(1+1/(200*.9*.9)))**2)).ceil())

    def prepare_limits(self,deltas=None,lims_from_disk=None,rng=None,out_limits=None,cusp_cf=None,max_length = 0,threshold = None, assume_E1 = False, level = 1):
        V=None
        #assert(len(deltas)==2)
        if not lims_from_disk is None:
            V=[]
            try:
                for lim in lims_from_disk:
                    V.extend([Limit(*l[0],conj=l[1]) for l in load(lim)])
            except IOError: pass
        if V is None:
            if threshold is None:
                threshold = 10**(-3)
            if level == 1:
                V,img_part = quotients_and_limits(self.base, self.tau_0, self.v0, self.v0approx, self.v1,cusp_cf, max_length,threshold,assume_E1,level = level)
            else:
                if isinstance(self.tau_0,list):
                    opt_V = None
                    opt_img = 0
                    for ii in range(len(self.tau_0)):
                        tau0 = self.tau_0[ii]
                        gtau = self.gtau[ii]
                        V,img_part = quotients_and_limits_withlevel(self.base,tau0,self.v0,self.v0approx,self.v1,gtau,level)
                        if img_part > opt_img:
                            opt_V = V
                            opt_img = img_part
                    V,img_part = opt_V,opt_img
                else:
                    V,img_part = quotients_and_limits_withlevel(self.base,self.tau_0,self.v0,self.v0approx,self.v1,self.gtau,level)
                e = self.e
                '''
                badmatrix = matrix(self.base,2,2,[1,0,3*level,1])
                lim = V[3]
                lower = (lim.x0,lim.x1)
                upper = (lim.y0,lim.y1)
                l1 = lower
                l2 = (lim.x0,self.I/self.v1(level).sqrt())
                u1 = (lim.y0,self.I/self.v1(level).sqrt())
                u2 = upper

                newl1 = self.act(badmatrix,l1)
                newu1 = self.act(badmatrix,u1)
                V[3] = Limit(newl1[0],newl1[1],newu1[0],newu1[1])
                V.append(Limit(l2[0],l2[1],u2[0],u2[1]))
                '''
            if out_limits is None:
                mystr='/tmp/limits_initial_%s_%s.sobj'%(self.base.discriminant(),self.ext.relative_discriminant().norm())
            else:
                mystr=out_limits
            save([(v.unpack(),v.conj) for v in V],mystr)
            print('Saved %s limits to %s.'%(len(V),mystr))

        if len(V)==0:
            raise RuntimeError
        print('Number of limits=%s'%len(V))
        print('Original Min Imag part:')
        min_imag_part_original=min([l.imaginary_part_approx() for l in V])
        print("%s,%s"%(min_imag_part_original,img_part))


        self.VV=V
        self.WW=[]
        if rng is None:
            rng=range(len(V))
        ww_info=0
        for ii in rng:
            new_WW,ww_info=self.break_limits([V[ii]],deltas,ii,ww_info)
            self.WW.append(new_WW)
        print('Original Min Imag part:')
        print(min_imag_part_original)
        limits=[]
        for filename in self.WW:
        	limits.extend([Limit(*l[0],conj=l[1]) for l in load(filename)])
        print('Min Imag part after prepare_limits:')
        tmp=min([oo]+[l.imaginary_part_approx() for l in limits])
        print(tmp)
        return tmp

    #this function computes the ATR point
    def ATR_point(self,f,m0,m1,preci=LOW_PREC_BITS,limits=None,parallel=True,use_quads=False,n_cpus=None,expand_factor=512):
        Cf=ComplexField(2*preci)
        if limits is None:
            limits=[]
            for filename in self.WW:
                limits.extend([Limit(*l[0],conj=l[1]) for l in load(filename)])
        #print('Min Imag part:')
        #print(min([l.imaginary_part_approx() for l in limits]))

        print('Starting integration...')
        gc.disable()
        swap=self.swapped_embeddings
        Points=[]
        Npoints=0
        Points0=[]
        Points1=[]
        Points2=[]
        Points3=[]

        for P in limits:
            conjug=P.conj
            P=P.unpack()
            if P[1]!=infinity:
                Points.extend([(Cf(P[0]),Cf(P[1]),conjug,1),(Cf(P[2]),Cf(P[1]),conjug,-1)])
                Points0.extend([Cf(P[0]),Cf(P[2])])
                Points1.extend([Cf(P[1]),Cf(P[1])])
                Points2.extend([conjug,conjug])
                Points3.extend([1,-1])
                Npoints+=2
            if P[3]!=infinity:
                Points.extend([(Cf(P[0]),Cf(P[3]),conjug,-1),(Cf(P[2]),Cf(P[3]),conjug,1)])
                Points0.extend([Cf(P[0]),Cf(P[2])])
                Points1.extend([Cf(P[3]),Cf(P[3])])
                Points2.extend([conjug,conjug])
                Points3.extend([-1,1])
                Npoints+=2
        t=walltime()
        if parallel:
            if n_cpus is None:
                n_cpus=sage.parallel.ncpus.ncpus()
            print('n_cpus=%s'%n_cpus)
            N=ceil(len(Points)/(n_cpus*expand_factor))
            M=ceil(QQ(m1-m0)/(n_cpus*expand_factor))
            inputs=[]
            for ii in range(n_cpus*expand_factor):
                # Each process deals with some of the points (and all the coefficients)
                # inputs.append((Points[N*ii:N*(ii+1)],self,f,nloop,preci,swap,False))
                # Each process deals with some of the coefficients (and all the points)
                inputs.append((Points0,Points1,Points2,Points3,self,f,m0+M*ii,m0+min([m1-m0,M*(ii+1)]),preci,swap, False,Npoints))
                # inputs.append((Points,self,f,nloop[M*ii:M*(ii+1)],preci,swap,False, Npoints))
            total=Cf(0)
            print('len(inputs)=',len(inputs))
            I=py_integrate(inputs)
            ii=0
            out_vec=[]
            for res in I:
                ii+=1
                print('Done with input %s/%s'%(ii,len(inputs)))
                print(res[1])
                out_vec.append(Cf((res[1])))
            total=sum(sorted(out_vec,key=lambda x:x.abs(),reverse=True))
        else:
            if use_quads==True:
                assert 0
                # total=integrate_quads(Points,self,f,nloop,preci,swap)
            else:
                if preci==53:
                    total=Cf(integrate_doubles(Points0,Points1,Points2,Points3,self,f,m0,m1,preci,swap, Npoints))
                else:
                    total=Cf(integrate_hiprec(Points0,Points1,Points2,Points3,self,f,m0,m1,preci,swap, Npoints))
        gc.enable()
        print('Total (wall) time = %s'%walltime(t))
        return total

    def delta(self,tolerance=1,C=None):
        if C is None:
            return tolerance*self.ddelta
        else:
            return self.epssqrt/(C*(1+self.e0approx))

    def fundom_rep(self,x):
        A=self._changebasismatrixapprox
        return ((A[0,0]*x[0]+A[0,1]*x[1]).floor(),(A[1,0]*x[0]+A[1,1]*x[1]).floor())

    # This is an implementation of Lemma 3.6_1 of Freitag
    # Given x in R^2, and a positive epsilon, finds
    # integers c,d in O with c nonzero and such that
    # ||cx+d|| <= epsilon, ||c|| <= C/epsilon
    # with some C which is only depending on the number field.
    def approximate(self,Z,epsilon=None):
        initial_time=cputime()
        if epsilon is None:
            epsilon=self.find_epsilon(z0,z1)
        x=(Z[0].real(),Z[1].real())
        y=(Z[0].imag(),Z[1].imag())
        e=self.e

        x0approx,x1approx=RealField(prec_bits)(x[0]),RealField(prec_bits)(x[1])
        y0approx,y1approx=RealField(prec_bits)(y[0]),RealField(prec_bits)(y[1])
        y0sq,y1sq=y0approx**2,y1approx**2

        M=(1/epsilon).ceil()
        w=self.base.maximal_order().ring_generators()[0]
        matches=dict([])
        for c1,nrm in self.SmallInts():
            if cputime(initial_time)>TIME_LIMIT:
                print('Reached TIME_LIMIT of %s seconds...'%(TIME_LIMIT))
                raise RuntimeError
            cx=(self.v0approx(c1)*x0approx,self.v1approx(c1)*x1approx)
            dtmp0,dtmp1=self.fundom_rep(cx)
            d1=-dtmp0-dtmp1*w
            # So cx+d1 belongs to the fundamental domain
            cxpd0,cxpd1=cx[0]+self.v0approx(d1),cx[1]+self.v1approx(d1)
            r=((M*cxpd0).round(),(M*cxpd1).round())
            try:
                v=matches[r]
                assert len(v)<100
                for c0,d0 in v:
                    c,d=c0-c1,d0-d1
                    a,mb=bezout(self.base,d,c)
                    g=mb*c+a*d
                    if g.norm().abs()==1:
                        c0approx,c1approx=self.v0approx(c),self.v1approx(c)
                        d0approx,d1approx=self.v0approx(d),self.v1approx(d)
                        cxpd0=c0approx*x0approx+d0approx
                        cxpd1=c1approx*x1approx+d1approx
                        if (cxpd0**2+y0sq*c0approx**2)*(cxpd1**2+y1sq*c1approx**2)<1:
                            return c,d #,epsilon*max([abs(self.v0approx(c)),abs(self.v1approx(c))])
                    else:
                        c2=c/g
                        d2=d/g
                        min_denom=1
                        min_n=None
                        for n in range(-5,6):
                            c=self.e**n*c2
                            d=self.e**n*d2
                            c0approx,c1approx=self.v0approx(c),self.v1approx(c)
                            d0approx,d1approx=self.v0approx(d),self.v1approx(d)
                            cxpd0=c0approx*x0approx+d0approx
                            cxpd1=c1approx*x1approx+d1approx
                            denom=(cxpd0**2+y0sq*c0approx**2)*(cxpd1**2+y1sq*c1approx**2)
                            if denom<min_denom:
                                min_denom=denom
                                min_n=n
                        if min_denom<1:
                            c=self.e**min_n*c2
                            d=self.e**min_n*d2
                            return c,d
                matches[r].append((c1,d1))
            except KeyError:
                matches[r]=[(c1,d1)]
        return None

    # Acts by the diagonal matrix e,e^-1 to bring the
    # imaginary parts together
    def equalize_imags(self,Z):
        # assert Z[0].imag()>0 and Z[1].imag()>0
        e00=self.e0approx
        if e00>1:
            em,e0m,sign=-self.e,-self.e0,1
        else:
            e00=1/e00
            em,e0m,sign=-1/self.e,-1/self.e0,-1

        assert(self.v0approx(self.e)>0)
        x0,x1=Z[0],Z[1]
        n=0
        Rf=RealField(prec_bits)
        while Rf(x0.imag())*e00<Rf(x1.imag()):
            x0=(x0*e0m).conjugate()
            x1/=-e0m
            n+=1
        while Rf(x0.imag())>Rf(x1.imag())*e00:
            x0=(x0/e0m).conjugate()
            x1*=-e0m
            n-=1
        return sign*n,[x0,x1],(-1)**(n%2)


    def act(self,M,Z):
        tup=(act_single(M,Z[0],self.v0),act_single(M,Z[1],self.v1))
        return tup

    # Find an epsilon to approximate. Ensures that at least the imaginary
    # parts are multiplied by lambda.
    def find_epsilon(self,xx,lam=2):
        # print('xx=',xx)

        Rf=RealField(prec_bits)
        tmp=(2**LOW_PREC_BITS*((self.C**2*Rf(xx[0].imag())*Rf(xx[1].imag()))**(.25))).round()/2**LOW_PREC_BITS
        # print('Would like eps =',RealField()(tmp))
        if tmp<1/50:
            A=self.C**2*(xx[0].imag()**2+xx[1].imag()**2)
            B=(self.C**2*xx[0].imag()*xx[1].imag())**2
            if A>=1 or B>=((1-A)/2)**2:
                raise RuntimeError
            tmp=((100*(((1/lam-A)+RealField()(((1/lam-A)**2-4*B)).sqrt())/2)**(1/4)).floor())/100
            if tmp<1/100:
                raise RuntimeError
        elif tmp>=1:
            raise RuntimeError
            #tmp=99/100
        # print('eps =',tmp)
        return tmp

    def move_to_siegel_domain(self,lim):
        #print('Entering move_to_siegel_domain...')
        Rf=RealField(prec_bits)
        n_iters=1
        curr_delta=lim.imaginary_part_approx()
        limconj=1
        g=Matrix(self.base,2,2,1)
        opt_delta=curr_delta
        opt_g=g
        opt_conj=limconj
        while n_iters<20 and curr_delta<self.delta():
            g,curr_delta,newconj=self.reduction(lim,G=g)
            limconj*=newconj
            if curr_delta>opt_delta:
                opt_delta=curr_delta
                opt_g=g
                opt_conj=limconj
            n_iters+=1
        if opt_delta<self.delta():
            print('Updating delta...')
            self.ddelta=min([self.ddelta,RationalField()(9/10)*RationalField()(((2**LOW_PREC_BITS*opt_delta).ceil())/2**LOW_PREC_BITS)])
            print('New delta =',RealField()(self.ddelta))
        #print('Done with move_to_siegel_domain.')
        return opt_g,opt_delta,opt_conj

    # Find a better representative for Z in H^2.
    # By that we mean one such that its imaginary part is reasonable
    def reduction(self,initial_lim,G=None):
        # print('Entering reduction...')
        initial=(initial_lim.x0,initial_lim.x1)
        O=self.base.ring_of_integers()
        w=O.ring_generators()[0]
        e=-self.e
        conjugate=1

        if G is None:
            G=Matrix(self.base,2,2,1)
            Z1=initial
        else:
            Z1=self.act(G,initial)
        n,Z,conj=self.equalize_imags(Z1)
        conjugate*=conj
        exxn=e**n
        MM=Matrix(self.base,2,2,[exxn,0,0,1])*G
        try:
            fin0,fin1=self.act(MM,initial)
            bound=RealField(prec_bits)(fin0.imag()*fin1.imag()).sqrt()
            if bound>self.delta():
                return MM,bound,conjugate
            try:
                #print('Entering approximate...')
                c,d=self.approximate(Z,epsilon=1/2)
                #print('Done with approximate')
                MM=matrix_from_bottom_row(self.base,c,d)*MM

                tmp=self.act(MM,initial)
                n,Z,conj=self.equalize_imags(tmp)
                conjugate*=conj
                exxn=e**n
                MM=Matrix(self.base,2,2,[exxn,0,0,1])*MM
            except RuntimeError: pass
            fin0,fin1=self.act(MM,initial)
            # print('Done with reduction')
            bound=(RealField(prec_bits)(fin0.imag()*fin1.imag())).sqrt()
        except RuntimeError:
            print('Maybe too much recursion depth...')
            raise RuntimeError
        return MM,bound,conjugate

    def geodesic(self,P,Q):
        a=(P.real()**2-Q.real()**2+P.imag()**2-Q.imag()**2)/(2*(P.real()-Q.real()))
        r2=(P.real()-a)**2+P.imag()**2
        return a,r2

    # def intersect_with_boundary(self,P,Q,t0):
    #     t=RationalField()(((RealField(prec_bits)(2**LOW_PREC_BITS*t0)).round())/2**LOW_PREC_BITS)
    #     return Q.real()+t*self.I

    def intersect_with_boundary(self,P,Q,t0):
        t=RationalField()(((RealField(prec_bits)(2**LOW_PREC_BITS*t0)).round())/2**LOW_PREC_BITS)
        P0=RealField(prec_bits)(P.real())
        Q0=RealField(prec_bits)(Q.real())
    
        if abs(P0-Q0)<2**(-10):
            m=(P.real()-Q.real())/(P.imag()-Q.imag())
            return P.real()+m*(t-P.imag())+t*self.I
        a,r2=self.geodesic(P,Q)
        tmp1=(r2-t**2)
        tmp=RationalField()(((2**LOW_PREC_BITS*((RealField(prec_bits)(tmp1)).sqrt())).round())/(2**LOW_PREC_BITS))
        if P0<Q0:
            realparts=[P.real(),Q.real()]
            r0=P0
            r1=Q0
        else:
            realparts=[Q.real(),P.real()]
            r0=Q0
            r1=P0
    
        x1=a+tmp
        x10=RealField(prec_bits)(x1)
        if r0<x10 and x10<r1:
            return x1+t*self.I
        else:
            return a-tmp+t*self.I

    def break_limits(self,vec,deltas,ii_info=0,ww_info=0):
        V=[]
        W=[]
        filename='/tmp/limits_%s_%s_%s_%s.sobj'%(self.base.discriminant(),self.ext.relative_discriminant().norm(),ii_info,self.embedding)
        # try:
        #     if os.path.isfile(filename):
        #         print('File corresponding to limit %s exists...skipping.'%(ii_info))
        #         tmp=load(filename)
        #         return filename,ww_info+len(tmp)
        # except IOError: pass
        if deltas is None:
            tolerance=RationalField()(9/10)
            d00=self.delta(C=self.estimated_C)*tolerance**2
            d11=self.delta(C=self.estimated_C)*tolerance
        else:
            d00=deltas[0]
            d11=deltas[1]
        assert d00<d11

        e00=self.e0approx
        if e00>1:
            em,e0m,sign=-self.e,-self.e0,1
        else:
            e00=1/e00
            em,e0m,sign=-1/self.e,-1/self.e0,-1
        e0mapprox=self.v0approx(em)

        part=RationalField()(100)/len(vec)
        for l in vec:
            if self.level != 1 or l.imaginary_part_approx()>d00:
                print("Don't break this limit...")
                W.append(l)
            else:
                print('We break it!')
                x0,x1,y0,y1=l.unpack()
                n=ZZ((max(map(lambda x:((((-RealField(prec_bits)(x.imag()).log()))/((-e0mapprox).log()))).round(),[x0,y0]))))
                e0xxn=(-e0mapprox)**n
                V.append((Limit(x0,x1,y0,self.I*e0xxn,conj=l.conj),part))
                if n%2==0:
                	W.append(Limit(x0*e0xxn,self.I,y0*e0xxn,y1,conj=l.conj))
                else:
                	W.append(Limit((-x0*e0xxn).conjugate(),self.I,(-y0*e0xxn).conjugate(),y1,conj=-l.conj))
        t=None
        flag=False
        ii=0
        progress=RationalField()(0)

        # original_volume = sum([l.volume() for l in V])
        # print('Original Volume:')
        # print("%s"%(original_volume))
        # print('Estimated number of limits')
        # print("%s"%self.estimate_number_limits(original_volume))

        while len(V)>0:
            ii+=1
            print('len(V)=%s, len(W)=%s, deltaAG=%s, done_limits=%s, progress=%s'%(len(V),2*len(W)+ww_info,self.minW,ii_info,RealField()(progress)))
            lim,part=V.pop()
            x0,x1,y0,y1=lim.unpack()
            limconj=lim.conj
            g,newt,conj=self.move_to_siegel_domain(lim)
            conj*=limconj
            x0,x1=self.act(g,(x0,x1))
            y0,y1=self.act(g,(y0,y1))
            x0app,x1app,y0app=ComplexField(LOWER_PREC_BITS)(x0),ComplexField(LOWER_PREC_BITS)(x1),ComplexField(LOWER_PREC_BITS)(y0)
            if y1 is Infinity:
                y1app=Infinity
            else:
                y1app=ComplexField(LOWER_PREC_BITS)(y1)
            newV=[]
            newW=[]
            d0sq=d00**2
            d1sq=d11**2
            if y1app is not Infinity and x0app.imag()*y1app.imag()<d1sq:
                t1=self.intersect_with_boundary(x1,y1,d1sq/x0app.imag())
                #print('t1img=',RealField()(t1.imag()))
                print('A', end=' ')
                t1app=ComplexField(LOWER_PREC_BITS)(t1)
                if y0app.imag()*t1app.imag()<d0sq:
                    print('1')
                    t0=self.intersect_with_boundary(x0,y0,d0sq/t1app.imag())
                    #print('t0img=',RealField()(t0.imag()))
                    newV=[(Limit(t0,x1,y0,t1,conj),part/3),(Limit(x0,t1,y0,y1,conj),part/3)]
                    newW=[Limit(x0,x1,t0,t1,conj)]
                    progress+=part/3
                else:
                    print('2')
                    newV=[(Limit(x0,t1,y0,y1,conj),part/2)]
                    newW=[Limit(x0,x1,y0,t1,conj)]
                    progress+=part/2
            elif x1app.imag()*y0app.imag()<d1sq:
                t0=self.intersect_with_boundary(x0,y0,d1sq/x1app.imag())
                #print('t0img=',RealField()(t0.imag()))
                print('B', end=' ')
                t0app=ComplexField(LOWER_PREC_BITS)(t0)
                if y1app is not Infinity and y1app.imag()*t0app.imag()<d0sq:
                    print('1')
                    t1=self.intersect_with_boundary(x1,y1,d0sq/t0app.imag())
                    #print('t1img=',RealField()(t1.imag()))
                    newV=[(Limit(t0,x1,y0,y1,conj),part/3),(Limit(x0,t1,t0,y1,conj),part/3)]
                    newW=[Limit(x0,x1,t0,t1,conj)]
                    progress+=part/3
                else:
                    print('2')
                    newV=[(Limit(t0,x1,y0,y1,conj),part/2)]
                    newW=[Limit(x0,x1,t0,y1,conj)]
                    progress+=part/2
            elif y1app is not Infinity and y0app.imag()*y1app.imag()<d0sq:
                t0=self.intersect_with_boundary(x0,y0,d0sq/y1app.imag())
                #print('t0img=',RealField()(t0.imag()))
                print('C', end=' ')
                t0app=ComplexField(LOWER_PREC_BITS)(t0)
                if t0app.imag()*x1app.imag()<d1sq:
                    print('1')
                    newV=[(Limit(t0,x1,y0,y1,conj),part/2),(Limit(x0,x1,t0,y1,conj),part/2)]
                    newW=[]
                    progress+=0
                else:
                    print('2')
                    newV=[(Limit(t0,x1,y0,y1,conj),part/2)]
                    newW=[Limit(x0,x1,t0,y1,conj)]
                    progress+=part/2
            else:
                print('D')
                newV=[]
                newW=[Limit(x0,x1,y0,y1,conj)]
                progress+=part
            V.extend(newV)
            self.minW=min([self.minW]+[RDF(ww.imaginary_part_approx()) for ww in newW])
            W.extend(newW)
        Wconj=[]
        for l in W:
            #limits of the first integral
            xt=(self.e0*l.x0,self.e1*l.x1.conjugate())
            if l.y1 is Infinity:
                Wconj.append(Limit(self.e0*l.x0,self.e1*l.x1.conjugate(),self.e0*l.y0,Infinity,conj=l.conj))
            else:
                Wconj.append(Limit(self.e0*l.x0,self.e1*l.x1.conjugate(),self.e0*l.y0,self.e1*l.y1.conjugate(),conj=l.conj))
        W.extend(Wconj)

        WW=[(w.unpack(),w.conj) for w in W]
        save(WW,filename)
        return filename,len(WW)+ww_info

    def update_cached_ints(self):
        r=self.base.ring_of_integers().ring_generators()[0]
        D=self.base.discriminant()
        sqrtD=RealField()(D).sqrt()
        Nold=self.cached_ints_N
        N=Nold+D
        for a in range(2*Nold):
            amod2=a%2
            a2=QQ(a/2)
            for b in filter(lambda x:x%2==amod2,range((2*Nold/sqrtD).ceil(),(2*N/sqrtD).ceil())):
                b2=QQ(b/2)
                nrm=(a2+b2*sqrtD).floor()
                #c1=a2-b2
                #d1=b*r
                try:
                    self.cached_ints[nrm].append((a2,b2))
                except KeyError,AttributeError:
                    self.cached_ints[nrm]=[(a2,b2)]

        for a in range(2*Nold,2*N):
            amod2=a%2
            a2=QQ(a/2)
            for b in filter(lambda x:x%2==amod2,range((2*N/sqrtD).ceil())):
                b2=QQ(b/2)
                nrm=(a2+b2*sqrtD).floor()
                #c1=a2-b2
                #d1=b*r
                try:
                    self.cached_ints[nrm].append((a2,b2))
                except KeyError,AttributeError:
                    self.cached_ints[nrm]=[(a2,b2)]
        self.cached_ints_N=N

    def SmallInts(self):
        r=self.base.ring_of_integers().ring_generators()[0]
        mycache=self.cached_ints
        ii=0
        while True:
            if ii>=self.cached_ints_N:
                # print('updating cache...')
                self.update_cached_ints()
                # print('done!')
            try:
                for a,b in mycache[ii]:
                    for aa,bb in list(set([(a,b),(-a,b),(a,-b),(-a,-b)])):
                        yield aa-bb+2*bb*r,ii
            except KeyError: pass
            ii+=1

    def estimate_C(self,eps):
        r=self.base.ring_of_integers().ring_generators()[0]
        w0=self.v0approx(r)
        w1=self.v0approx(r)
        B=self.base.discriminant().sqrt().round()
        Np=(2*(w0**2+w1**2)).sqrt()
        return 2*eps*max([x[1] for x in islice(self.SmallInts(),ZZ((Np/(eps**2)).ceil()))])

@parallel('fork')
def py_integrate(Points0, Points1,Points2,Points3, atrfield,aprs,m0,m1,preci,swap,use_quads, Npoints):
    if preci != 53 and use_quads==False:
        return integrate_hiprec(Points0, Points1,Points2,Points3,atrfield,aprs,m0,m1,preci, swap, Npoints)
    if use_quads==True:
        return integrate_quads(Points0, Points1,Points2,Points3,atrfield,aprs,m0,m1,preci, swap, Npoints)
    else:
        return integrate_doubles(Points0, Points1 ,Points2,Points3,atrfield,aprs,m0,m1,preci, swap, Npoints)


@fork
def compute_coeffs_p(F,max_norm,p_range,Ainv,level,v0,v1,e2,single):
    def equalize(x):
        x0,x1=v0(x),v1(x)
        if not (x0>0 and x1>0):
            x=-x
            x0,x1=-x0,-x1
        e0=v0(e2).abs()
        while x0*e0<x1:
            x*=e2
            x0,x1=x0*e0,x1/e0
        while x0>x1*e0:
            x/=e2
            x0,x1=x0/e0,x1*e0
        return x
    if single:
        new_aprs=[]
        for indx,p in p_range:
            eltlist=[equalize(xx) for xx in F.elements_of_norm(p)]
            new_aprs.append((indx,[(x,get_ap(F,Ainv,x),v0(x),v1(x)) for x in eltlist],p))
        return new_aprs
    else:
        return compute_coeffs_p_cython(F,max_norm,p_range,Ainv,level,v0,v1,e2)

@parallel(ncpus=NCPUS)
def compute_coefficients(F,min_norm, max_norm,max_norm_powers, i0, Ainv, level, v0, v1, e2, prec,single):
    aprs=[]
    pending_primes=[]
    ii=0
    ind=i0
    for p in prime_range(min_norm,max_norm):
        if single and kronecker_symbol(F.discriminant(),p) == -1:
            continue
        pending_primes.append((i0+ii,p))
        ii+=1
        if ii%1000==0:
            aprs.extend(compute_coeffs_p(F, max_norm_powers, pending_primes,Ainv,level,v0,v1,e2,single))
            pending_primes=[]
    aprs.extend(compute_coeffs_p(F,max_norm_powers, pending_primes,Ainv,level,v0,v1,e2,single))
    return aprs

class TestDL:
    def __init__(self,E,prec=14,working_prec=53,level = 1, hmf = None):
        self.working_prec=working_prec
        self.prec=prec
        self.E=E
        self.F=E.base_ring()
        self.Jtau=None
        self.hmf=hmf
        self.max_norm=[0,0]
        self.level = level

    def calculate_HMF(self,max_norm,parallel = True):
        print('Creating HilbertModularForm...')
        self.hmf=HilbertModularForm(self.E,max_norm,level = self.level, expand_factor=8, n_cpus = None if parallel == True else 1)

    def calculate_limits(self,K,max_length = None, assume_E1 = False): # Defaults to 5
        delta=K.prepare_limits(max_length=max_length, assume_E1 = assume_E1,level = self.level)
        max_norm = ceil(RealField()((K.base.discriminant()*self.prec**2)/(30*delta**2)))
        print('Updated max_norm to %s'%(RR(max_norm)))
        return max_norm

    def calculate_Jtau(self,K,max_norm,parallel=True):
        F=self.F
        print('Computing up to max_norm =',RR(max_norm))
        if self.hmf is None or self.hmf.max_norm_aprs<max_norm:
            print("Need more coefficients for the HMF. We'll calculate them now...")
            self.calculate_HMF(max_norm,parallel)
        return K.ATR_point(self.hmf,1,max_norm,self.working_prec,parallel=parallel,use_quads=False,expand_factor=8)

    def check_conjecture(self,x_coord,max_length=5,max_norm=None,alpha = None, Jtau = None, parallel = True,prec=15, force_embedding = None, assume_E1 = False, integrate = True, extra_parameters = None):
        # Initialize the ATR field
        Fx=PolynomialRing(self.F,names='Y')
        Y=Fx.gen()
        if alpha is None:
            self.alpha = find_alpha(self.E,x_coord)
        else:
            self.alpha = alpha
        print('alpha=',self.alpha)
        hK=self.F.extension(Y**2-self.alpha,names='beta').class_number()
        print('Class number =',hK)
        if hK>2:
            raise NotImplementedError, "Only deal with class number up to 2"
        if Jtau is None:
            self.Jtau=0
        else:
            self.Jtau=Jtau
        self.K=[]
        for h in range(hK):
            alpha = self.alpha
            if not force_embedding is None:
                K = ATRField(self.F,alpha,embedding=force_embedding, extra_parameters = extra_parameters,level = self.level)
            else:
                K = ATRField(self.F,alpha,embedding=h, extra_parameters = extra_parameters,level = self.level)
            self.K.append(K)
            if Jtau is None:
                # Calculate the limits
                print('Calculating limits...')
                this_max_norm=self.calculate_limits(K,max_length, assume_E1 = assume_E1)
                self.max_norm[h]=this_max_norm
                # print('Done')
                if max_norm is None:
                    max_norm=this_max_norm
                # Calculate the ATR point
                if integrate == True:
                    print('Calculating the ATR point...')
                    print('Jtau[%s]=%s'%(h,new_Jtau))
                else:
                    new_Jtau = 0
                self.Jtau+=new_Jtau

        #self.Jtau=Jtau

        # Check that the computed point is in the same line as the one given by x_coord
        beta=self.K[0].ext.gen()
        EK=self.E.base_extend(self.F.embeddings(self.K[0].ext)[0])
        P=EK.lift_x(x_coord)

        #these are the embeddings of K that we want to consider...v0K is complex and v1K is real. They extend K1.ext.v0 and K1.ext.v1 (the two embeddings of F)
        emb=self.K[0].ext.complex_embeddings(prec_bits)

        found = False
        for ii in range(len(emb)):
            if abs(emb[ii](self.F.gen())-self.K[0].v0(self.F.gen()))<10**(-10):
                found = True
                self.v0K=emb[ii]
        assert found

        found = False
        for ii in range(len(emb)):
            if abs(emb[ii](self.F.gen())-self.K[0].v1(self.F.gen()))<10**(-10):
                found = True
                self.v1K=emb[ii]
        assert found


        #the period lattice of EK with respect to these embeddings. It is the same as the period lattice of E with respect to K1.ext.v0 and K1.ext.v1, although for the lattice of EK sage returns the basis in the reversed order...
        self.LK0 = EK.period_lattice(self.v0K)
        self.LK1 = EK.period_lattice(self.v1K)

        B=self.LK0.basis(prec_bits)
        if abs(B[0].real())>10**(-10):
            lambda_0_plus_K=B[0].real()
        else:
            lambda_0_plus_K=B[1].real()
        if abs(B[0].imag())>10**(-10):
            lambda_0_minus_K=B[0].imag()*ComplexField(prec_bits).gen()
        else:
            lambda_0_minus_K=B[1].imag()*ComplexField(prec_bits).gen()

        B=self.LK1.basis(prec_bits)
        if abs(B[0].real())>10**(-10):
            lambda_1_plus_K=B[0].real()
        else:
            lambda_1_plus_K=B[1].real()
        if abs(B[0].imag())>10**(-10):
            lambda_1_minus_K=B[0].imag()*ComplexField(prec_bits).gen()
        else:
            lambda_1_minus_K=B[1].imag()*ComplexField(prec_bits).gen()

        self.lambda_0_minus_K=lambda_0_minus_K
        self.lambda_0_plus_K=lambda_0_plus_K
        self.lambda_1_minus_K=lambda_1_minus_K
        self.lambda_1_plus_K=lambda_1_plus_K
        print('lambda_0_minus_K=',lambda_0_minus_K)
        print('lambda_0_plus_K=',lambda_0_plus_K)
        print('lambda_1_minus_K=',lambda_1_minus_K)
        print('lambda_1_plus_K=',lambda_0_plus_K)

        #this is the complex number in C/Lambda_0 corresponding to the nontorsion point	
        self.z=self.LK0.elliptic_logarithm(P,prec = prec_bits,reduce = True)
        #a small verification that z is the correct value
        try:
            exp_tildeP=self.LK0.elliptic_exponential(self.z)
            print('Verifying that z is the correct value...')
            print('Checking that the real part of z is indeed torsion. This should be a small linear relation:')
            print(gp.lindep([self.z.real(),self.lambda_0_plus_K],25))
        except ArithmeticError:
            print('Could not check that z was OK...continuing nevertheless')
        self.check_lindeps(prec=prec)

    def check_lindeps(self,prec,Jtau = None,z = None):
        if Jtau is None:
            Jtau=self.Jtau
        if z is None:
            z=self.z
        #finally the verification
        print('Verification of the imaginary part of the point')
        print(gp.lindep([(Jtau/self.lambda_1_plus_K).imag(),z.imag(),self.lambda_0_minus_K.imag()],prec))
        print('Verification of the real part of the point (should be torsion)')
        print(gp.lindep([(Jtau/self.lambda_1_plus_K).real(),self.lambda_0_plus_K],prec))
