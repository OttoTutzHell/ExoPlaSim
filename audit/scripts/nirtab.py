"""Offline replica of mknirtab (radmod.f90): stellar re-weighting of near-IR H2O absorption.
Usage: nirtab.py [TEFF ...]   (blackbody temperatures in K; any values)"""
import numpy as np, sys
xk=np.array([0.0010,0.0133,0.0422,0.1334,0.4217,1.334,5.623,31.62,177.8,1000.0])
hk=np.array([.20673,.08236,.01074,.03497,.01157,.00360,.03011,.01133,.00411,.02260,.01143,.00421,
             .01336,.01240,.00389,.00696,.01258,.00326,.00441,.01381,.00499,.00115,.00650,.00465,
             .00026,.00244,.00245,.00000,.00094,.00145]).reshape(10,3).T   # hk[band,k]
edges=np.array([0.70e-6,1.22e-6,2.27e-6,10e-6])
def planck(w,T): return w**-5/np.expm1(0.0143877735383/(w*T))
def fractions(T):
    w=np.logspace(np.log10(316.036116751e-9),-4,20000); b=planck(w,T)
    seg=0.5*(b[1:]+b[:-1])*np.diff(w); mid=0.5*(w[1:]+w[:-1]); tot=seg.sum()
    return np.array([seg[(mid>=edges[i])&(mid<edges[i+1])].sum()/tot for i in range(3)]), seg[mid<0.75e-6].sum()/tot
def alh(y): return 2.9*y/((1+141.5*y)**0.635+5.925*y)
def achou(u,r): return float(((r[:,None]*hk)*(1-np.exp(-xk[None,:]*u))).sum())
def astar(y,r):
    lo,hi=-12.,6.
    for _ in range(60):
        m=0.5*(lo+hi)
        if achou(10**m,np.ones(3))<alh(y): lo=m
        else: hi=m
    return achou(10**(0.5*(lo+hi)),r)
if __name__=="__main__":
    fg,z1g=fractions(5772.)
    print("solar band fractions, blackbody 5772 K:",fg.round(4)," Chou table sums:",hk.sum(1).round(4))
    temps=[float(a) for a in sys.argv[1:]] or [5772.,3800.,3042.,2825.,2600.,2566.]
    for T in temps:
        f,z1=fractions(T); r=f/fg
        print("T=%5.0f  frac<0.75um %.3f  band ratios %s  A(y)/A_LH(y) at y=0.01,0.1,1,10 cm: %s   A(1cm)=%.3f, as share of band 2: %.3f"%(
            T,z1,r.round(2),' '.join('%.2f'%(astar(y,r)/alh(y)) for y in (0.01,0.1,1,10)),astar(1,r),astar(1,r)/(1-z1)))
