import numpy as np, sys
from scipy.special import sph_harm_y
B=sys.argv[1].rstrip('/')+'/'   # directory holding <NAME>_T21_surf_0129.sra etc.; NAME and gravity below
def sra(f,nlat,nlon):
    v=[]
    with open(f) as fh:
        fh.readline()
        for l in fh: v+=[float(x) for x in l.split()]
    return np.array(v).reshape(nlat,nlon)
def run(res,N,nlat,nlon):
    topo=sra(B+'Ki_%s_surf_0129.sra'%res,nlat,nlon)/7.99
    mask=sra(B+'Ki_%s_surf_0172.sra'%res,nlat,nlon)
    x,w=np.polynomial.legendre.leggauss(nlat); x=x[::-1]; w=w[::-1]   # north first
    theta=np.arccos(x); phi=2*np.pi*np.arange(nlon)/nlon
    TH,PH=np.meshgrid(theta,phi,indexing='ij')
    W=w[:,None]*np.ones((1,nlon))*(2*np.pi/nlon)
    coef={}
    for n in range(N+1):
        for m in range(-n,n+1):
            Y=sph_harm_y(n,m,TH,PH); coef[(n,m)]=(topo*np.conj(Y)*W).sum(),Y
    def synth(filt):
        out=np.zeros_like(topo,dtype=complex)
        for (n,m),(c,Y) in coef.items(): out+=filt(n)*c*Y
        return out.real
    sea=mask<0.5; aw=W
    def stats(name,z):
        s=z[sea]; a=aw[sea]
        print('%s %-34s peak %5.0f m | sea rms %4.0f min %5.0f max %5.0f | sea area |z|>200m %4.1f%% | land rms err %4.0f'%(res,name,z.max(),np.sqrt((s*s*a).sum()/a.sum()),s.min(),s.max(),100*(a*(abs(s)>200)).sum()/a.sum(),np.sqrt((((z-topo)**2*aw)[~sea]).sum()/aw[~sea].sum())))
        return z
    print(res,'raw cell-mean peak %.0f m'%topo.max())
    stats('truncation only',synth(lambda n:1.0))
    e=lambda n:np.exp(-8*(n/N)**8)
    stats('exp filter once (dynamics sees)',synth(e))
    stats('exp filter twice (output zs)',synth(lambda n:e(n)**2))
    stats('Lanczos sigma',synth(lambda n:np.sinc(n/(N+1))))
    stats('Cesaro',synth(lambda n:1-n/(N+1)))
    return
run('T21',21,32,64); run('T42',42,64,128)
