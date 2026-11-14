import sys, numpy as np, netCDF4 as nc
def gm(x,lat): 
    w=np.cos(np.radians(lat))[:,None]*np.ones(x.shape[-2:]); return float((x*w).sum()/w.sum())
def load(f):
    d=nc.Dataset(f); lat=d['lat'][:]; lon=d['lon'][:]
    v={k:np.array(d[k][:]).mean(0) for k in d.variables if d[k].ndim==3}
    return lat,lon,v
def summary(f,sslon=180.):
    lat,lon,v=load(f); LON,LAT=np.meshgrid(lon,lat)
    mu=np.cos(np.radians(LAT))*np.cos(np.radians(LON-sslon)); w=np.cos(np.radians(LAT))
    dw=np.where(mu>1e-6,1.0,np.where(mu<-1e-6,0.0,0.5)); r={}   # terminator columns count half to each side
    r['ts']=gm(v['ts'],lat); r['ts_day']=float((v['ts']*w*dw).sum()/(w*dw).sum()); r['ts_night']=float((v['ts']*w*(1-dw)).sum()/(w*(1-dw)).sum())
    r['ts_max']=float(v['ts'].max()); r['ts_min']=float(v['ts'].min())
    r['toa']=gm(v['ntr'],lat); r['sfc']=gm(v['hfns'],lat)
    r['rst']=gm(v['rst'],lat); r['rss']=gm(v['rss'],lat); r['olr']=-gm(v['rlut'],lat)
    r['sic']=gm(v['sic'],lat); r['sit']=gm(v['sit'],lat); r['pr']=gm(v['pr'],lat)*86400*360*1000
    r['clt']=gm(v['clt'],lat); r['prw']=gm(v['prw'],lat)
    openw=(v['sic']<0.5); r['open_frac']=float((w*openw).sum()/w.sum())
    return r,(lat,lon,v,mu)
if __name__=='__main__':
    for f in sys.argv[1:]:
        r,_=summary(f); print(f.split('/')[-2],f[-8:-3],' '.join('%s=%.2f'%(k,x) for k,x in r.items()))
