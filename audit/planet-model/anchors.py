"""Summarise ExoPlaSim Ki-Engir runs as anchors for the planet model, and export two fields."""
import sys, json, glob, numpy as np, netCDF4 as nc
RUNS=[("Teegarden's Star",3034,7.22e-4,0.097,6.097,"tg_stock",""),("Teegarden's Star",3034,7.22e-4,0.097,5.677,"tg_f1200",""),("Teegarden's Star",3034,7.22e-4,0.097,5.346,"tg_f1300","still warming"),
      ("TRAPPIST-1",2566,5.53e-4,0.0898,5.539,"tr_f1000",""),("TRAPPIST-1",2566,5.53e-4,0.0898,5.156,"tr_f1100",""),("TRAPPIST-1",2566,5.53e-4,0.0898,5.0,"tr_p500",""),("TRAPPIST-1",2566,5.53e-4,0.0898,4.831,"tr_f1200","still warming")]
out=[]; fields={}
for star,T,L,M,per,run,note in RUNS:
    fs=sorted(glob.glob('runs/%s/MOST.0*.nc'%run))[-5:]; acc={}
    for f in fs:
        d=nc.Dataset(f); lat=np.array(d['lat'][:]); lon=np.array(d['lon'][:])
        for k in ('ts','tas','pr','rst','sic','lsm'): acc[k]=acc.get(k,0)+np.array(d[k][:]).mean(0)/len(fs)
    LON,LAT=np.meshgrid(lon,lat); w=np.cos(np.radians(LAT)); mu=np.cos(np.radians(LAT))*np.cos(np.radians(LON)); el=np.degrees(np.arcsin(np.clip(mu,-1,1)))
    day=np.where(mu>1e-6,1.0,np.where(mu<-1e-6,0.0,0.5)); g=lambda x,m=None:(float((x*w).sum()/w.sum()) if m is None else float((x*w*m).sum()/(w*m).sum()))
    land=acc['lsm']>0.5; a=(M*(per/365.25)**2)**(1/3); flux=1361*L/a**2
    e=dict(star=star,starT=T,L=L,M=M,period=per,flux=flux,G=g(acc['ts']),D=g(acc['ts'],day),N=g(acc['ts'],1-day),maxC=float(acc['ts'].max()-273.15),
           alb=1-g(acc['rst'])/(flux/4),asr=g(acc['rst']),iceFree=g((acc['sic']<0.5)*1.0,(~land)*1.0),
           l1=g(acc['ts'],((el>=15)&(el<30)&land)*1.0)-273.15,l2=g(acc['ts'],((el>=30)&(el<60)&land)*1.0)-273.15,
           nightFrozen=g((acc['ts']<273.15)*1.0,((el<0)&land)*1.0),landRain=g(acc['pr'],land*1.0)*86400*360*1000,note=note,run=run,w=0.5 if note else 1,name=run)
    out.append(e); fields[run]=dict(ts=acc['ts'],tas=acc['tas'],pr=acc['pr']*86400*360*1000,lat=lat,lon=lon)
    print('%-9s %s %5.3fd F %4.0f | G %.1f D %.1f N %.1f max %.0fC alb %.3f asr %.0f iceFree %.2f land %.1f/%.1f nightFrozen %.2f landRain %.0f %s'%(run,star[:4],per,flux,e['G'],e['D'],e['N'],e['maxC'],e['alb'],e['asr'],e['iceFree'],e['l1'],e['l2'],e['nightFrozen'],e['landRain'],note))
json.dump(out,open('pm/anchors.json','w'),indent=1)
for r in sys.argv[1:]: np.savez('pm/gcm_%s.npz'%r,**fields[r])
