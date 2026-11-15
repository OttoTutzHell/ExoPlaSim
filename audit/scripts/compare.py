"""Before/after table and figure for benchmark branches.
Usage: compare.py OUTPREFIX NYEARS label=rundir [label=rundir ...]   (first run is the reference)
Averages the last NYEARS annual files of each run. Optional sslon=DEG (default 180, as thai_hab1.py) and flux=W/m2 (default 900)."""
import sys, glob, numpy as np, netCDF4 as nc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SSLON=180.0
def load(rundir,nyears):
    files=sorted(f for f in glob.glob(rundir+"/MOST.[0-9]*.nc"))[-nyears:]
    acc={}
    for f in files:
        d=nc.Dataset(f); lat=np.array(d["lat"][:]); lon=np.array(d["lon"][:])
        for k in d.variables:
            if d[k].ndim==3: acc[k]=acc.get(k,0)+np.array(d[k][:]).mean(0)/len(files)
        d.close()
    return lat,lon,acc,len(files)
def stats(lat,lon,v,flux):
    LON,LAT=np.meshgrid(lon,lat); w=np.cos(np.radians(LAT))
    mu=np.cos(np.radians(LAT))*np.cos(np.radians(LON-SSLON)); ss=mu>0.9
    # hemisphere weights: grid columns lying exactly on the terminator count half to each side
    day=np.where(mu>1e-6,1.0,np.where(mu<-1e-6,0.0,0.5)); night=1.0-day
    def g(x,m=None):
        m=np.ones_like(w) if m is None else m*1.0
        return float((x*w*m).sum()/(w*m).sum())
    mmyr=86400*360*1000.
    r={}
    r["T_s global (K)"]=g(v["ts"]); r["T_s dayside (K)"]=g(v["ts"],day); r["T_s nightside (K)"]=g(v["ts"],night)
    r["T_s max (K)"]=float(v["ts"].max()); r["T_s min (K)"]=float(v["ts"].min())
    r["T_s substellar cap, <26 deg (K)"]=g(v["ts"],ss)
    r["ice-free ocean fraction"]=g((v["sic"]<0.5)*1.0)
    r["sea-ice thickness, global mean (m)"]=g(v["sit"])
    r["planetary albedo"]=1-g(v["rst"])/(flux/4)
    r["absorbed SW, TOA (W/m2)"]=g(v["rst"]); r["absorbed SW, surface (W/m2)"]=g(v["rss"])
    r["absorbed SW, atmosphere (W/m2)"]=g(v["rst"])-g(v["rss"])
    r["atmosphere share of absorbed SW"]=1-g(v["rss"])/g(v["rst"])
    r["dayside surface SW (W/m2)"]=g(v["rss"],day)
    r["OLR (W/m2)"]=-g(v["rlut"])
    r["TOA net (W/m2)"]=g(v["ntr"]); r["surface net (W/m2)"]=g(v["hfns"])
    r["precipitation, global (mm/yr)"]=g(v["pr"])*mmyr
    r["precipitation, substellar cap (mm/yr)"]=g(v["pr"],ss)*mmyr
    r["evaporation, global (mm/yr)"]=abs(g(v["evap"]))*mmyr
    r["column water vapour (kg/m2)"]=g(v["prw"]); r["cloud cover"]=g(v["clt"]); r["cloud cover, dayside"]=g(v["clt"],day)
    if "ohdf" in v: r["ocean transport heating, mean |.| (W/m2)"]=g(abs(v["ohdf"])); r["ocean transport heating, global mean (W/m2)"]=g(v["ohdf"])
    return r
if __name__=="__main__":
    out=sys.argv[1]; ny=int(sys.argv[2]); flux=900.0
    runs=[a.split("=",1) for a in sys.argv[3:]]
    for k,v in list(runs):          # optional settings: sslon=DEG flux=W/m2
        if k=="sslon": SSLON=float(v); runs.remove([k,v])
        if k=="flux": flux=float(v); runs.remove([k,v])
    data=[(lab,)+load(d,ny) for lab,d in runs]
    st=[stats(lat,lon,v,flux) for _,lat,lon,v,_ in data]
    keys=[]; [keys.append(k) for s in st for k in s if k not in keys]
    lines=["| quantity | "+" | ".join(l for l,*_ in data)+" |","|---|"+"---|"*len(data)]
    for k in keys:
        fmt="%.3f" if ("fraction" in k or "albedo" in k or "share" in k or "cloud" in k or "global mean (W" in k) else "%.1f"
        lines.append("| "+k+" | "+" | ".join((fmt%s[k]) if k in s else "—" for s in st)+" |")
    lines.append("| years averaged | "+" | ".join(str(n) for *_,n in data)+" |")
    open(out+".md","w").write("\n".join(lines)+"\n"); print("\n".join(lines))
    TMIN=10*np.floor(min(d[3]["ts"].min() for d in data)/10); TMAX=10*np.ceil(max(d[3]["ts"].max() for d in data)/10)
    n=len(data); fig,ax=plt.subplots(3,n,figsize=(4.2*n,8.2),squeeze=False)
    for i,(lab,lat,lon,v,_) in enumerate(data):
        a=ax[0,i]; c=a.pcolormesh(lon,lat,v["ts"],vmin=TMIN,vmax=TMAX,cmap="RdYlBu_r",shading="auto")
        a.contour(lon,lat,v["sic"],levels=[0.5],colors="k",linewidths=1.2); a.set_title(lab+"\nsurface temperature (K), ice edge",fontsize=9)
        fig.colorbar(c,ax=a,orientation="horizontal",pad=0.14,fraction=0.05)
        a=ax[1,i]; c=a.pcolormesh(lon,lat,v["pr"]*86400*360*1000,vmin=0,vmax=4000,cmap="Blues",shading="auto"); a.set_title("precipitation (mm/yr)",fontsize=9)
        fig.colorbar(c,ax=a,orientation="horizontal",pad=0.14,fraction=0.05)
        a=ax[2,i]; ref=data[0][3]
        if i==0:
            c=a.pcolormesh(lon,lat,v["rss"],vmin=0,vmax=500,cmap="magma",shading="auto"); a.set_title("surface absorbed SW (W/m2)",fontsize=9)
        else:
            c=a.pcolormesh(lon,lat,v["ts"]-ref["ts"],vmin=-15,vmax=15,cmap="RdBu_r",shading="auto"); a.set_title("surface temperature minus "+data[0][0]+" (K)",fontsize=9)
        fig.colorbar(c,ax=a,orientation="horizontal",pad=0.14,fraction=0.05)
    for a in ax.ravel(): a.set_xticks([0,90,180,270]); a.set_yticks([-60,0,60])
    fig.tight_layout(); fig.savefig(out+".png",dpi=110); print("wrote",out+".png")
