"""THAI Hab 1-like tidally locked aquaplanet (Fauchez et al. 2020 protocol values),
used as the published benchmark configuration. Usage: thai_hab1.py NAME YEARS [key=value ...]"""
import sys, os, exoplasim as exo
name=sys.argv[1]; years=int(sys.argv[2])
other={}; kw={}
res="T21"; ncpu=4; restart=None
for a in sys.argv[3:]:
    k,v=a.split("=",1)
    if k=="res": res=v
    elif k=="ncpus": ncpu=int(v)
    elif k=="restart": restart=os.path.abspath(v)
    elif "@" in k: other[k]=v
    else: kw[k]=float(v)
os.makedirs("runs",exist_ok=True); os.chdir("runs")
m=exo.Model(workdir=name,modelname=name,resolution=res,layers=10,ncpus=ncpu,precision=4,
            outputtype=".nc",hyperthreading=False,mpi_opts="--bind-to none --oversubscribe")
cfg=dict(flux=900.0,startemp=2600.0,rotationperiod=6.1,synchronous=True,substellarlon=180.0,
         eccentricity=0.0,obliquity=0.0,fixedorbit=True,radius=0.91,gravity=9.12,
         pN2=1.0-400e-6,pCO2=400e-6,aquaplanet=True,ozone=False,timestep=30.0,
         physicsfilter="gp|exp|sp",tlcontrast=100.0,snapshots=None,mldepth=50.0,
         restartfile=restart,otherargs=other)
cfg.update(kw)
m.configure(**cfg)
m.exportcfg()
m.run(years=years,crashifbroken=True,clean=True)
