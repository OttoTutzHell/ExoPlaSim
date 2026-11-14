"""Ki-Engir T21 branch with the vault's f1800_rc070 settings, for demonstrating the new options.
Usage: kiengir_t21.py NAME YEARS BOUNDARYDIR RESTART [ncpus=N] [key=value ...]
key=value pairs go to Model.configure(); keys containing '@' go to otherargs."""
import sys, os, exoplasim as exo
name=sys.argv[1]; years=int(sys.argv[2]); B=os.path.abspath(sys.argv[3]); rst=os.path.abspath(sys.argv[4])
other={"rcritmod@rainmod_namelist":"0.70","tswr1@radmod_namelist":"0.3","tswr2@radmod_namelist":"0.18","tswr3@radmod_namelist":"0.002"}
kw={}; ncpu=2; recompile=False
for a in sys.argv[5:]:
    k,v=a.split("=",1)
    if k=="ncpus": ncpu=int(v)
    elif k=="recompile": recompile=bool(int(v))   # needed once: the executable cache ignores precision
    elif "@" in k: other[k]=v
    else: kw[k]=float(v)
os.makedirs("runs",exist_ok=True); os.chdir("runs")
m=exo.Model(workdir=name,modelname=name,resolution="T21",layers=10,ncpus=ncpu,precision=8,outputtype=".nc",
            hyperthreading=False,recompile=recompile,mpi_opts="--bind-to none --oversubscribe")
cfg=dict(flux=1800.0,startemp=2879.0,rotationperiod=6.097,synchronous=True,substellarlon=0.0,
         fixedorbit=True,eccentricity=0.0,obliquity=0.0,radius=5861.0/6371.0,gravity=7.99,
         pN2=0.999,pCO2=0.001,mldepth=50.0,seaice=True,
         landmap=B+"/Ki_T21_surf_0172.sra",topomap=B+"/Ki_T21_surf_0129.sra",
         timestep=30.0,physicsfilter="gp|exp|sp",tlcontrast=100.0,snapshots=None,restartfile=rst,otherargs=other)
cfg.update(kw); m.configure(**cfg); m.exportcfg()
for y in range(years): m.run(years=1,crashifbroken=True)
