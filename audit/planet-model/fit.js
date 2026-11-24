const fs=require('fs'); eval(fs.readFileSync('kidata.js','utf8').replace('const KI_DATA','var KI_DATA'));
const {createKiModel}=require('./model.js');
const m=createKiModel(KI_DATA);
const A=JSON.parse(fs.readFileSync('anchors.json','utf8'));
function run(a,over,n){ Object.assign(m.P,{useOrbit:1,spectral:1,starT:a.starT,starL:a.L,starM:a.M,period:a.period,spec:a.starT===2566?{fvis:0.015,ratio:[1.046,3.192,4.331]}:null},over); m.reset({T0:a.G}); for(let i=0;i<n;i++) m.step(i<2?0.5:2,true); return m.stats(); }
function cost(over,n=12,show=false){ let c=0; for(const a of A){ const s=run(a,over,n);
  const e=[(s.global-a.G)/3,(s.day-a.D)/3,(s.night-a.N)/4,(s.albedo-a.alb)/0.03,(s.iceFree-a.iceFree)/0.15,((s.land1530-273.15)-a.l1)/5,((s.land3060-273.15)-a.l2)/5];
  c+=a.w*e.reduce((q,x)=>q+x*x,0);
  if(show) console.log(`${a.name} F ${s.S0.toFixed(0)} | G ${s.global.toFixed(1)} (${a.G}) D ${s.day.toFixed(1)} (${a.D}) N ${s.night.toFixed(1)} (${a.N}) alb ${s.albedo.toFixed(3)} (${a.alb}) land ${(s.land1530-273.15).toFixed(1)}/${(s.land3060-273.15).toFixed(1)} (${a.l1}/${a.l2}) max ${(s.max-273.15).toFixed(0)}C iceFree ${s.iceFree.toFixed(2)} colW ${s.colWater.toFixed(0)}`);} return c; }
const keys=process.argv[2].split(','), B={albCloud:[0.3,0.95],cloudMax:[0.4,0.98],cloudBase:[0.1,0.6],cloudSens:[0,0.15],tauDry:[0.1,0.8],tauWv:[0.1,1.2],Datm:[1,12],wvScale:[0.3,1.5],cloudLW:[0,0.6],Docean:[0,3]};
let best=JSON.parse(process.argv[3]||'{}'); for(const k of keys) if(!(k in best)) best[k]=m.P[k];
let bc=cost(best); console.log('start',bc.toFixed(1));
let step={}; for(const k of keys) step[k]=(B[k][1]-B[k][0])*0.15;
const iters=+process.argv[4]||6;
for(let it=0;it<iters;it++){ for(const k of keys){ for(const sg of [1,-1]){ const t={...best}; t[k]=Math.max(B[k][0],Math.min(B[k][1],best[k]+sg*step[k])); if(t[k]===best[k]) continue; const c=cost(t); if(c<bc){bc=c;best=t;break;} } } for(const k of keys) step[k]*=0.6; console.log(it,bc.toFixed(1),JSON.stringify(best)); }
cost(best,30,true); fs.writeFileSync('fit.json',JSON.stringify(best));
