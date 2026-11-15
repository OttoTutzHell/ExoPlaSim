import sys
d=sys.argv[1]
t=open(d+'/template.html').read()
t=t.replace('/*DATA*/',open(d+'/kidata.js').read()).replace('/*MODEL*/',open(d+'/model.js').read())
open(d+'/ki-engir-planet-model.html','w').write(t)
print(len(t))
