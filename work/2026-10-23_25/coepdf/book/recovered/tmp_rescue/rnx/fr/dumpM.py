exec(open('recon_n3.py').read().split("cols=['purple'")[0])
import json;out={'l':list(map(float,l))}
for nm,Y,mb,k in (('A',A0,AFl,0),('B',B0,BLl,1)):
  m=mb.copy()
  for n,p in steps2[:3]:m=m+T(n,p,th,k)
  r=100*(Y-m)/Y;out[nm]=[float(np.mean(r[np.abs((l-c+180)%360-180)<=25])) for c in l]
json.dump(out,open('/tmp/M25.json','w'))
