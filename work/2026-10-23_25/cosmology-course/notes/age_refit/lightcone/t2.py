import numpy as np; from lc import *
L,tn=make('exp',[-1,-1,0,np.log10(1.5),1.0,-0.3]); o=obs(L,tn)
r=ray(L,tn,0.0,tmin=1e-4*tn,rtol=1e-8); lr=np.log(L.rho_rel(r['t'],r['x'])/o['rho'])
for k in np.linspace(0,len(r['t'])-1,15).astype(int): print(r['t'][k]/tn, r['x'][k], np.exp(r['x'][k]), lr[k], r['lnz1'][k], L.R(r['t'][k],np.array([r['x'][k]]))[0])
print(L.rho_rel(tn*1e-3,np.array([0.,2,4,6])), L.m(np.array([0.,2,4,6])), L.m(np.array([0.,2,4,6]),1))
