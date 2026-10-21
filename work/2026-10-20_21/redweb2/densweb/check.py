import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
CC=[-0.041,-0.035,0.1353,0.0];Rb=1.5*(1+1/2.937);aC=np.array([1,0,np.tan(np.radians(15))]);aC/=np.linalg.norm(aC);PC=42.0*aC
shp=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)))
rng=np.random.default_rng(0);n=4000000;u=rng.random(n);r=(0.5**3+u*(20**3-0.5**3))**(1/3);v=rng.normal(size=(n,3));v/=np.linalg.norm(v,axis=1)[:,None]
dU=np.linalg.norm(r[:,None]*v-PC,axis=1);E=np.linspace(23,61,20);c=0.5*(E[1:]+E[:-1]);Vs=np.histogram(dU,E)[0]
fig,ax=plt.subplots(1,2,figsize=(12,4.5));o=np.load('orig.npz');sets=[('seed-2000 web (orig)',o['arr_0'],o['arr_1'])]+[(f'densweb s{s}',*[np.load(f'densweb_s{s}.npz')[k] for k in 'XM']) for s in (1,2,3)]
ref=shp(c/42);ref/=np.average(ref,weights=Vs);out={}
for nm,X,M in sets:
  d=np.linalg.norm(X-PC,axis=1);pos=M>0
  for a,(w,lab) in zip(ax,[(M*pos,'positive mass (clusters+walls)'),(np.ones(len(M)),'body number')]):
    h=np.histogram(d,E,weights=w)[0]/Vs;h/=np.average(h,weights=Vs);a.plot(c,h,'o-',ms=3,label=nm);a.set_title('shell-avg '+lab+' density / mean')
  hn=np.histogram(d,E,weights=M)[0]/Vs;out[nm]=(hn/np.average(np.histogram(d,E,weights=M*pos)[0]/Vs,weights=Vs)).tolist()
for a in ax:a.plot(c,ref,'k-',lw=2,label='red ρ(r) shape (normalised)');a.set_xlabel('distance from CoE r (Gpc)');a.legend(fontsize=8);a.set_ylim(.9,1.1)
plt.tight_layout();plt.savefig('/workspace/redweb2/density_check.png',dpi=120)
import json;json.dump(dict(r=c.tolist(),rho=ref.tolist(),net_over_pos=out),open('density_check.json','w'))
print(np.round(ref,4))
