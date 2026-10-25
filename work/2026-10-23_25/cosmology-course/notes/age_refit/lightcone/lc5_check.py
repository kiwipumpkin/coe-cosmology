"""Visible sky = concentric circle centred on us NOW (radius rho=c(t_now - t*)), made of elements that co-move with our group.
Per direction: exact-GR past null ray from us now; emission event = point where the element on the ray, carried to t_now along its own
LTB worldline, sits on that circle (|x_elem(t_now) - x_us(t_now)| = rho, areal flat distance). Path drawn as actually flown (elements and
us displacing radially during flight). Path length = c(t_now - t_e) (coordinate light travel), also areal path length. z exact GR."""
import json,pickle,numpy as np; from lc2 import *
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
r=json.load(open('res_exp.json')); L,tn=make2('exp',r['x'],1e-3); ts=tstar(L,tn); X0=np.array([0.]); RO=lambda t: L.R(t,X0)[0]
o=obs(L,tn); lam=o['Hbar']*F.C/(70/F.KMSMPC); Ron=RO(tn); rho=tn-ts
def g(t,y):
    Rn=L.R(tn,np.array([y[0]]))[0]; return np.sqrt(Rn**2+Ron**2-2*Rn*Ron*np.cos(y[1]))-rho
g.terminal=True
res=[]
for th in np.arange(0,360,30.):
    rr=ray(L,tn,np.radians(th),tmin=ts,rtol=1e-9,budget=40000,evx=g); q=dict(th=th,stop=rr['stop'],rr=rr)
    if rr['stop']=='circle_hit':
        Rr=L.R(rr['t'],rr['x']); X=Rr*np.cos(rr['phi']); Y=Rr*np.sin(rr['phi']); al=float(np.sum(np.hypot(np.diff(X),np.diff(Y))))
        q.update(te=rr['t'][-1],z=float(np.expm1(rr['lnz1'][-1])),light_travel=tn-rr['t'][-1],areal_path=al,beyond=bool(rr['nc'].max()>0),xe=rr['x'][-1],phie=rr['phi'][-1])
    res.append(q); print(th,q['stop'],q.get('te',0)/tn,q.get('z'),flush=True)
pickle.dump(dict(ts=ts,tn=tn,rho=rho,res=[{k:v for k,v in q.items() if k!='rr'} for q in res]),open('lc5_check.pkl','wb'))
fig,ax=plt.subplots(figsize=(19,13)); ax.set_aspect('equal'); c=np.linspace(0,2*np.pi,300)
ax.plot((Ron+rho*np.cos(c))*lam,rho*np.sin(c)*lam,c='orange',lw=4,label=f'visible circle around us now, radius c(t_now−t*) = {rho*lam:.2f} Gpc')
tt=np.linspace(ts,tn,100); ax.plot([RO(t)*lam for t in tt],0*tt,'k-',lw=3,label='our path'); ax.plot(Ron*lam,0,'k*',ms=28,label='us now')
for q in res:
    p=q['rr']; Rr=L.R(p['t'],p['x']); X=Rr*np.cos(p['phi'])*lam; Y=Rr*np.sin(p['phi'])*lam
    if 'z' not in q: ax.plot(X,Y,'r--',lw=1.5); ax.plot(X[-1],Y[-1],'rx',ms=14,mew=3); continue
    ax.plot(X,Y,'-',c='C0',lw=2.5); ax.plot(X[-1],Y[-1],'o',c='C0',ms=12)          # emission event (where the element was at t_e)
    Rn=L.R(tn,np.array([q['xe']]))[0]*lam; ax.plot([X[-1],Rn*np.cos(q['phie'])],[Y[-1],Rn*np.sin(q['phie'])],':',c='green',lw=2)  # element's radial displacement since
    ax.plot(Rn*np.cos(q['phie']),Rn*np.sin(q['phie']),'s',c='green',ms=10)
    ax.annotate(f"{q['th']:.0f}°: z={q['z']:.2f}\npath {q['light_travel']*lam:.2f} Gpc",(X[-1],Y[-1]),fontsize=13,xytext=(6,6),textcoords='offset points')
ax.plot(0,0,'rx',ms=25,mew=4); ax.text(0.15,-0.5,'CoE',fontsize=22,color='r')
nb=sum(q.get('beyond',False) for q in res); ne=sum('z' in q for q in res)
ax.set_title(f'exact GR, exp best fit (U={o["U"]:.2f}c): blue = photon path as flown, ● emission event, ■ that element now (on the circle)\nemitting directions {ne}/12, beyond CoE {nb}; red x = ray reached t* without meeting the circle',fontsize=19)
ax.set_xlabel('Gpc (areal R cos φ)',fontsize=20); ax.tick_params(labelsize=18); ax.legend(fontsize=14,loc='lower left'); ax.grid()
fig.tight_layout(); fig.savefig('lc5_check.png',dpi=80); print('FIN',flush=True)
