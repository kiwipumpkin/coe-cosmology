"""ONE emitting circle C(t_e): centred on O(t_e), radius c(t_e - t*) (areal flat distance), all points emit at the same t_e.
For each arrival direction theta, solve the arrival time t_a at which our exact-GR past null ray, traced back to t_e, lands on that circle.
If all t_a were equal (= now), the single circle is what we see now. t_e chosen so the 0-deg direction arrives exactly now (t_a=tn)."""
import json,pickle,numpy as np; from lc2 import *
from scipy.optimize import brentq
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
r=json.load(open('res_exp.json')); L,tn=make2('exp',r['x'],1e-3); ts=tstar(L,tn); X0=np.array([0.]); RO=lambda t: L.R(t,X0)[0]
o=obs(L,tn); lam=o['Hbar']*F.C/(70/F.KMSMPC)
def at(ta,th,te):
    rr=ray(L,ta,np.radians(th),tmin=te,rtol=1e-9,budget=40000)
    if rr['stop']!='reached_tmin': return None,rr
    Rr=L.R(te,np.array([rr['x'][-1]]))[0]; Ro=RO(te); d=np.sqrt(Rr**2+Ro**2-2*Rr*Ro*np.cos(rr['phi'][-1])); return d-(te-ts),rr
# t_e from the 0-deg direction arriving now
te=brentq(lambda t: at(tn,0.,t)[0],0.3*tn,0.95*tn,xtol=1e-6*tn); print('t_e/tn',te/tn,flush=True)
res=[]
for th in np.arange(0,360,30.):
    f=lambda ta: (at(ta,th,te)[0] if at(ta,th,te)[0] is not None else np.nan)
    grid=np.linspace(te*1.0005,3*tn,25); v=np.array([f(t) for t in grid]); i=np.where(np.isfinite(v[:-1])&np.isfinite(v[1:])&(v[:-1]*v[1:]<=0))[0]
    if len(i): ta=brentq(f,grid[i[0]],grid[i[0]+1],xtol=1e-6*tn); _,rr=at(ta,th,te); z=float(np.expm1(rr["lnz1"][-1])); res.append(dict(th=th,ta=ta,z=z,T=3000/(1+z),rr=rr,beyond=bool(rr["nc"].max()>0)))
    else: res.append(dict(th=th,ta=None))
    print(th,res[-1].get('ta',None) and res[-1]['ta']/tn,res[-1].get('z'),flush=True)
pickle.dump(dict(te=te,ts=ts,tn=tn,res=[{k:v for k,v in q.items() if k!='rr'} for q in res]),open('lc4_check.pkl','wb'))
fig,ax=plt.subplots(figsize=(18,12)); ax.set_aspect('equal'); thc=np.linspace(0,2*np.pi,300)
ax.plot((RO(te)+(te-ts)*np.cos(thc))*lam,(te-ts)*np.sin(thc)*lam,c='orange',lw=4,label=f'emitting circle C(t_e), t_e={te/tn:.3f} t_now, radius c(t_e−t*)')
ax.plot(RO(te)*lam,0,'o',c='orange',ms=16,label='us at t_e (circle centre)')
tt=np.linspace(te,max([q['ta'] for q in res if q.get('ta')]+[tn]),100); ax.plot([RO(t)*lam for t in tt],0*tt,'k-',lw=3,label='our path')
ax.plot(RO(tn)*lam,0,'k*',ms=26,label='us now')
cm=plt.cm.coolwarm
for q in res:
    if not q.get('ta'): continue
    p=q['rr']; Rr=L.R(p['t'],p['x']); X=Rr*np.cos(p['phi'])*lam; Y=Rr*np.sin(p['phi'])*lam
    c=cm((q['ta']/tn-0.5)/1.5); ax.plot(X,Y,'-',c=c,lw=2.5); ax.plot(X[-1],Y[-1],'o',c=c,ms=12); ax.plot(RO(q['ta'])*lam,0,'|',c=c,ms=30,mew=4)
    ax.annotate(f"{q['th']:.0f}°: arrives {q['ta']/tn:.2f} t_now",(X[-1],Y[-1]),fontsize=13,xytext=(5,5),textcoords='offset points')
ax.plot(0,0,'rx',ms=25,mew=4); ax.text(0.15,-0.5,'CoE',fontsize=22,color='r')
ax.set_title('One circle, one emission time t_e: exact-GR null paths from each circle point to our worldline\n(ticks on our path = where each path actually reaches us; colour = arrival time)',fontsize=20)
ax.set_xlabel('Gpc (areal R cos φ)',fontsize=20); ax.tick_params(labelsize=18); ax.legend(fontsize=14,loc='lower left'); ax.grid()
fig.tight_layout(); fig.savefig('lc4_check.png',dpi=80); print('FIN',flush=True)
