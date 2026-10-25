import sys,json,numpy as np;sys.path.insert(0,'..');from modes import ring,build,addv
from scipy.optimize import brentq
VEXP=0.5
def ring_lorentz(v,N=360):
    g=1/np.sqrt(1-v*v);th=np.radians(np.arange(N)*360/N);b=np.array([v,0,0]);out=[]
    for t in th:
        kp=-np.array([np.cos(t),np.sin(t),0]);kx=(kp[0]+v)/(1+v*kp[0]);ky=kp[1]/(g*(1+v*kp[0]));k=np.array([kx,ky,0]);k/=np.linalg.norm(k);n=-k
        def rest(s):
            x=n*s;tl=-s;tp=g*(tl-v*x[0]);xp=np.array([g*(x[0]-v*tl),x[1],0]);return tp,xp
        f=lambda s:np.linalg.norm(rest(s)[1])-(1+VEXP*rest(s)[0])
        s=brentq(f,1e-9,100,xtol=1e-15,rtol=1e-15);tp,xp=rest(s);m=xp/np.linalg.norm(xp)
        ve=addv(b,VEXP*m);ge=1/np.sqrt(1-ve@ve);out.append(g*(1-b@k)/(ge*(1-ve@k)))
    return np.array(out)
st=lambda x:(lambda d:(float(np.sqrt(np.mean(d**2))),float(d.max()-d.min())))(x/x.mean()-1)
V=[0.01,0.1,0.3,0.5,0.9,0.99];rows=[]
for v in V:
    a=st(ring(**build(3,v,0,True)));b=st(ring_lorentz(v));rows.append(dict(v=v,old_rms=a[0],old_range=a[1],lorentz_rms=b[0],lorentz_range=b[1]));print(rows[-1])
lv=np.log(V);la=np.log([r['old_rms'] for r in rows])
p=np.polyfit(lv[:3],la[:3],1);pall=np.polyfit(lv,la,1)
out=dict(rows=rows,powerlaw_rms_low_speed=dict(exponent=p[0],coeff=float(np.exp(p[1]))),powerlaw_rms_all=dict(exponent=pall[0],coeff=float(np.exp(pall[1]))))
json.dump(out,open('mode3_speed.json','w'),indent=1);print(out['powerlaw_rms_low_speed'],out['powerlaw_rms_all'])
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':15})
fig,ax=plt.subplots(figsize=(11,7))
ax.loglog(V,[r['old_rms'] for r in rows],'o-',lw=2.5,label=f'(a) shape fixed in lab: RMS ≈ {np.exp(p[1]):.3g}·v^{p[0]:.2f} (fit v≤0.3)')
ax.loglog(V,[max(r['lorentz_rms'],1e-18) for r in rows],'s-',lw=2.5,label='(b) Lorentz-contracted (round in rest frame)')
ax.set_xlabel('system speed (c)');ax.set_ylabel('ring ΔT/T RMS');ax.legend(loc='lower center',bbox_to_anchor=(.5,1),frameon=False);fig.tight_layout();fig.savefig('mode3_speed.png',dpi=110)
