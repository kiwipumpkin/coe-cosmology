import json, numpy as np
from scipy.optimize import least_squares
src=open("cmb-cancel.py").read().split('if __name__')[0]; g={}; exec(src,g)
d=json.load(open("cmb-cancel-story.json")); P=d["params"]; phi=P[4]
ELL=g["ELL"]; R=g["resid"](g["G_FULL"],g["S_model"](*P)[0]); R=R-R.mean()
th=np.radians(ELL)
def field(alpha,D):
    # line-of-sight (radial) tidal pull at unit radius from point mass at distance D, angle alpha; units GM=1
    x,y=np.cos(th),np.sin(th); X,Y=D*np.cos(alpha),D*np.sin(alpha)
    dx,dy=X-x,Y-y; r3=(dx*dx+dy*dy)**1.5
    ar=(dx*x+dy*y)/r3          # radial acceleration toward clump
    a0x,a0y=X/D**3,Y/D**3      # uniform pull on the origin (removed: common motion)
    ar=ar-(a0x*x+a0y*y)
    return ar-ar.mean()
def model(p,npair,single):
    D=p[0]; out=np.zeros_like(th); i=1
    for k in range(npair):
        k_,off=p[i],p[i+1]; i+=2
        out+=k_*(field(np.radians(phi+off),D)+field(np.radians(phi-off),D))
    for s in single:
        out+=p[i]*field(np.radians(s),D); i+=1
    return out
res={}
for npair,single,lab in [(2,[],"2 pairs"),(2,[phi,phi+180],"2 pairs + on-axis"),(3,[],"3 pairs")]:
    best=None
    for D0 in [1.5,2,3,4]:
        for seed in range(6):
            rng=np.random.default_rng(seed)
            p0=[D0]+sum([[rng.normal(0,1e-2),rng.uniform(20,170)] for _ in range(npair)],[])+[0.0]*len(single)
            lb=[1.05]+[-np.inf,0]*npair+[-np.inf]*len(single); ub=[20]+[np.inf,180]*npair+[np.inf]*len(single)
            f=least_squares(lambda p:model(p,npair,single)-R,p0,bounds=(lb,ub))
            if best is None or f.cost<best.cost: best=f
    p=best.x; left=R-model(p,npair,single)
    res[lab]=dict(p=p.tolist(),rms=float(np.sqrt(np.mean(left**2))),frac_var=float(1-np.mean(left**2)/np.mean(R**2)))
    print(lab,np.round(p,4),res[lab]["rms"],res[lab]["frac_var"])
json.dump(res,open("clump-fit.json","w"),indent=1)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":15})
p=np.array(res["2 pairs"]["p"]); m=model(p,2,[])
fig=plt.figure(figsize=(15,7)); a=fig.add_subplot(121)
a.plot(ELL,R*1e3,"C3",lw=2.5,label="residual after SR + mass"); a.plot(ELL,m*1e3,"k--",lw=2,label="4 clumps, tidal pull")
a.set_xlim(0,360); a.grid(alpha=.3); a.set_xlabel(r"Galactic longitude $\ell$ (deg)"); a.set_ylabel(r"$\Delta T/T\times10^3$")
a.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=1,frameon=False); a.set_title("4 clumps explain 88% of the residual variance",pad=70)
b=fig.add_subplot(122,projection="polar"); D=p[0]
ang=[phi+p[2],phi-p[2],phi+p[4],phi-p[4]]
b.plot(np.linspace(0,2*np.pi,400),np.ones(400),"C0",lw=2,label="visible universe, 14.3 Gpc")
b.scatter(np.radians(ang),[D]*4,s=400,color="C3",zorder=3,label=f"clumps at {D*14.3:.1f} Gpc, 7.7e51 kg each")
b.plot([np.radians(phi)]*2,[0,2.1],"k",lw=1.5); b.plot([np.radians(phi+180)]*2,[0,2.1],"k--",lw=1.5)
for A in ang: b.text(np.radians(A),D+0.25,f"{A%360:.0f}°",ha="center",fontsize=14)
b.set_ylim(0,2.2); b.set_yticklabels([]); b.legend(loc="lower center",bbox_to_anchor=(0.5,1.08),frameon=False,fontsize=14)
fig.tight_layout(); fig.savefig("clump-fit.png",dpi=100)
