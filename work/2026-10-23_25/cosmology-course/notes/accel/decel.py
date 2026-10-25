import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp; from scipy.optimize import brentq
v0=0.3; to=3.0; Mp=np.array([-2.0,0.0]); th=np.radians(np.arange(0,361,1.0))
def traj(p0,mu):
    def f(t,y):
        d=y[:2]-Mp; r=np.linalg.norm(d); return [y[2],y[3],*(-mu*d/r**3)]
    return solve_ivp(f,(0,to),[p0[0],p0[1],v0,0.0],dense_output=True,rtol=1e-11,atol=1e-12).sol
plt.rcParams.update({"font.size":18}); fig,ax=plt.subplots(figsize=(13,8))
for mu in (0.0,0.02,0.05,0.1):
    O=traj(np.zeros(2),mu); yo=O(to); xo,vo=yo[:2],yo[2:]
    d=[]
    for t in th:
        E=traj(np.array([np.cos(t),np.sin(t)]),mu)
        te=brentq(lambda s: np.linalg.norm(xo-E(s)[:2])-(to-s),0,to-1e-6)
        ye=E(te); n=(xo-ye[:2])/np.linalg.norm(xo-ye[:2]); ve=ye[2:]
        ge,go=1/np.sqrt(1-ve@ve),1/np.sqrt(1-vo@vo)
        dop=go*(1-vo@n)/(ge*(1-ve@n))
        re=np.linalg.norm(ye[:2]-Mp); ro=np.linalg.norm(xo-Mp)
        d.append(dop*np.sqrt((1-2*mu/re)/(1-2*mu/ro)))
    d=np.array(d); r=d/d.mean()-1; m=np.fft.rfft(r[:-1])/360
    print(mu,"v_obs %.4f"%vo[0],"dipole %.4f quad %.4f max %.4f min %.4f"%(2*abs(m[1]),2*abs(m[2]),r.max(),r.min()))
    ax.plot(np.degrees(th),r,lw=3,label=f"GM/(Rc²) = {mu}")
ax.axhline(0,color="gray"); ax.set_xlabel("direction of emitter (deg); start 0.3c toward 0°, mass at 180°"); ax.set_ylabel("ΔT/T")
ax.set_xlim(0,360); ax.grid(alpha=.3); ax.legend(loc="upper center",bbox_to_anchor=(0.5,1.13),ncol=4,fontsize=15)
plt.tight_layout(); plt.savefig("decel.png",dpi=120)
