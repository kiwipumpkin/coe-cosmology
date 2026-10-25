import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
v=0.2; th=np.radians(np.arange(0,361,1.0)); M=np.array([-2.0,0.0])
plt.rcParams.update({"font.size":18}); fig,ax=plt.subplots(figsize=(13,8))
for mu in (0.05,0.1,0.2,0.3):
    d=[]
    for t in th:
        g=lambda te: np.hypot(-(np.cos(t)+v*te),-np.sin(t))-(-te)
        te=brentq(g,-10,0)
        xe=np.array([np.cos(t)+v*te,np.sin(t)]); xo=np.zeros(2); n=(xo-xe)/np.linalg.norm(xo-xe)
        dop=(1-v*n[0])/(1-v*n[0])   # same velocity: cancels exactly
        re=np.linalg.norm(xe-M); ro=np.linalg.norm(xo-M)
        d.append(dop*np.sqrt((1-2*mu/re)/(1-2*mu/ro)))
    d=np.array(d); r=d/d.mean()-1
    m=np.fft.rfft(r[:-1])/360; print(mu,"dipole %.4f quad %.4f max %.4f min %.4f"%(2*abs(m[1]),2*abs(m[2]),r.max(),r.min()))
    ax.plot(np.degrees(th),r,lw=3,label=f"GM/(Rc²) = {mu}")
ax.axhline(0,color="gray"); ax.set_xlabel("direction of emitter (deg); moving toward 0° at 0.2c, mass at 180°"); ax.set_ylabel("ΔT/T")
ax.set_xlim(0,360); ax.grid(alpha=.3); ax.legend(loc="upper center",bbox_to_anchor=(0.5,1.13),ncol=4,fontsize=15)
plt.tight_layout(); plt.savefig("mass180.png",dpi=120)
