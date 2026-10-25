import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
X=lambda t,a: (np.sqrt(1+(a*t)**2)-1)/a
V=lambda t,a: a*t/np.sqrt(1+(a*t)**2)
th=np.radians(np.arange(0,361,1.0)); to=2.0
plt.rcParams.update({"font.size":18}); fig,ax=plt.subplots(figsize=(13,8))
for a in (0.05,0.1,0.2,0.5):
    T=[];TE=[]
    for t in th:
        g=lambda te:(X(to,a)-np.cos(t)-X(te,a))**2+np.sin(t)**2-(to-te)**2
        te=brentq(g,-0.0,to-1e-9) if g(0)<0 else 0.0
        dx=X(to,a)-np.cos(t)-X(te,a); dy=-np.sin(t); L=np.hypot(dx,dy); nx=dx/L
        ve,vo=V(te,a),V(to,a); ge,go=1/np.sqrt(1-ve*ve),1/np.sqrt(1-vo*vo)
        T.append(go*(1-vo*nx)/(ge*(1-ve*nx))); TE.append(te)
    T=np.array(T); d=T/T.mean()-1; print("te front/side/back",TE[0],TE[90],TE[180],"d",d[0],d[90],d[180])
    m=np.fft.rfft(d[:-1])/360; print(a,"dipole %.4f quad %.5f"%(2*abs(m[1]),2*abs(m[2])),"aL/c2",a)
    ax.plot(np.degrees(th),d,lw=3,label=f"a = {a} c²/L")
ax.axhline(0,color="gray"); ax.set_xlabel("direction of emitter from acceleration axis (deg)"); ax.set_ylabel("ΔT/T"); ax.set_xlim(0,360); ax.grid(alpha=.3)
ax.legend(loc="upper center",bbox_to_anchor=(0.5,1.13),ncol=4,fontsize=16)
plt.tight_layout(); plt.savefig("accel_levels.png",dpi=120)
