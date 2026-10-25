import sys,pickle; sys.path.insert(0,'.'); from grcore import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=pickle.load(open('fast.pkl','rb')); x=np.degrees(th)
for v0,to in ((0.9,11.0),(0.99,106.0)):
    f,ax=plt.subplots(figsize=(14,8))
    for u,st,l in ((0,'-','not expanding'),(0.9999,'--','expanding ≈c')):
        r=norm(R[(to,0.1,u,v0)][1][0]); d,q,_,mn=stats(r)
        ax.plot(x,r,st,color='C0' if u==0 else 'C3',lw=3,label=f"{l}: dipole {d:.4f}, quad {q:.4f}, coldest {mn:.4f}")
    ax.axhline(0,color='gray'); ax.set_xlim(0,360); ax.grid(alpha=.3); ax.set_ylabel("ΔT/T",fontsize=16); ax.tick_params(labelsize=14)
    ax.set_xlabel(f"direction (deg); exact GR; start {v0}c toward 0°, GM/(Rc²) = 0.1 at 180°, viewed at t = {to:g}",fontsize=15)
    ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.02),fontsize=14); plt.tight_layout(); plt.savefig(f'speeds_gr_v{v0}.png',dpi=120); plt.close()
