import sys,pickle; sys.path.insert(0,'.'); from grcore import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=pickle.load(open('speeds.pkl','rb'))
fig,ax=plt.subplots(figsize=(14,8.5)); x=np.degrees(th)
for v0,col in ((0.1,'C0'),(0.3,'C3'),(0.5,'C2')):
    for u,st,l in ((0,'-','not expanding'),(0.9999,'--','expanding ≈c')):
        g=R[(0.1,u,v0)][1][0]
        if np.isnan(g).any():
            r=g/np.nanmean(g)-1; lab=f"start {v0}c, {l}: 148°–212° not yet visible at t = 2 (no dipole/quad)"
        else:
            r=norm(g); d,q,_,_=stats(r); lab=f"start {v0}c, {l}: dipole {d:.4f}, quad {q:.4f}"
        ax.plot(x,r,st,color=col,lw=3,label=lab)
ax.axhline(0,color='gray'); ax.set_xlim(0,360); ax.grid(alpha=.3); ax.set_ylabel("ΔT/T",fontsize=16); ax.tick_params(labelsize=14)
ax.set_xlabel("direction (deg); exact GR; motion toward 0°, GM/(Rc²) = 0.1 at 180°, viewed at t = 2",fontsize=15)
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.02),ncol=2,fontsize=12); plt.tight_layout(); plt.savefig('decel_overlay_speeds_gr.png',dpi=120)
