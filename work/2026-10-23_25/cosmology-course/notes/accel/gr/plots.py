import sys,pickle; sys.path.insert(0,'.'); from grcore import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=pickle.load(open('results.pkl','rb')); x=np.degrees(th)
def fig(keys,labs,xl,fn):
    f,ax=plt.subplots(figsize=(14,8.5))
    for i,(k,l) in enumerate(zip(keys,labs)):
        w,g=R[k]; gr=norm(g[0])
        for r,st,n in ((w,'-','weak'),(gr,'--','GR')):
            if r is None: continue
            d,q,_,_=stats(r); ax.plot(x,r,st,color=f'C{i}',lw=3,label=f"{l}, {n}: dipole {d:.4f}, quad {q:.4f}")
    ax.axhline(0,color='gray'); ax.set_xlim(0,360); ax.grid(alpha=.3); ax.set_ylabel("ΔT/T",fontsize=16); ax.tick_params(labelsize=14)
    ax.set_xlabel(xl,fontsize=15); ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.02),ncol=2,fontsize=12); plt.tight_layout(); plt.savefig(fn,dpi=120); plt.close()
fig([('thrust',m,0,0.2) for m in (0.05,0.1,0.2,0.3)],[f"GM/(Rc²) = {m}" for m in (0.05,0.1,0.2,0.3)],"direction (deg); 0.2c toward 0° with thrust, not expanding, mass at 180°; weak solid, GR dashed","gr_const.png")
fig([('thrust',m,0.1,0.2) for m in (0.05,0.1,0.2,0.3)],[f"GM/(Rc²) = {m}" for m in (0.05,0.1,0.2,0.3)],"direction (deg); 0.2c toward 0° with thrust, expanding at 0.1c, mass at 180°; weak solid, GR dashed","gr_const_expand.png")
fig([('coast',0.1,u,0.3) for u in (0,0.15,0.9999)],["not expanding","expanding 0.15c","expanding ≈c"],"direction (deg); start 0.3c, GM/(Rc²) = 0.1, slowed by mass, t = 3; weak solid, GR dashed","gr_decel.png")
fig([('coast',0.1,0,v) for v in (0.1,0.3,0.5)],[f"start {v}c" for v in (0.1,0.3,0.5)],"direction (deg); not expanding, GM/(Rc²) = 0.1, slowed by mass, t = 3; weak solid, GR dashed","gr_decel_speeds.png")
