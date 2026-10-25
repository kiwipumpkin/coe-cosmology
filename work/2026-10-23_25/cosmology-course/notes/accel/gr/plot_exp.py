import sys,pickle; sys.path.insert(0,'.'); from grcore import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=pickle.load(open('boost.pkl','rb')); x=np.degrees(th)
def P(sets,fn,xl):
    f,ax=plt.subplots(figsize=(14,8))
    for (to,v0),col in zip(sets,('C0','C3','C2')):
        r=norm(R[(to,0.1,0.9999,v0)][1][0]); d,q,_,mn=stats(r)
        ax.plot(x,r,color=col,lw=3,label=f"start {v0}c (viewed t = {to:g}): dipole {d:.4f}, quad {q:.4f}, coldest {mn:.4f}")
    ax.axhline(0,color='gray'); ax.set_xlim(0,360); ax.grid(alpha=.3); ax.set_ylabel("ΔT/T",fontsize=16); ax.tick_params(labelsize=14)
    ax.set_xlabel(xl,fontsize=15); ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.02),fontsize=14); plt.tight_layout(); plt.savefig(fn,dpi=120); plt.close()
P([(2.3,0.1),(2.3,0.3),(2.3,0.5)],'expanding_speeds_slow.png',"direction (deg); exact GR, shell expanding ≈c (round in its own frame), GM/(Rc²) = 0.1 at 180°")
P([(11.0,0.9),(106.0,0.99)],'expanding_speeds_fast.png',"direction (deg); exact GR, shell expanding ≈c (round in its own frame), GM/(Rc²) = 0.1 at 180°")
