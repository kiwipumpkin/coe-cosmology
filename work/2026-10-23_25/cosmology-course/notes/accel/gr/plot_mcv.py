import sys,pickle; sys.path.insert(0,'.'); from grcore import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=pickle.load(open('mcv_2.0.pkl','rb')); x=np.degrees(th)
f,ax=plt.subplots(figsize=(14,8))
for u,st,l,c in ((0.0,'-','shell carried by the expansion','C0'),(0.9999,'--','shell also expanding ≈c through space','C3')):
    r=norm(R[(0.1,u,2.0)][0]); d,q,_,mn=stats(r); ax.plot(x,r,st,color=c,lw=3,label=f"{l}: dipole {d:.4f}, quad {q:.4f}, coldest {mn:.4f}")
ax.axhline(0,color='gray'); ax.set_xlim(0,360); ax.grid(alpha=.3); ax.set_ylabel("ΔT/T",fontsize=16); ax.tick_params(labelsize=14)
ax.set_xlabel("direction (deg); exact GR, point mass in expanding universe (McVittie), H R/c = 0.3, GM/(Rc²) = 0.1 at 180°, viewed 2R/c later",fontsize=13)
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.02),fontsize=14); plt.tight_layout(); plt.savefig('mcvittie_3c.png',dpi=120)
