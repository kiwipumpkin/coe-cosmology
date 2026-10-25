import sys,pickle; sys.path.insert(0,'.'); from grcore import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
S=pickle.load(open('sweep_static.pkl','rb')); S[0.0]=pickle.load(open('sweep_v0.pkl','rb'))
V=[0.0,0.1,0.3,0.5,0.9,0.99]; cm=plt.cm.viridis; x=np.degrees(th)
f,ax=plt.subplots(figsize=(14,8.5))
for i,v in enumerate(V):
    to,res=S[v]; r=norm(res[(0.1,0.9999)]); d,q,_,mn=stats(r)
    ax.plot(x,r,color=cm(i/(len(V)-1)*0.9),lw=3,label=f"{v:g}c (t={to:.1f}): dip {d:.4f}, quad {q:.4f}")
    print(v,'%.1f'%to,'%.4f %.4f %.4f'%(d,q,mn))
ax.axhline(0,color='gray'); ax.set_xlim(0,360); ax.grid(alpha=.3); ax.set_ylabel("ΔT/T",fontsize=16); ax.tick_params(labelsize=14)
ax.set_xlabel("direction (deg); exact GR, shell expanding ≈c and moving toward 0° at the start speed shown, GM/(Rc²) = 0.1 at 180°; each viewed at its earliest full-sky time",fontsize=12)
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.02),ncol=3,fontsize=12,title="start speed (view time)"); plt.tight_layout(); plt.savefig('combined_levels.png',dpi=120)
