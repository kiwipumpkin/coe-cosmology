import sys,pickle; sys.path.insert(0,'.'); from runall import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def job2(c):
    mu,u,v0=c; return c,weak_coast(mu,u,v0,2.0),gr_coast(mu,u,v0,2.0)
if __name__=='__main__':
    C=[(mu,u,v0) for v0 in (0.1,0.3,0.5) for u in (0,0.9999) for mu in (0.0,0.1)]
    with mp.Pool(8) as P: R=dict((c,(w,g)) for c,w,g in P.map(job2,C))
    pickle.dump(R,open('speeds.pkl','wb'))
    for c,(w,g) in sorted(R.items()):
        gr=norm(g[0]); print(c,"weak d %.4f q %.4f min %.4f"%(*stats(w)[:2],stats(w)[3]),"| GR d %.4f q %.4f min %.4f"%(*stats(gr)[:2],stats(gr)[3]),"failed",g[1],"horizon",g[2],flush=True)
    fig,ax=plt.subplots(figsize=(14,8.5)); x=np.degrees(th)
    for v0,col in ((0.1,'C0'),(0.3,'C3'),(0.5,'C2')):
        for u,st,l in ((0,'-','not expanding'),(0.9999,'--','expanding ≈c')):
            r=norm(R[(0.1,u,v0)][1][0]); d,q,_,_=stats(r)
            ax.plot(x,r,st,color=col,lw=3,label=f"start {v0}c, {l}: dipole {d:.4f}, quad {q:.4f}")
    ax.axhline(0,color='gray'); ax.set_xlim(0,360); ax.grid(alpha=.3); ax.set_ylabel("ΔT/T",fontsize=16); ax.tick_params(labelsize=14)
    ax.set_xlabel("direction (deg); exact GR; motion toward 0°, GM/(Rc²) = 0.1 at 180°, viewed at t = 2",fontsize=15)
    ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.02),ncol=2,fontsize=13); plt.tight_layout(); plt.savefig('decel_overlay_speeds_gr.png',dpi=120)
