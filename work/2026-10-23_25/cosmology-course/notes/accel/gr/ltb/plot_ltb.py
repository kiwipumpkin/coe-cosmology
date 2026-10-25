import sys,pickle; sys.path.insert(0,'..'); import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from grcore import stats, norm
d=pickle.load(open('ltb_0.0_0.1_0.3_0.5_0.9_0.99.pkl','rb')); V=sorted(k for k in d if d[k][1])
th=np.arange(0,361); cm=plt.cm.viridis
for core,fn,ttl in [(0.0946,'combined_levels_ltb.png','with central core (0.0946 of critical mass)'),(0.0,'combined_levels_ltb_nocore.png','no core: uniform expanding space (baseline)')]:
    fig,ax=plt.subplots(figsize=(9,5.5))
    for i,v in enumerate(V):
        x=norm(d[v][1][core][0]); dp,q,_,_=stats(x)
        ax.plot(th,x,color=cm(i/max(1,len(V)-1)),lw=1.8,label=f"{v:g}c (t={d[v][0]:.1f}): dip {dp:.4f}, quad {q:.4f}")
    ax.axhline(0,color='k',lw=.5); ax.set_xlim(0,360); ax.set_xticks(range(0,361,45))
    ax.set_xlabel('LTB a(r,t): space expanding (3c recession) + motion through space, CoE at 180 deg'); ax.set_ylabel('ΔT/T')
    ax.set_title(ttl,fontsize=10,pad=58)
    ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=8,frameon=False)
    fig.tight_layout(); fig.savefig('../'+fn,dpi=150); plt.close(fig)
for v in V:
    a=norm(d[v][1][0.0946][0]); b=norm(d[v][1][0.0][0]); print(v,d[v][0],["%.4f"%s for s in stats(a)[:2]]+["%.4f"%a.min()],["%.4f"%s for s in stats(b)[:2]]+["%.4f"%b.min()],["%.4f"%s for s in stats(a-b)[:2]]+["%.4f"%(a-b).min()])
