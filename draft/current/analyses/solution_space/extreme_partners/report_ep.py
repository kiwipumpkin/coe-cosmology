import json,glob,os,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from ep_core import *
plt.rcParams.update({"font.size":14,"axes.titlesize":15,"legend.fontsize":12})
ORDER=[("over_1p","1"),("over_2p","2"),("over_3p","3"),("under_1p","1"),("under_2p","2"),("under_3p","3"),("free_1p","1"),("free_2p","2"),("free_3p","3"),
       ("over_2p+4C","2\n+4M"),("over_3p+4C","3\n+4M"),("over_2p+4C4V","2\n+4M\n+4V"),("under_2p+4C","2\n+4M"),("under_3p+4C","3\n+4M"),("under_2p+4C4V","2\n+4M\n+4V"),
       ("free_2p+4C","2\n+4M"),("free_3p+4C","3\n+4M"),("free_2p+4C4V","2\n+4M\n+4V")]
COL={"over":"#2471a3","under":"#c0392b","free":"#7f8c8d"}
LAB={"over":"we're overdense, partners are giant voids","under":"we're underdense, partners are giant masses","free":"partners free (either sign)"}
best={}; stats={}
for nm,_ in ORDER:
    d=[]
    for f in (f"raw/{nm}.json",f"raw_b/{nm}.json"):
        if os.path.exists(f): d+=json.load(open(f))
    if not d: continue
    ok=[r for r in d if r["viol"]<=1e-3 and abs(r["net"])<=NET_TOL]
    if ok: best[nm]=min(ok,key=lambda r:r["left"]); stats[nm]=(len(ok),len(d),float(np.median([r["left"] for r in ok])))
REF=[(237.3,"start 237×"),(35.6,"old: 4 masses 36×"),(4.46,"old: 4 masses\n+ 4 voids 4.5×")]
fig,ax=plt.subplots(figsize=(14,7.5)); xs=[]; xl=[]
for i,(nm,lab) in enumerate(ORDER):
    if nm not in best: continue
    sg=nm.split("_")[0]; r=best[nm]
    ax.scatter(i,r["left"],s=110,color=COL[sg],zorder=3,label=LAB[sg] if sg not in [h.get_label() for h in ax.collections] else None)
    ax.text(i,r["left"]*1.25,f"{r['left']:.0f}×" if r["left"]>=10 else f"{r['left']:.1f}×",ha="center",fontsize=11)
    xs.append(i); xl.append(lab)
seen=set()
for c in ax.collections:
    pass
for v,t in REF:
    ax.axhline(v,color="0.6",ls="--",lw=1); ax.text(len(ORDER)-0.3,v*1.04,t,va="center",fontsize=11,color="0.35")
ax.axhline(1,color="green",lw=2); ax.text(len(ORDER)-0.3,1.04,"CMB flat target 1×",va="center",fontsize=11,color="green")
ax.axvline(8.5,color="0.85",lw=1)
ax.text(4,0.55,"giant partners only",ha="center",fontsize=12,color="0.3"); ax.text(13,0.55,"giant partners + ordinary masses (M) / voids (V)",ha="center",fontsize=12,color="0.3")
ax.set_yscale("log"); ax.set_ylim(0.4,600); ax.set_xlim(-0.6,len(ORDER)+2.6)
ax.set_yticks([1,2,5,10,20,50,100,200]); ax.set_yticklabels(["1","2","5","10","20","50","100","200"])
ax.set_xticks(xs); ax.set_xticklabels(xl,fontsize=11); ax.set_xlabel("number of partners (+ M = ordinary masses, V = ordinary voids)")
ax.set_ylabel("best leftover (× floor, log scale)"); ax.set_title("Our sphere over/underdense + giant partners: best leftover")
h,l=ax.get_legend_handles_labels(); u=dict(zip(l,h)); ax.legend(u.values(),u.keys(),loc="upper center",bbox_to_anchor=(0.5,-0.2),ncol=3,frameon=False)
fig.tight_layout(); fig.savefig("a_best_by_setup_logy.png",dpi=110); plt.close(fig)
def lb(X): d=np.linalg.norm(X,axis=1); return np.degrees(np.arctan2(X[:,1],X[:,0]))%360,np.degrees(np.arcsin(X[:,2]/d)),d
def draw_map(nmb,fname,extra=''):
    b=best[nmb]
    # map of best overall
    Xp=np.array(b["Xp"]); Mp=np.array(b["Mp"]); Rp=np.array(b["Rp"]); Xo=np.array(b["Xo"]).reshape(-1,3); Mo=np.array(b["Mo"]); Rvo=np.array(b["Rvo"])
    fig,(a1,a2)=plt.subplots(1,2,figsize=(16,7.5)); t=np.linspace(0,2*np.pi,300)
    a2.fill(14.3*np.cos(t),14.3*np.sin(t),color="#2471a3" if b["delta_us"]<0 else "#c0392b",alpha=0.10)
    a2.plot(14.3*np.cos(t),14.3*np.sin(t),"k",lw=1)
    a2.plot(C_G[0]+A_IN_G*np.cos(t),C_G[1]+A_IN_G*np.sin(t),color="0.5",lw=1); a2.plot(C_G[0]+72.5*np.cos(t),C_G[1]+72.5*np.sin(t),color="0.5",lw=1,ls=":")
    ax_l=np.radians(AXIS); a2.plot([-62*np.cos(ax_l),62*np.cos(ax_l)],[-62*np.sin(ax_l),62*np.sin(ax_l)],color="0.7",ls="--",lw=1)
    l,bb,d=lb(Xp)
    for i in range(len(Mp)):
        c="#2471a3" if Mp[i]<0 else "#c0392b"
        a2.add_patch(plt.Circle(Xp[i,:2],Rp[i],facecolor=c,alpha=0.25,edgecolor=c,lw=2.5))
        a2.text(Xp[i,0],Xp[i,1],f"P{i+1}",ha="center",va="center",fontsize=12,fontweight="bold")
        a1.scatter(l[i],bb[i],s=200+60*Rp[i],marker="o",facecolors=c,alpha=0.35,edgecolors=c,lw=2.5); a1.text(l[i],bb[i],f"P{i+1}",ha="center",va="center",fontsize=11,fontweight="bold")
    if len(Mo):
        l,bb,d=lb(Xo)
        for i in range(len(Mo)):
            if Rvo[i]>0:
                a1.scatter(l[i],bb[i],s=60+240*Rvo[i]/4,facecolors="none",edgecolors="#2471a3",lw=2); a2.add_patch(plt.Circle(Xo[i,:2],Rvo[i],color="#2471a3",alpha=0.35))
            else:
                s=50+200*Mo[i]/2e52; a1.scatter(l[i],bb[i],s=s,c="#c0392b"); a2.scatter(Xo[i,0],Xo[i,1],s=s,c="#c0392b",zorder=3)
    a2.plot(0,0,"k+",ms=14); a2.set_aspect("equal"); a2.set_xlim(-62,62); a2.set_ylim(-62,62); a2.set_xlabel("Gpc")
    a2.set_title("Top view (us = +, shaded circle = our 14.3 Gpc sphere,\ngrey arc = CoE interior edge, dotted = shell outer edge)",fontsize=12)
    a1.set_xlim(360,0); a1.set_ylim(-90,90); a1.set_xticks(range(0,361,45)); a1.axhline(0,color="0.8",lw=1)
    a1.axvline(AXIS,color="0.7",ls="--",lw=1); a1.axvline(AXIS+180,color="0.7",ls="--",lw=1)
    a1.axvspan(270,300,color="0.9",zorder=-1); a1.text(285,82,"Migkas",ha="center",fontsize=11)
    a1.set_xlabel("galactic longitude ℓ (deg)"); a1.set_ylabel("latitude b (deg)")
    a1.set_title("Sky positions (P = giant partner; red = extra mass, blue = void;\nsmall red dots = ordinary masses, rings = ordinary voids)",fontsize=12)
    fig.suptitle(f"Best extreme-partner layout ({nmb}): our sphere δ = {b['delta_us']:+.2f}, leftover {b['left']:.1f}×"+extra); fig.tight_layout()
    fig.savefig(fname,dpi=100); plt.close(fig)
nmb=min(best,key=lambda n:best[n]["left"]); draw_map(nmb,"d_best_partner_map.png")
giant=[n for n in best if max(best[n]["Rp"])>=8]
if giant:
    ng=min(giant,key=lambda n:best[n]["left"]); draw_map(ng,"d_best_giant_partner_map.png"," (best layout using a giant partner)")
out={}
for nm,_ in ORDER:
    if nm not in best: continue
    r=best[nm]; Xp=np.array(r["Xp"]); dc=np.linalg.norm(Xp-C_G,axis=1); d=np.linalg.norm(Xp,axis=1); Rp=np.array(r["Rp"]); dp=np.array(r["dp"])
    l,bb,_=lb(Xp); Mo=np.array(r["Mo"]); Rvo=np.array(r["Rvo"])
    lo,hi=SIGN[r["sign"]][1]; dlo,dhi=SIGN[r["sign"]][0]
    bounds=dict(delta_us_at_bound=bool(abs(r["delta_us"]-dlo)<1e-3 or abs(r["delta_us"]-dhi)<1e-3),
        partner_delta_at_bound=int(((np.abs(dp-lo)<1e-3)|(np.abs(dp-hi)<1e-3)).sum()),partner_R_at_15=int((Rp>=14.9).sum()),partner_R_at_1=int((Rp<=1.02).sum()),
        partner_touch_sphere=int((d-Rp<=14.35).sum()),partner_touch_CoE=int((dc-Rp<=A_IN_G+0.05).sum()),partner_touch_shell_out=int((dc+Rp>=72.45).sum()),
        ord_mass_cap=int(((Rvo==0)&(Mo>=0.98*MHI)).sum()) if len(Mo) else 0,ord_void_cap=int((Rvo>=3.98).sum()) if len(Mo) else 0)
    out[nm]=dict(left=r["left"],ok_runs=stats[nm][0],runs=stats[nm][1],median=stats[nm][2],delta_us=r["delta_us"],M_sphere=r["M_sphere"],
        partners=[dict(l=float(l[i]),b=float(bb[i]),d_us=float(d[i]),d_CoE=float(dc[i]),R=float(Rp[i]),delta=float(dp[i]),M=float(r["Mp"][i])) for i in range(len(Rp))],
        M_partners=float(np.sum(r["Mp"])),M_ordinary=float(Mo.sum()) if len(Mo) else 0.0,M_ord_pos=float(Mo[Mo>0].sum()) if len(Mo) else 0.0,M_ord_neg=float(Mo[Mo<0].sum()) if len(Mo) else 0.0,net=r["net"],bounds=bounds,
        top_modes=(np.argsort(r["amps"])[::-1][:3]+1).tolist())
json.dump(out,open("ep_results.json","w"),indent=1)
for nm,v in out.items():
    print(f"{nm:14s} left={v['left']:7.2f} (ok {v['ok_runs']}/{v['runs']}, med {v['median']:.1f}) d_us={v['delta_us']:+.3f} Msph={v['M_sphere']:.2e} Mpart={v['M_partners']:.2e} Mord={v['M_ordinary']:.2e} net={v['net']:.1e} modes={v['top_modes']}")
    for p in v["partners"]: print("     P l=%.0f b=%+.0f d=%.1f dCoE=%.1f R=%.1f delta=%+.2f M=%.2e"%(p["l"],p["b"],p["d_us"],p["d_CoE"],p["R"],p["delta"],p["M"]))
    print("     bounds",v["bounds"])
print("BEST",nmb)
