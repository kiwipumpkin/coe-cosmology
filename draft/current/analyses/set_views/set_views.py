"""Views of the latest set (no new fits): inside-only, outside-only (no wall, smooth limit), all-open (smooth limit)."""
import os,sys,json,glob,numpy as np
os.environ["NCW"]="smooth"; os.environ["AO"]="smooth"
CS="/workspace/clump_scan"; SS=f"{CS}/solution_space"
for p in (SS,f"{CS}/inside_visible",f"{SS}/no_core_wall",f"{SS}/all_open"): sys.path.insert(0,p)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from ss_core import R,FLOOR,TH,rms,amps,C_G,A_IN_G,RV
_c=os.getcwd(); from inside_core import InModel; from ncw_core import NCW; from ao_core import AO; os.chdir(_c)
plt.rcParams.update({"font.size":13,"axes.titlesize":13})
KEY=[4,6,8,10,11,12]
def best_by_N(runs):
    b={}
    for r in runs:
        if r["viol"]>1e-3: continue
        N=len(r["kinds"])
        if N in KEY and (N not in b or r["left"]<b[N]["left"]): b[N]=r
    return b
ins=json.load(open(f"{CS}/inside_visible/results_cap2e52.json")); B_in=best_by_N([c["best_run"] for c in ins.values()])
B_out=best_by_N([r for f in glob.glob(f"{SS}/no_core_wall/raw_smooth/*.json") for r in json.load(open(f))])
B_ao=best_by_N([r for f in glob.glob(f"{SS}/all_open/raw_smooth/*.json") for r in json.load(open(f))])
SC={"inside":("inside our sphere only",B_in,lambda k:InModel(k,2e52,"local"),15),
    "outside":("outside only (no wall, smooth-density limit)",B_out,NCW,45),
    "all_open":("all open (smooth-density limit)",B_ao,AO,45)}
COL={"inside":"#c0392b","outside":"#7f8c8d","all_open":"#1f77b4"}
res={}
for s,(lab,B,mk,lim) in SC.items():
    res[s]={}
    for N,r in sorted(B.items()):
        m=mk(r["kinds"]); p=np.array(r["p"]); w=m.pattern(p); left=rms(R-w)/FLOOR
        assert abs(left-r["left"])<1e-3*max(1,r["left"]), (s,N,left,r["left"])
        v=R-w; th=TH; rec=np.zeros(360)
        for mm in range(0,13): rec+=(np.mean(v) if mm==0 else 2*np.mean(v*np.cos(mm*th))*np.cos(mm*th)+2*np.mean(v*np.sin(mm*th))*np.sin(mm*th))
        res[s][N]=dict(r=r,resid=v,amps=amps(v)/FLOOR,left=left,hi=float(np.sqrt(2)*rms(v-rec)/FLOOR))
# 1. fit vs N
fig,ax=plt.subplots(figsize=(10,6.5))
for s,mkr,ls in (("inside","o","-"),("outside","s","--"),("all_open","D","-")):
    N=sorted(res[s]); ax.plot(N,[res[s][n]["left"] for n in N],ls,marker=mkr,color=COL[s],lw=2.5,ms=8,label=SC[s][0])
ax.axhline(1,color="green",lw=2); ax.text(12.3,1.1,"CMB flat\ntarget 1×",color="green",fontsize=12,va="bottom")
ax.axhline(237.3,color="0.5",ls=":",lw=1.3); ax.text(12.3,237.3*1.1,"start 237×",color="0.4",fontsize=12,va="bottom")
ax.set_yscale("log"); ax.set_ylim(0.01,600); ax.set_xlim(3.5,13.4); ax.set_xticks(KEY)
ax.set_yticks([0.01,0.1,1,10,100]); ax.set_yticklabels(["0.01","0.1","1","10","100"])
ax.set_xlabel("number of objects N"); ax.set_ylabel("best leftover (× 10⁻⁵ floor, log scale)")
ax.set_title("Fit quality vs number of objects"); ax.legend(loc="lower left",frameon=False,fontsize=12)
fig.tight_layout(); fig.savefig("fit_vs_N.png",dpi=110); plt.close(fig)
# 2. top views
t=np.linspace(0,2*np.pi,400)
def mixstr(k): c=k.count("C"); v=k.count("V"); return f"{c} masses" + (f" + {v} void" + ("s" if v>1 else "") if v else "")
for s,(lab,B,mk,lim) in SC.items():
    fig,axs=plt.subplots(2,3,figsize=(16,11))
    for ax,N in zip(axs.ravel(),KEY):
        ax.set_aspect("equal"); ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
        ax.plot(RV*np.cos(t),RV*np.sin(t),"k",lw=1.2)
        ax.plot(C_G[0]+A_IN_G*np.cos(t),C_G[1]+A_IN_G*np.sin(t),color="0.7",ls="--",lw=1)
        ax.plot(0,0,"k+",ms=14,mew=2)
        if N not in res[s]: ax.set_title(f"N={N}: not run"); continue
        r=res[s][N]["r"]; k=np.array(list(r["kinds"])); X=np.array(r["X"]); M=np.array(r["M"]); Rv=np.array(r["Rv"])
        for i in range(len(k)):
            if k[i]=="C": ax.scatter(X[i,0],X[i,1],s=40+260*M[i]/2e52,c="#c0392b",zorder=3,edgecolors="k",lw=0.5)
            else: ax.add_patch(plt.Circle(X[i,:2],Rv[i],color="#2471a3",alpha=0.4))
        ax.set_title(f"N={N}: {mixstr(r['kinds'])}\nleftover {res[s][N]['left']:.2f}×")
        ax.tick_params(labelsize=11)
    for ax in axs[1]: ax.set_xlabel("Gpc")
    fig.suptitle(f"Top views, {lab}\n(us = +, black circle = our 14.3 Gpc sphere, grey dashed = old wall; red = masses sized by mass, blue = voids at true radius;\nobjects above/below the plane are projected onto it)",fontsize=14)
    fig.tight_layout(rect=(0,0,1,0.93)); fig.subplots_adjust(top=0.86); fig.savefig(f"topviews_{s}.png",dpi=95); plt.close(fig)
# 3a. modes
m_=np.arange(1,13); start_amps=amps(R)/FLOOR
for s,(lab,B,mk,lim) in SC.items():
    fig,axs=plt.subplots(2,3,figsize=(16,9.5),sharey=True)
    for ax,N in zip(axs.ravel(),KEY):
        if N not in res[s]: ax.set_title(f"N={N}: not run"); continue
        a=res[s][N]["amps"]; ax.bar(m_,a,color=COL[s]); hi=res[s][N]["hi"]; ax.bar(13.6,hi,color="0.55",width=1.2)
        ax.axhline(1,color="green",lw=2); ax.set_yscale("log"); ax.set_ylim(1e-3,300); ax.set_xticks(list(m_)+[13.6]); ax.set_xticklabels([str(i) for i in m_]+['13+'])
        ax.set_title(f"N={N}: leftover {res[s][N]['left']:.2f}×"); ax.tick_params(labelsize=11)
    for ax in axs[:,0]: ax.set_ylabel("amplitude (× floor)")
    for ax in axs[1]: ax.set_xlabel("mode m")
    fig.suptitle(f"What's left after the fit, per mode — {lab}\n(green line = CMB flat 1×; grey bar '13+' = all modes above 12 combined)",fontsize=15)
    fig.tight_layout(); fig.subplots_adjust(top=0.88); fig.savefig(f"residual_modes_{s}.png",dpi=95); plt.close(fig)
# 3b. ring
deg=np.degrees(TH)
for s,(lab,B,mk,lim) in SC.items():
    Ns=sorted(res[s]); bestN=min(Ns,key=lambda n:res[s][n]["left"]); pick=[n for n in (4,8) if n in Ns]+([bestN] if bestN not in (4,8) else [])
    fig,axs=plt.subplots(len(pick)+1,1,figsize=(12,3.0*(len(pick)+1)),sharex=True)
    axs[0].plot(deg,R/FLOOR,color="0.3",lw=2); axs[0].set_title("Starting residual (237×)",fontsize=13)
    for ax,N in zip(axs[1:],pick):
        y=res[s][N]["resid"]/FLOOR; ax.plot(deg,y,color=COL[s],lw=2)
        ax.axhline(1,color="green",lw=1.2,ls="--"); ax.axhline(-1,color="green",lw=1.2,ls="--")
        ax.set_title(f"After the fit, N={N} ({mixstr(res[s][N]['r']['kinds'])}): leftover {res[s][N]['left']:.2f}×   (green dashed = ±1)",fontsize=13)
    for ax in axs: ax.axhline(0,color="0.8",lw=1); ax.set_ylabel("ΔT/T (10⁻⁵)",fontsize=12)
    axs[-1].set_xlabel("angle around the ring (galactic longitude ℓ, degrees)"); axs[-1].set_xticks(range(0,361,45)); axs[-1].set_xlim(0,360)
    fig.suptitle(f"Leftover pattern around the CMB ring — {lab}\n(each panel has its own y scale)",fontsize=14)
    fig.tight_layout(); fig.subplots_adjust(top=0.91); fig.savefig(f"residual_ring_{s}.png",dpi=95); plt.close(fig)
# notes
for s in res:
    for N in sorted(res[s]):
        a=res[s][N]["amps"]; top=np.argsort(a)[::-1][:3]
        print(s,N,res[s][N]["r"]["kinds"],"left %.3f"%res[s][N]["left"],"top m",[(int(i+1),round(float(a[i]),3)) for i in top],"m1-6 max %.3f"%a[:6].max(),"13+ %.3f"%res[s][N]["hi"])
print("start amps",np.round(start_amps,1))
